//! One time step: advect, detect boundaries, subduct/collide, create crust at gaps,
//! evolve crust (erosion, hotspots, rift thinning), force balance, rifting, suturing, cleanup.
use crate::geom::*;
use crate::world::*;
use rand::Rng;
use rayon::prelude::*;
use std::collections::HashMap;

pub fn step(w: &mut World) {
    w.stats = Stats::default();
    let dt = w.p.dt;

    // 1. Rigid-plate advection (and accumulate each plate's finite rotation for export).
    for pc in w.parcels.iter_mut() {
        if pc.alive {
            let om = w.plates[pc.plate as usize].omega;
            pc.pos = rotate(pc.pos, om, dt);
        }
    }
    for (i, pl) in w.plates.iter().enumerate() {
        if pl.alive { w.rot[i] = mat_mul(rot_matrix(pl.omega, dt), w.rot[i]); }
    }
    w.rebuild_hash();
    // 2. Boundary detection.
    w.binfo = detect(w);
    // 3. Subduction, collision, accretion.
    interact(w);
    // 4. New crust at divergent gaps (Rule III: ridges are passive - crust breaks because it is pulled).
    fill_gaps(w);
    // 5. Crust evolution.
    relax(w);
    // 6. Force balance -> plate velocities.
    crate::forces::update_omegas(w);
    // 7. Continental rifting.
    rifting(w);
    // 7b. Arc detachment / slab rollback / back-arc basins.
    backarc(w);
    // 8. Suturing and removal of consumed plates.
    merge_and_cleanup(w);
    plate_stats(w);
    w.rebuild_hash();
    w.t += dt;
}

fn detect(w: &World) -> Vec<Option<BInfo>> {
    // Wide enough to see a passive margin (200 km); force and consumption rules apply their own radii.
    let r = w.reach(1.5, 200.0);
    (0..w.parcels.len())
        .into_par_iter()
        .map(|i| {
            let pi = &w.parcels[i];
            if !pi.alive { return None; }
            let mut best: Option<(u32, f64)> = None;
            w.hash.query(pi.pos, r, |j| {
                let pj = &w.parcels[j as usize];
                if pj.alive && pj.plate != pi.plate {
                    let d = dist(pi.pos, pj.pos);
                    if d < r && best.map_or(true, |(_, bd)| d < bd) { best = Some((j, d)); }
                }
            });
            best.map(|(j, d)| {
                let pj = &w.parcels[j as usize];
                let n = tangent_toward(pi.pos, pj.pos);
                let vi = surface_velocity(w.plates[pi.plate as usize].omega, pi.pos);
                let vj = surface_velocity(w.plates[pj.plate as usize].omega, pj.pos);
                BInfo { other: pj.plate, oidx: j, n, conv: dot(sub(vi, vj), n), dist: d }
            })
        })
        .collect()
}

fn interact(w: &mut World) {
    let t = w.t;
    let s = w.s;
    let dt = w.p.dt;
    let r_col = 0.8 * s;
    let r_deep = 0.45 * s;
    let r_contact = 1.5 * s;
    let r_arc = w.reach(1.5, 170.0);
    let r_fl = w.reach(2.5, 280.0);
    let r_virus = w.km(w.p.virus_km);
    let mut kills: Vec<usize> = vec![];
    let mut arcs: Vec<usize> = vec![];
    let mut thick_add: Vec<(usize, f64)> = vec![];
    let mut transfers: Vec<(usize, u32)> = vec![];
    let mut sutures: Vec<usize> = vec![];
    let mut absorbed: Vec<usize> = vec![];
    let mut obducted: Vec<usize> = vec![];
    let mut comp_step: HashMap<(u32, u32), (f64, u32)> = HashMap::new();

    // Active trench points (lower-plate contact parcels of established subduction zones): a new
    // convergent contact next to one of these starts subducting at once - Scotese's "subduction virus".
    let mut trench_hash = SpatialHash::new(r_virus.max(1.5 * s));
    {
        let pts = w.parcels.iter().enumerate().filter_map(|(i, pc)| {
            let b = w.binfo.get(i).copied().flatten()?;
            if !pc.alive || b.conv <= CONV_EPS || b.dist >= r_contact { return None; }
            let key = (pc.plate.min(b.other), pc.plate.max(b.other));
            match w.polarity.get(&key) { Some(&sp) if sp == pc.plate => Some((i as u32, pc.pos)), _ => None }
        });
        trench_hash.build(pts);
    }

    for i in 0..w.parcels.len() {
        let b = match w.binfo[i] { Some(b) => b, None => continue };
        let pi = w.parcels[i];
        if !pi.alive || b.dist >= r_contact { continue; }
        let j = b.oidx as usize;
        let pj = w.parcels[j];
        if !pj.alive { continue; }
        let (a, bp) = (pi.plate, pj.plate);
        if a == bp { continue; }
        let key = (a.min(bp), a.max(bp));
        // Consume on any convergence once overlapping, or on deep overlap even at near-zero convergence
        // (oblique / transform contacts must not interpenetrate).
        let closing = (b.conv > 0.0 && b.dist < r_col) || (b.conv > -CONV_EPS && b.dist < r_deep);

        let sub = match w.polarity.get(&key) {
            Some(&sp) => sp,
            None => {
                if b.conv <= CONV_EPS { continue; }
                // Convergence without a subduction zone: the contact takes up compression (resisted in the
                // force balance) until subduction initiates. Rule IV: subduction is hard to start.
                let e = comp_step.entry(key).or_insert((0.0, 0));
                e.0 += b.conv;
                e.1 += 1;
                let lower = match (pi.kind, pj.kind) {
                    (Kind::Oceanic, Kind::Continental) => a,
                    (Kind::Continental, Kind::Oceanic) => bp,
                    (Kind::Oceanic, Kind::Oceanic) => if pi.birth <= pj.birth { a } else { bp },
                    (Kind::Continental, Kind::Continental) => {
                        if w.plates[a as usize].n <= w.plates[bp as usize].n { a } else { bp }
                    }
                };
                let lower_p = if lower == a { &pi } else { &pj };
                let shortened = w.pair_compress.get(&key).copied().unwrap_or(0.0);
                let old_enough = lower_p.kind == Kind::Oceanic && t - lower_p.birth >= w.p.init_age;
                let cc = pi.kind == Kind::Continental && pj.kind == Kind::Continental;
                let mut virus = false;
                {
                    let parcels = &w.parcels;
                    trench_hash.query(pi.pos, r_virus, |k| { if !virus && dist(parcels[k as usize].pos, pi.pos) < r_virus { virus = true; } });
                }
                let start = cc || virus || b.dist < r_deep || (old_enough && shortened >= 20.0) || shortened >= w.p.init_short;
                if !start { continue; }
                w.polarity.insert(key, lower);
                w.pair_compress.remove(&key);
                w.stats.initiations += 1;
                lower
            }
        };
        if sub != a { continue; }

        // Arc volcanism on the upper plate above every active trench segment (continuous along the arc,
        // not only where a parcel happens to be consumed this step).
        if pi.kind == Kind::Oceanic && b.conv > CONV_EPS {
            let hash = &w.hash;
            let parcels = &w.parcels;
            hash.query(pj.pos, r_arc, |k| {
                let pk = &parcels[k as usize];
                if pk.alive && pk.plate == bp && dist(pk.pos, pj.pos) < r_arc { arcs.push(k as usize); }
            });
        }
        if !closing { continue; }

        if pi.kind == Kind::Oceanic {
            // Oceanic lithosphere of the lower plate is consumed.
            kills.push(i);
        } else if b.conv <= CONV_EPS {
            // Continental crust grinding along a transform / near-static contact: nothing is consumed.
            continue;
        } else if pj.kind == Kind::Continental && w.pair_absorbed.get(&key).copied().unwrap_or(0) >= w.rows_lock {
            // The collision belt has locked up after enough shortening. The whole connected continental
            // block of the lower plate is accreted to the upper plate in one event and the boundary
            // re-forms behind it (Rule IV: the trench jumps over the colliding terrane).
            let block = connected_continent(w, i, a);
            for &k in &block { w.parcels[k].plate = bp; }
            w.stats.accreted += block.len();
            sutures.push(i);
            sutures.push(j);
            w.pair_absorbed.remove(&key);
        } else if pj.kind == Kind::Continental {
            // Continent-continent collision: the arriving parcel is absorbed (shortening). Its crustal
            // volume thickens the surrounding belt on both plates (volume-conserving); the contact is a suture.
            let hash = &w.hash;
            let parcels = &w.parcels;
            let mut near: Vec<(usize, f64)> = vec![];
            let (r_belt, w_belt) = (w.reach(3.0, 350.0), w.reach(1.4, 160.0));
            hash.query(pi.pos, r_belt, |k| {
                let pk = &parcels[k as usize];
                if k as usize != i && pk.alive && pk.kind == Kind::Continental {
                    let d = dist(pk.pos, pi.pos);
                    if d < r_belt { near.push((k as usize, (-(d / w_belt).powi(2)).exp())); }
                }
            });
            if near.is_empty() { continue; }
            let wsum: f64 = near.iter().map(|x| x.1).sum::<f64>().max(1e-9);
            for (k, wk) in &near { thick_add.push((*k, pi.thick * wk / wsum)); sutures.push(*k); }
            absorbed.push(i);
            *w.pair_absorbed.entry(key).or_insert(0) += 1;
        } else if pi.thick < 33.0 && !attached_to_continent(w, i) {
            // Thin, isolated continental slivers (old arcs, stretched shelf) go down with the slab at
            // intra-oceanic trenches - island arcs do not ride across oceans (Rule VIII).
            kills.push(i);
            w.stats.cont_lost += 1;
        } else {
            // Continent arriving at an intra-oceanic trench: it docks onto the upper plate (terrane accretion);
            // the oceanic crust it overrides is obducted (removed). Subduction then resumes behind it (Rule IV).
            let hash = &w.hash;
            let parcels = &w.parcels;
            hash.query(pi.pos, 0.8 * s, |k| {
                let pk = &parcels[k as usize];
                if pk.alive && pk.plate == bp && pk.kind == Kind::Oceanic && dist(pk.pos, pi.pos) < 0.8 * s { obducted.push(k as usize); }
            });
            transfers.push((i, bp));
            sutures.push(i);
        }
    }
    for (key, (sum, n)) in comp_step {
        if n > 0 { *w.pair_compress.entry(key).or_insert(0.0) += sum / n as f64 * dt; }
    }

    let kinfo: Vec<(V3, u32)> = kills.iter().map(|&i| (w.parcels[i].pos, w.parcels[i].plate)).collect();
    for &i in &kills { w.parcels[i].alive = false; }
    // Flexural deepening on the lower plate: strongest at the trench, fading to the outer rise.
    let mut trench_marks: Vec<(usize, f64)> = vec![];
    {
        let hash = &w.hash;
        let parcels = &w.parcels;
        for &(pos, pl) in &kinfo {
            hash.query(pos, r_fl, |k| {
                let pk = &parcels[k as usize];
                if pk.alive && pk.plate == pl {
                    let d = dist(pk.pos, pos);
                    if d < r_fl { trench_marks.push((k as usize, (1.0 - d / r_fl).powi(2))); }
                }
            });
        }
    }
    for (k, wgt) in trench_marks {
        let pc = &mut w.parcels[k];
        let prev = if t - pc.trench_t < 5.0 { pc.trench_w } else { 0.0 };
        pc.trench_w = prev.max(wgt);
        pc.trench_t = t;
    }
    for &i in &absorbed { w.parcels[i].alive = false; }
    for &i in &obducted { w.parcels[i].alive = false; }
    for k in arcs { w.parcels[k].arc_t = t; }
    for (k, v) in thick_add { w.parcels[k].thick += v; }
    for &(i, pl) in &transfers { w.parcels[i].plate = pl; }
    for i in sutures { w.parcels[i].suture_t = t; }
    w.stats.subducted = kills.len() + obducted.len();
    w.stats.accreted += transfers.len();
    w.stats.absorbed = absorbed.len();
}

/// Flood-fill the connected continental block of plate `plate` that contains parcel `start`
/// (links: continental parcels of the same plate within 1.5 spacings).
fn connected_continent(w: &World, start: usize, plate: u32) -> Vec<usize> {
    let s = w.s;
    let mut seen = std::collections::HashSet::new();
    let mut stack = vec![start];
    let mut out = vec![];
    seen.insert(start);
    while let Some(i) = stack.pop() {
        out.push(i);
        let pi = w.parcels[i].pos;
        let mut next = vec![];
        w.hash.query(pi, 1.5 * s, |k| {
            let pk = &w.parcels[k as usize];
            if pk.alive && pk.plate == plate && pk.kind == Kind::Continental && dist(pk.pos, pi) < 1.5 * s { next.push(k as usize); }
        });
        for k in next { if seen.insert(k) { stack.push(k); } }
    }
    out
}

/// True if parcel `i` has a neighbouring continental parcel of normal thickness within 1.5 spacings.
fn attached_to_continent(w: &World, i: usize) -> bool {
    let pi = &w.parcels[i];
    let s = w.s;
    let mut found = false;
    w.hash.query(pi.pos, 1.5 * s, |k| {
        if found || k as usize == i { return; }
        let pk = &w.parcels[k as usize];
        if pk.alive && pk.kind == Kind::Continental && pk.thick >= 33.0 && dist(pk.pos, pi.pos) < 1.5 * s { found = true; }
    });
    found
}

fn fill_gaps(w: &mut World) {
    let s = w.s;
    let r_gap = 0.9 * s;
    let r_search = 2.0 * s;
    let mut cands: Vec<(usize, f64, u32)> = w
        .grid
        .par_iter()
        .enumerate()
        .filter_map(|(gi, &g)| {
            let mut best: Option<(u32, f64)> = None;
            w.hash.query(g, r_search, |j| {
                let pj = &w.parcels[j as usize];
                if pj.alive {
                    let d = dist(pj.pos, g);
                    if d < r_search && best.map_or(true, |(_, bd)| d < bd) { best = Some((j, d)); }
                }
            });
            match best {
                Some((j, d)) if d > r_gap => {
                    // Never create crust inside a closing (convergent) boundary.
                    if let Some(b) = w.binfo[j as usize] {
                        if b.conv > CONV_EPS && b.dist < 1.5 * s { return None; }
                    }
                    // Ownership by weighted majority of nearby parcels (coherent patches, not interleaving).
                    let mut votes: Vec<(u32, f64)> = Vec::with_capacity(4);
                    w.hash.query(g, r_search, |k| {
                        let pk = &w.parcels[k as usize];
                        if pk.alive {
                            let dk = dist(pk.pos, g);
                            if dk < r_search {
                                let wgt = (-(dk / s).powi(2)).exp();
                                match votes.iter_mut().find(|v| v.0 == pk.plate) { Some(v) => v.1 += wgt, None => votes.push((pk.plate, wgt)) }
                            }
                        }
                    });
                    let owner = votes.iter().max_by(|a, b| a.1.partial_cmp(&b.1).unwrap()).map(|v| v.0).unwrap_or(w.parcels[j as usize].plate);
                    Some((gi, d, owner))
                }
                _ => None,
            }
        })
        .collect();
    cands.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    let mut placed: Vec<V3> = vec![];
    let t = w.t;
    for (gi, _, pl) in cands {
        let g = w.grid[gi];
        if placed.iter().any(|&q| dist(q, g) < 0.85 * s) { continue; }
        placed.push(g);
        w.parcels.push(Parcel {
            pos: g, plate: pl, kind: Kind::Oceanic, birth: t, thick: 7.0, volc: 0.0,
            trench_t: NEVER, suture_t: NEVER, hot_t: NEVER, arc_t: NEVER, rift_t: NEVER, trench_w: 0.0, alive: true,
        });
    }
    w.stats.created = placed.len();
}

fn relax(w: &mut World) {
    let dt = w.p.dt;
    let t = w.t;
    let s = w.s;
    let r_hot = w.reach(0.8, 90.0);
    let p = w.p.clone();
    for (h, d) in w.hotspots.iter_mut().zip(w.hot_drift.iter()) { *h = normalize(add(*h, scale(*d, dt))); }
    let hs = w.hotspots.clone();
    // Erosion as diffusion of crustal thickness between neighbouring continental parcels
    // (volume-conserving: mountains shed into their forelands and shelves), plus a slow loss term.
    // Scaled by (113 km / spacing)^2 so the physical diffusivity does not depend on resolution.
    let kappa = (dt / p.erosion_tau * (113.0 / (s * R_KM)).powi(2)).min(0.45);
    let r_d = 1.6 * s;
    let r_margin = w.reach(2.2, 250.0);
    // Cold cratonic lithosphere (original crust, not orogenically thickened) barely erodes; young
    // orogens and rifted margins relax. This keeps shield relief while mountain belts decay.
    let mobility = |pc: &Parcel| -> f64 { if pc.birth < -2500.0 && pc.thick < 40.0 && pc.thick > 32.0 { 0.12 } else { 1.0 } };
    let deltas: Vec<f64> = (0..w.parcels.len()).into_par_iter().map(|i| {
        let pc = &w.parcels[i];
        if !pc.alive || pc.kind != Kind::Continental { return 0.0; }
        let (mut acc, mut wsum) = (0.0, 0.0);
        w.hash.query(pc.pos, r_d, |k| {
            let pk = &w.parcels[k as usize];
            if k as usize != i && pk.alive && pk.kind == Kind::Continental {
                let d = dist(pk.pos, pc.pos);
                if d < r_d { let wg = (-(d / s).powi(2)).exp(); acc += wg * (pk.thick - pc.thick); wsum += wg; }
            }
        });
        if wsum > 0.0 { kappa * acc / wsum.max(1.0) * mobility(pc) } else { 0.0 }
    }).collect();
    let binfo = &w.binfo;
    w.parcels.par_iter_mut().enumerate().for_each(|(i, pc)| {
        if !pc.alive { return; }
        if pc.kind == Kind::Continental {
            pc.thick += deltas[i];
            // Slow net loss of high-standing crust to the oceans (sediment leaving the continent).
            if pc.thick > 38.0 { pc.thick -= (pc.thick - 38.0) / (4.0 * p.erosion_tau) * dt; }
            // Stretching at a divergent boundary thins the margin, most at the coast, tapering inland
            // (a graded passive margin: shelf, slope, hinge).
            if let Some(Some(b)) = binfo.get(i) {
                if b.conv < -CONV_EPS && b.dist < r_margin {
                    pc.thick -= p.thin_coeff * (-b.conv) * dt * (1.0 - b.dist / r_margin).max(0.15);
                    if pc.thick < 20.0 { pc.thick = 20.0; }
                }
            }
            if pc.thick > 70.0 { pc.thick = 70.0; }
        }
        // Arc volcanism on the upper plate: thickens continental crust, builds island arcs on oceanic crust.
        if pc.arc_t == t {
            if pc.kind == Kind::Continental { pc.thick += p.arc_rate * dt; } else { pc.volc += 2.5 * p.arc_rate * dt; }
        }
        pc.volc -= pc.volc / p.volc_tau * dt;
        pc.volc = pc.volc.clamp(0.0, 6.0);

        // Hotspots (Rule X): build volcanic piles, leave a thermal/weakness imprint.
        for h in &hs {
            if dist(*h, pc.pos) < r_hot {
                pc.hot_t = t;
                pc.volc += if pc.kind == Kind::Oceanic { p.hot_rate } else { 0.3 * p.hot_rate } * dt;
            }
        }
    });
    // Failed rifts fill: an intraplate strip of ocean floor enclosed by continental crust subsides under
    // sediment into a shallow basin (aulacogen) - it becomes thin continental crust, not a deep slot.
    {
        let r_n = 1.6 * s;
        let binfo = &w.binfo;
        let fill: Vec<usize> = (0..w.parcels.len()).into_par_iter().filter(|&i| {
            let pc = &w.parcels[i];
            if !pc.alive || pc.kind != Kind::Oceanic || t - pc.birth < 15.0 { return false; }
            if let Some(Some(_)) = binfo.get(i) { return false; } // at a live plate boundary
            let (mut n_all, mut n_cont) = (0usize, 0usize);
            let mut resultant = [0.0; 3];
            w.hash.query(pc.pos, r_n, |k| {
                let pk = &w.parcels[k as usize];
                if k as usize != i && pk.alive && dist(pk.pos, pc.pos) < r_n {
                    n_all += 1;
                    if pk.kind == Kind::Continental && pk.plate == pc.plate {
                        n_cont += 1;
                        resultant = add(resultant, tangent_toward(pc.pos, pk.pos));
                    }
                }
            });
            // Enclosed, not merely bordered: continental neighbours must lie on opposing sides
            // (short resultant), otherwise every passive margin would creep seaward row by row.
            n_cont >= 4 && n_cont * 2 >= n_all && norm(resultant) / (n_cont as f64) < 0.45
        }).collect();
        for i in fill {
            let pc = &mut w.parcels[i];
            pc.kind = Kind::Continental;
            pc.thick = 28.0;
            pc.volc = 0.0;
            pc.birth = t;
            pc.trench_t = NEVER;
        }
    }
    // Accretionary / forearc growth at continental arcs: the margin slowly builds seaward (juvenile crust).
    {
        let mut grow: Vec<usize> = vec![];
        {
            let hash = &w.hash;
            let parcels = &w.parcels;
            for (i, pc) in parcels.iter().enumerate() {
                if !(pc.alive && pc.kind == Kind::Continental && pc.arc_t == t && pc.thick >= 33.0) { continue; }
                if w.rng_f64_hash(i) >= p.accrete_rate * dt { continue; }
                let mut best: Option<(usize, f64)> = None;
                hash.query(pc.pos, 1.4 * s, |k| {
                    let pk = &parcels[k as usize];
                    if pk.alive && pk.kind == Kind::Oceanic && pk.plate == pc.plate {
                        let d = dist(pk.pos, pc.pos);
                        if d < 1.4 * s && best.map_or(true, |(_, bd)| d < bd) { best = Some((k as usize, d)); }
                    }
                });
                if let Some((k, _)) = best { grow.push(k); }
            }
        }
        for k in grow {
            let pk = &mut w.parcels[k];
            if pk.kind != Kind::Oceanic { continue; }
            pk.kind = Kind::Continental;
            pk.thick = 33.5;
            pk.volc = 0.0;
            pk.birth = t;
            w.stats.cont_grown += 1;
        }
    }
    // A mature island arc is new continental crust (continents grow at arcs). Conversion happens in
    // contiguous segments: a parcel with a big pile converts together with its arc neighbours.
    let mut mature: Vec<usize> = vec![];
    {
        let hash = &w.hash;
        let parcels = &w.parcels;
        for (i, pc) in parcels.iter().enumerate() {
            if !(pc.alive && pc.kind == Kind::Oceanic && pc.arc_t == t && pc.volc >= 5.0) { continue; }
            let mut seg = vec![i];
            hash.query(pc.pos, 1.3 * s, |k| {
                let pk = &parcels[k as usize];
                if k as usize != i && pk.alive && pk.kind == Kind::Oceanic && pk.plate == pc.plate && pk.volc >= 2.0 && dist(pk.pos, pc.pos) < 1.3 * s {
                    seg.push(k as usize);
                }
            });
            if seg.len() >= (3.0 * w.n_scale.sqrt()).round() as usize { mature.extend(seg); }
        }
    }
    for i in mature {
        let pc = &mut w.parcels[i];
        if pc.kind != Kind::Oceanic { continue; }
        w.stats.cont_grown += 1;
        let depth = (2600.0 + 350.0 * (t - pc.birth).max(0.0).sqrt()).min(5700.0);
        pc.kind = Kind::Continental;
        pc.thick = (32.8 + (pc.volc * 1000.0 - depth) / 180.0).clamp(30.0, 38.0);
        pc.volc = 0.0;
        pc.birth = t;
    }
}

fn rifting(w: &mut World) {
    advance_rifts(w);
    let dt = w.p.dt;
    let np = w.plates.len();
    // reference pull: slab pull on a 3000 km trench
    let pull_ref = w.p.k_slab * (3000.0 / R_KM);
    for a in 0..np {
        let pl = w.plates[a].clone();
        if !pl.alive || (pl.n as f64) < 300.0 * w.n_scale || (pl.n_cont as f64) < 150.0 * w.n_scale { continue; }
        let cf = pl.n_cont as f64 / pl.n as f64;
        if cf < 0.15 { continue; }
        // Tension from the forces actually pulling on the plate: slab pull on its own slabs plus slab
        // suction toward the trenches it overrides (Rule V). A stagnant plate still creeps toward failure.
        let pulled = (0.35 + (pl.slab + pl.suction) / pull_ref).clamp(0.35, 2.5);
        let size = (pl.n_cont as f64 / (1500.0 * w.n_scale)).sqrt().clamp(0.5, 2.0);
        let rate = pulled * size * (0.4 + 0.6 * cf) / 320.0;
        w.plates[a].tension += rate * dt;
        // Sutures and hot spots weaken the lithosphere (Rules X, XI).
        let weak = pl.n_weak as f64 / pl.n_cont.max(1) as f64;
        let threshold = w.p.rift_threshold / (1.0 + 3.0 * weak);
        let busy = w.rifts.iter().any(|r| r.plate == a as u32);
        if !busy && w.plates[a].tension > threshold && w.rng.gen::<f64>() < w.p.rift_rate * dt {
            if nucleate_rift(w, a) { w.plates[a].tension = 0.0; }
        }
    }
}

/// Mean weakness of the continental crust around `p` (recent sutures, hotspot passages, thin crust).
fn weakness_at(w: &World, p: V3) -> f64 {
    let t = w.t;
    let r = 1.2 * w.s;
    let (mut sc, mut n) = (0.0, 0.0);
    w.hash.query(p, r, |k| {
        let pk = &w.parcels[k as usize];
        if pk.alive && pk.kind == Kind::Continental && dist(pk.pos, p) < r {
            n += 1.0;
            if t - pk.suture_t < 500.0 { sc += 2.0; }
            if t - pk.hot_t < 30.0 { sc += 1.5; }
            if pk.thick < 31.0 { sc += 1.0; }
        }
    });
    if n > 0.0 { sc / n } else { 0.0 }
}

/// Start a rift at a weak continental parcel, running perpendicular to the plate's motion.
fn nucleate_rift(w: &mut World, a: usize) -> bool {
    let t = w.t;
    let idxs: Vec<usize> = w.parcels.iter().enumerate()
        .filter(|(_, pc)| pc.alive && pc.plate == a as u32 && pc.kind == Kind::Continental)
        .map(|(i, _)| i).collect();
    if idxs.is_empty() { return false; }
    let weights: Vec<f64> = idxs.iter().map(|&i| {
        let pc = &w.parcels[i];
        1.0 + if t - pc.suture_t < 500.0 { 5.0 } else { 0.0 } + if t - pc.hot_t < 30.0 { 3.0 } else { 0.0 }
    }).collect();
    let total: f64 = weights.iter().sum();
    let mut pick = w.rng.gen::<f64>() * total;
    let mut nucleus = idxs[0];
    for (k, &i) in idxs.iter().enumerate() {
        pick -= weights[k];
        if pick <= 0.0 { nucleus = i; break; }
    }
    let c = w.parcels[nucleus].pos;
    let vc = surface_velocity(w.plates[a].omega, c);
    let vn = norm(vc);
    let vhat = if vn > 1e-6 { scale(vc, 1.0 / vn) } else { any_tangent(c) };
    let u = normalize(cross(c, vhat));
    w.rifts.push(ActiveRift { plate: a as u32, nucleus: c, normal: vhat, path: vec![c], tip: [c, c], dir: [u, scale(u, -1.0)], done: [false, false], born: t });
    mark_rift(w, c);
    true
}

/// A rift valley forms where the tip passes: flag the parcels and drop them by ~1.5 km of crust.
fn mark_rift(w: &mut World, p: V3) {
    let t = w.t;
    let r = 0.9 * w.s;
    let mut ks = vec![];
    w.hash.query(p, r, |k| {
        let pk = &w.parcels[k as usize];
        if pk.alive && pk.kind == Kind::Continental && dist(pk.pos, p) < r { ks.push(k as usize); }
    });
    for k in ks {
        let pc = &mut w.parcels[k];
        if t - pc.rift_t > 5.0 { pc.thick = (pc.thick - 1.5).max(24.0); }
        pc.rift_t = t;
    }
}

/// Advance every propagating rift tip along the weakest nearby path, roughly perpendicular to plate
/// motion; when both tips have left the continent, split the plate along the path.
fn advance_rifts(w: &mut World) {
    let dt = w.p.dt;
    let t = w.t;
    let step = w.km(w.p.rift_prop_v) * dt;
    let mut i = 0;
    while i < w.rifts.len() {
        let plate = w.rifts[i].plate as usize;
        if !w.plates[plate].alive { w.rifts.remove(i); continue; }
        for e in 0..2 {
            if w.rifts[i].done[e] { continue; }
            let tip = w.rifts[i].tip[e];
            let dir = w.rifts[i].dir[e];
            let side = normalize(cross(tip, dir));
            let mut best: Option<(f64, V3, V3)> = None;
            for k in -4i32..=4 {
                let ang = k as f64 * 15f64.to_radians();
                let d = normalize(add(scale(dir, ang.cos()), scale(side, ang.sin())));
                let np = move_along(tip, d, step);
                let nd = normalize(sub(d, scale(np, dot(d, np))));
                let weak = weakness_at(w, np);
                let v = surface_velocity(w.plates[plate].omega, np);
                let vn = norm(v);
                let vh = if vn > 1e-6 { scale(v, 1.0 / vn) } else { w.rifts[i].normal };
                let jitter = 0.3 * w.rng.gen::<f64>();
                let score = 2.0 * ang.cos() + weak - dot(nd, vh).abs() + jitter;
                if best.map_or(true, |(bs, _, _)| score > bs) { best = Some((score, np, nd)); }
            }
            let (_, np, nd) = best.unwrap();
            let (mut inside, mut foreign) = (false, false);
            {
                let r = 1.2 * w.s;
                let r2 = 1.0 * w.s;
                w.hash.query(np, r, |k| {
                    let pk = &w.parcels[k as usize];
                    if !pk.alive { return; }
                    let d = dist(pk.pos, np);
                    if d < r && pk.plate == plate as u32 && pk.kind == Kind::Continental { inside = true; }
                    if d < r2 && pk.plate != plate as u32 { foreign = true; }
                });
            }
            w.rifts[i].tip[e] = np;
            w.rifts[i].dir[e] = nd;
            w.rifts[i].path.push(np);
            if !inside || foreign { w.rifts[i].done[e] = true; } else { mark_rift(w, np); }
        }
        if w.rifts[i].done[0] && w.rifts[i].done[1] {
            let r = w.rifts.remove(i);
            if split_along(w, &r) { w.stats.rifts += 1; }
            continue;
        }
        if t - w.rifts[i].born > 80.0 { w.rifts.remove(i); continue; } // stalled: a failed rift
        i += 1;
    }
}

/// Split plate `r.plate` along the rift path (extended straight through any oceanic crust to the plate
/// edge). The smaller side becomes a new plate and the halves are pushed apart.
fn split_along(w: &mut World, r: &ActiveRift) -> bool {
    let a = r.plate as usize;
    let s = w.s;
    let t = w.t;
    let mut path = r.path.clone();
    // this plate's current boundary parcels: the extension through the ocean steers toward the nearest
    // one, so the cut takes the shortest plausible route to an existing boundary (a transform link)
    let bnd: Vec<V3> = w.parcels.iter().enumerate()
        .filter(|(i, pc)| pc.alive && pc.plate == a as u32 && matches!(w.binfo.get(*i), Some(Some(b)) if b.dist < 1.5 * s))
        .map(|(_, pc)| pc.pos).collect();
    for e in 0..2 {
        let mut p = r.tip[e];
        let mut d = r.dir[e];
        for _ in 0..600 {
            if !bnd.is_empty() {
                let mut best = (f64::MAX, bnd[0]);
                for &q in &bnd { let dq = dist(q, p); if dq < best.0 { best = (dq, q); } }
                if best.0 > 1.5 * s {
                    let tb = tangent_toward(p, best.1);
                    d = normalize(add(scale(d, 0.75), scale(tb, 0.25)));
                }
            }
            p = move_along(p, d, 0.8 * s);
            d = normalize(sub(d, scale(p, dot(d, p))));
            let mut on_plate = false;
            w.hash.query(p, 1.0 * s, |k| { let pk = &w.parcels[k as usize]; if pk.alive && pk.plate == a as u32 && dist(pk.pos, p) < 1.0 * s { on_plate = true; } });
            if !on_plate { break; }
            path.push(p);
        }
    }
    let mut barrier: std::collections::HashSet<usize> = std::collections::HashSet::new();
    for &p in &path {
        w.hash.query(p, 0.9 * s, |k| { let pk = &w.parcels[k as usize]; if pk.alive && pk.plate == a as u32 && dist(pk.pos, p) < 0.9 * s { barrier.insert(k as usize); } });
    }
    let n0 = normalize(cross(r.nucleus, r.dir[0]));
    let seed_pt = move_along(r.nucleus, n0, 1.8 * s);
    let mut seed: Option<(usize, f64)> = None;
    w.hash.query(seed_pt, 2.0 * s, |k| {
        let pk = &w.parcels[k as usize];
        if pk.alive && pk.plate == a as u32 && !barrier.contains(&(k as usize)) {
            let d = dist(pk.pos, seed_pt);
            if seed.map_or(true, |(_, bd)| d < bd) { seed = Some((k as usize, d)); }
        }
    });
    let Some((seed, _)) = seed else { return false };
    let mut side_a: std::collections::HashSet<usize> = std::collections::HashSet::new();
    let mut stack = vec![seed];
    side_a.insert(seed);
    while let Some(i) = stack.pop() {
        let pi = w.parcels[i].pos;
        let mut next = vec![];
        w.hash.query(pi, 1.5 * s, |k| {
            let pk = &w.parcels[k as usize];
            if pk.alive && pk.plate == a as u32 && !barrier.contains(&(k as usize)) && dist(pk.pos, pi) < 1.5 * s { next.push(k as usize); }
        });
        for k in next { if side_a.insert(k) { stack.push(k); } }
    }
    let side_b: Vec<usize> = w.parcels.iter().enumerate()
        .filter(|(i, pc)| pc.alive && pc.plate == a as u32 && !side_a.contains(i) && !barrier.contains(i))
        .map(|(i, _)| i).collect();
    let min_n = (50.0 * w.n_scale).round() as usize;
    if side_a.len() < min_n || side_b.len() < min_n { return false; }
    let a_is_new = side_a.len() <= side_b.len();
    let mut moved: std::collections::HashSet<usize> = if a_is_new { side_a.clone() } else { side_b.iter().copied().collect() };
    // barrier parcels join whichever side their nearest non-barrier neighbour is on
    let mut adopt = vec![];
    for &k in &barrier {
        let pk = w.parcels[k].pos;
        let mut best: Option<(usize, f64)> = None;
        w.hash.query(pk, 2.5 * s, |q| {
            let pq = &w.parcels[q as usize];
            if pq.alive && pq.plate == a as u32 && !barrier.contains(&(q as usize)) {
                let d = dist(pq.pos, pk);
                if best.map_or(true, |(_, bd)| d < bd) { best = Some((q as usize, d)); }
            }
        });
        if let Some((q, _)) = best { if moved.contains(&q) { adopt.push(k); } }
    }
    for k in adopt { moved.insert(k); }
    let old_n = w.plates[a].n;
    let new_id = w.plates.len() as u32;
    let mut n_cont = 0;
    for &i in &moved {
        w.parcels[i].plate = new_id;
        if w.parcels[i].kind == Kind::Continental { n_cont += 1; }
    }
    let om = w.plates[a].omega;
    let mv = w.plates[a].mean_v;
    w.plates.push(Plate { omega: om, alive: true, tension: 0.0, n: moved.len(), n_cont, n_weak: 0, mean_v: mv, slab: 0.0, suction: 0.0, born: t });
    let parent_rot = w.rot[a];
    w.rot.push(parent_rot);
    // Initial separation (thermal uplift / plume push, Rule X) so ridge push can take over: 8 km/Myr
    // across the rift, shared between the halves in inverse proportion to their size.
    let dir_new = if a_is_new { n0 } else { scale(n0, -1.0) };
    let kick = scale(cross(r.nucleus, dir_new), 8.0 / R_KM);
    let frac = moved.len() as f64 / old_n.max(1) as f64;
    w.plates[new_id as usize].omega = add(om, scale(kick, 1.0 - frac));
    w.plates[a].omega = sub(om, scale(kick, frac));
    w.rift_pairs.insert(((a as u32).min(new_id), (a as u32).max(new_id)), t);
    w.plates[a].tension = 0.0;
    w.plates[a].n = old_n.saturating_sub(moved.len());
    w.plates[a].n_cont = w.plates[a].n_cont.saturating_sub(n_cont);
    true
}

/// Arc detachment and slab rollback (Rules VII-IX). Where an old, dense slab subducts, the arc can
/// detach from its upper plate as a small plate of its own; slab suction then rolls it trenchward at a
/// mantle-limited speed and a back-arc basin of new crust opens behind it. Back-arc basins stay small
/// (Rule VIII): the arc plate re-welds to its parent after `backarc_myr`, or when the slab gets young.
fn backarc(w: &mut World) {
    let t = w.t;
    let dt = w.p.dt;
    let s = w.s;
    // 1. retire arc plates
    let mut retire: Vec<(u32, u32)> = vec![];
    for (&arc, &(parent, _lower, born)) in &w.arc_plates {
        if !w.plates[arc as usize].alive { retire.push((arc, u32::MAX)); continue; }
        if t - born > w.p.backarc_myr || !w.plates[parent as usize].alive { retire.push((arc, parent)); }
    }
    for (arc, parent) in retire {
        w.arc_plates.remove(&arc);
        if parent == u32::MAX || !w.plates[parent as usize].alive || !w.plates[arc as usize].alive { continue; }
        for pc in w.parcels.iter_mut() { if pc.alive && pc.plate == arc { pc.plate = parent; } }
        w.plates[arc as usize].alive = false;
        w.stats.merges += 1;
    }
    // 2. candidate trenches: per (lower, upper) pair, contact parcels on the lower plate and mean slab age
    let mut pairs: HashMap<(u32, u32), (Vec<V3>, f64)> = HashMap::new();
    for (i, pc) in w.parcels.iter().enumerate() {
        if !pc.alive || pc.kind != Kind::Oceanic { continue; }
        let Some(Some(b)) = w.binfo.get(i) else { continue };
        if b.conv <= CONV_EPS || b.dist >= 1.5 * s { continue; }
        let key = (pc.plate.min(b.other), pc.plate.max(b.other));
        match w.polarity.get(&key) { Some(&sp) if sp == pc.plate => {}, _ => continue }
        let e = pairs.entry((pc.plate, b.other)).or_insert((vec![], 0.0));
        e.0.push(pc.pos);
        e.1 += t - pc.birth;
    }
    let min_contact = (20.0 * w.n_scale.sqrt()).round() as usize;
    let min_n = (40.0 * w.n_scale).round() as usize;
    let r_arc = w.km(w.p.backarc_km);
    let mut cands: Vec<((u32, u32), Vec<V3>)> = vec![];
    for ((lower, upper), (pts, age_sum)) in pairs {
        if pts.len() < min_contact { continue; }
        if age_sum / (pts.len() as f64) < w.p.rollback_age { continue; }
        if w.arc_plates.contains_key(&upper) { continue; }
        if w.arc_plates.values().any(|&(par, low, _)| par == upper && low == lower) { continue; }
        if !w.plates[upper as usize].alive { continue; }
        if w.rng.gen::<f64>() >= w.p.rollback_rate * dt { continue; }
        cands.push(((lower, upper), pts));
    }
    for ((lower, upper), pts) in cands {
        // the arc sliver: upper-plate parcels within backarc_km of the trench contact
        let mut set: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for &q in &pts {
            w.hash.query(q, r_arc, |k| {
                let pk = &w.parcels[k as usize];
                if pk.alive && pk.plate == upper && dist(pk.pos, q) < r_arc { set.insert(k as usize); }
            });
        }
        let upper_n = w.plates[upper as usize].n.max(1);
        if set.len() < min_n || set.len() * 2 > upper_n { continue; }
        let new_id = w.plates.len() as u32;
        let mut n_cont = 0;
        for &i in &set {
            w.parcels[i].plate = new_id;
            if w.parcels[i].kind == Kind::Continental { n_cont += 1; }
        }
        let om = w.plates[upper as usize].omega;
        let mv = w.plates[upper as usize].mean_v;
        w.plates.push(Plate { omega: om, alive: true, tension: 0.0, n: set.len(), n_cont, n_weak: 0, mean_v: mv, slab: 0.0, suction: 0.0, born: t });
        let parent_rot = w.rot[upper as usize];
        w.rot.push(parent_rot);
        w.polarity.insert((lower.min(new_id), lower.max(new_id)), lower);
        w.rift_pairs.insert((upper.min(new_id), upper.max(new_id)), t);
        w.arc_plates.insert(new_id, (upper, lower, t));
        w.plates[upper as usize].n = upper_n.saturating_sub(set.len());
        w.stats.backarcs += 1;
    }
}

fn resolve(remap: &[u32], mut x: u32) -> u32 {
    while remap[x as usize] != x { x = remap[x as usize]; }
    x
}

fn merge_and_cleanup(w: &mut World) {
    let s = w.s;
    let t = w.t;
    // Boundary census per plate pair: (boundary parcels, continent-continent contacts, sum |v_rel|).
    let mut pairs: HashMap<(u32, u32), (usize, usize, f64)> = HashMap::new();
    w.pair_ccf.clear();
    for (i, pc) in w.parcels.iter().enumerate() {
        if !pc.alive { continue; }
        let b = match w.binfo.get(i) { Some(Some(b)) => *b, _ => continue };
        if b.dist >= 1.5 * s || b.other == pc.plate { continue; }
        let pj = &w.parcels[b.oidx as usize];
        if pj.plate != b.other || !pj.alive { continue; }
        let key = (pc.plate.min(b.other), pc.plate.max(b.other));
        let e = pairs.entry(key).or_insert((0, 0, 0.0));
        e.0 += 1;
        if pc.kind == Kind::Continental && pj.kind == Kind::Continental { e.1 += 1; }
        let vi = surface_velocity(w.plates[pc.plate as usize].omega, pc.pos);
        let vj = surface_velocity(w.plates[pj.plate as usize].omega, pj.pos);
        e.2 += norm(sub(vi, vj));
    }
    // Suturing (Rule IV/XI): a collided, mostly continent-continent, nearly static boundary welds the plates.
    let mut merges: Vec<(u32, u32)> = vec![];
    let dt = w.p.dt;
    for (&(a, b), &(nb, ncc, vsum)) in &pairs {
        let ccf = (ncc as f64) / (nb as f64);
        w.pair_ccf.insert((a, b), ccf);
        if (nb as f64) < 10.0 * w.n_scale.sqrt() { continue; }
        let vrel = vsum / nb as f64;
        // Track how long this contact has been static: a boundary with no relative motion is not a
        // boundary (failed rift, locked collision), and the two plates are really one plate.
        let st = w.static_myr.entry((a, b)).or_insert(0.0);
        if vrel < 4.0 { *st += dt; } else { *st = 0.0; }
        let static_long = *st >= 40.0;
        // A freshly rifted pair gets time to separate before it can be declared a failed rift.
        if t - w.plates[a as usize].born < 80.0 || t - w.plates[b as usize].born < 80.0 { continue; }
        // Weld when the contact is a locked, mostly continental collision, or when it has simply stopped
        // moving for 40 Myr. Active partial collisions are handled locally by block accretion instead.
        let locked = w.polarity.contains_key(&(a, b)) && ccf > 0.6 && vrel < 20.0;
        if locked || static_long {
            if w.plates[a as usize].n <= w.plates[b as usize].n { merges.push((a, b)); } else { merges.push((b, a)); }
        }
    }
    let mut remap: Vec<u32> = (0..w.plates.len() as u32).collect();
    for (from, into) in merges {
        let f = resolve(&remap, from);
        let i = resolve(&remap, into);
        if f == i { continue; }
        remap[f as usize] = i;
        w.plates[f as usize].alive = false;
        w.stats.merges += 1;
    }
    if w.stats.merges > 0 {
        for i in 0..w.parcels.len() {
            if !w.parcels[i].alive { continue; }
            let pl = w.parcels[i].plate;
            let np = resolve(&remap, pl);
            if let Some(Some(b)) = w.binfo.get(i) {
                if b.other != pl && resolve(&remap, b.other) == np && b.dist < 1.5 * s { w.parcels[i].suture_t = t; }
            }
            w.parcels[i].plate = np;
        }
    }
    // Plates consumed to almost nothing (ridge subduction, Last Rule) dissolve into their neighbours.
    let mut counts = vec![0usize; w.plates.len()];
    for pc in &w.parcels { if pc.alive { counts[pc.plate as usize] += 1; } }
    let tiny_n = (30.0 * w.n_scale).round() as usize;
    let tiny: Vec<u32> = (0..w.plates.len()).filter(|&a| w.plates[a].alive && counts[a] < tiny_n).map(|a| a as u32).collect();
    if !tiny.is_empty() {
        let mut adopt: Vec<(usize, u32)> = vec![];
        {
            let hash = &w.hash;
            let parcels = &w.parcels;
            let plates = &w.plates;
            for (i, pc) in parcels.iter().enumerate() {
                if !pc.alive || !tiny.contains(&pc.plate) { continue; }
                let mut best: Option<(u32, f64)> = None;
                hash.query(pc.pos, 3.0 * s, |k| {
                    let pk = &parcels[k as usize];
                    if pk.alive && pk.plate != pc.plate && plates[pk.plate as usize].alive && !tiny.contains(&pk.plate) {
                        let d = dist(pk.pos, pc.pos);
                        if best.map_or(true, |(_, bd)| d < bd) { best = Some((pk.plate, d)); }
                    }
                });
                if let Some((pl, _)) = best { adopt.push((i, pl)); }
            }
        }
        for (i, pl) in adopt { w.parcels[i].plate = pl; }
    }
    let alive: Vec<bool> = w.plates.iter().map(|p| p.alive).collect();
    w.polarity.retain(|&(a, b), _| a != b && alive[a as usize] && alive[b as usize]);
    w.pair_absorbed.retain(|&(a, b), _| a != b && alive[a as usize] && alive[b as usize]);
    w.pair_compress.retain(|&(a, b), _| a != b && alive[a as usize] && alive[b as usize]);
    w.arc_plates.retain(|&arc, &mut (par, low, _)| alive[arc as usize] && alive[par as usize] && alive[low as usize]);
    w.static_myr.retain(|&(a, b), _| a != b && alive[a as usize] && alive[b as usize]);
    let (t_now, push_myr) = (w.t, w.p.rift_push_myr);
    w.rift_pairs.retain(|&(a, b), &mut tr| a != b && alive[a as usize] && alive[b as usize] && t_now - tr < push_myr + 1.0);
}

pub fn plate_stats(w: &mut World) {
    for pl in w.plates.iter_mut() { pl.n = 0; pl.n_cont = 0; pl.n_weak = 0; pl.mean_v = 0.0; }
    let t = w.t;
    let (mut total, mut cont) = (0usize, 0usize);
    let (mut vsum, mut vmax) = (0.0f64, 0.0f64);
    for pc in &w.parcels {
        if !pc.alive { continue; }
        let pl = &mut w.plates[pc.plate as usize];
        pl.n += 1;
        total += 1;
        if pc.kind == Kind::Continental {
            pl.n_cont += 1;
            cont += 1;
            if t - pc.suture_t < 500.0 || t - pc.hot_t < 30.0 { pl.n_weak += 1; }
        }
        let v = norm(surface_velocity(pl.omega, pc.pos));
        pl.mean_v += v;
        vsum += v;
        if v > vmax { vmax = v; }
    }
    for pl in w.plates.iter_mut() {
        if pl.n > 0 { pl.mean_v /= pl.n as f64; } else { pl.alive = false; }
    }
    w.stats.mean_v = vsum / total.max(1) as f64;
    w.stats.max_v = vmax;
    w.stats.cont_frac = cont as f64 / total.max(1) as f64;
    w.stats.n_plates = w.plates.iter().filter(|p| p.alive).count();
    w.stats.n_parcels = total;
}
