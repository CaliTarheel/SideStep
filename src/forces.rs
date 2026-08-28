//! Force balance (Scotese Rule I & II): plates move because they are pushed (ridges) or pulled (slabs),
//! resisted by mantle drag; continental keels add drag. No inertia: omega solves D*omega = torque.
use crate::geom::*;
use crate::world::*;
use std::collections::HashMap;

pub fn update_omegas(w: &mut World) {
    let np = w.plates.len();
    let s = w.s;
    let p = &w.p;
    let mut torque = vec![[0.0f64; 3]; np];
    let mut drag = vec![[[0.0f64; 3]; 3]; np];
    let mut slab = vec![0.0f64; np];
    let mut suction = vec![0.0f64; np];
    let area = s * s; // each parcel represents ~s^2 steradians
    // Boundary forces act on every parcel within `band` of the other plate. The band is km-fixed and
    // equal to the detection reach: with a lattice-fixed band the number of force-carrying rows grew
    // with resolution (200 km holds 1.8 rows at 113 km spacing but 3.6 at 56 km), which doubled the
    // force per unit boundary length at 160K and overheated the stress field.
    let band = w.reach(1.5, 200.0);
    let len = s * 0.55;

    // Normalise boundary forces by a smoothed boundary length: at fine spacing the boundary is rougher,
    // its measured length grows, and per-parcel forces would over-count. Count the contact parcels per
    // ~250 km cell and attenuate when a cell holds more than a straight boundary segment would.
    let cell = 250.0 / R_KM;
    let cdim = (2.0 / cell).ceil() as i32 + 2;
    let ckey = |p: V3| -> i32 {
        let c = |v: f64| ((v + 1.0) / cell).floor() as i32 + 1;
        c(p[0]) + cdim * (c(p[1]) + cdim * c(p[2]))
    };
    let mut occupancy: HashMap<(u32, u32, i32), f64> = HashMap::new();
    let mut cells: HashMap<(u32, i32), (V3, V3, f64)> = HashMap::new();
    for (i, pc) in w.parcels.iter().enumerate() {
        if !pc.alive { continue; }
        let Some(Some(b)) = w.binfo.get(i) else { continue };
        if b.other == pc.plate || b.dist >= band { continue; }
        *occupancy.entry((pc.plate, b.other, ckey(pc.pos))).or_insert(0.0) += 1.0;
    }
    // contact parcels a straight boundary puts in one cell across the whole band, at this spacing
    let expected = (cell / s) * (band / s);

    for (i, pc) in w.parcels.iter().enumerate() {
        if !pc.alive { continue; }
        let a = pc.plate as usize;
        let r = pc.pos;
        let c = if pc.kind == Kind::Continental { p.drag_cont } else { p.drag_ocean } * area;
        for u in 0..3 {
            for v in 0..3 {
                drag[a][u][v] += c * ((if u == v { 1.0 } else { 0.0 }) - r[u] * r[v]);
            }
        }
        let b = match w.binfo.get(i) { Some(Some(b)) => *b, _ => continue };
        if b.other == pc.plate || b.dist >= band { continue; }
        let pj = &w.parcels[b.oidx as usize];
        if pj.plate != b.other { continue; }
        let mut f = [0.0; 3];
        if b.conv > CONV_EPS {
            let key = (pc.plate.min(b.other), pc.plate.max(b.other));
            if let Some(&sub) = w.polarity.get(&key) {
                if sub == pc.plate {
                    if pc.kind == Kind::Oceanic {
                        // Slab pull: older (thicker, denser) slabs pull harder.
                        let age = (w.t - pc.birth).max(0.0);
                        let g = (age / 80.0).clamp(0.2, 1.0);
                        let m = p.k_slab * len * g;
                        f = add(f, scale(b.n, m));
                        let occ = occupancy.get(&(pc.plate, b.other, ckey(r))).copied().unwrap_or(1.0);
                        slab[a] += m * (expected / occ).min(1.0);
                    } else {
                        // Continent entering a trench resists (Rule IV: collision is what stops subduction).
                        f = add(f, scale(b.n, -p.k_coll * len));
                    }
                } else if pc.kind == Kind::Continental && pj.kind == Kind::Continental {
                    f = add(f, scale(b.n, -p.k_coll * len));
                } else if pj.kind == Kind::Oceanic {
                    // Slab suction: the upper plate is drawn toward its trench by the sinking slab's
                    // return flow. This is what lets a trench-ringed supercontinent pull itself apart.
                    let age = (w.t - pj.birth).max(0.0);
                    let g = (age / 80.0).clamp(0.2, 1.0);
                    let m = p.k_suction * len * g;
                    f = add(f, scale(b.n, m));
                    let occ = occupancy.get(&(pc.plate, b.other, ckey(r))).copied().unwrap_or(1.0);
                    suction[a] += m * (expected / occ).min(1.0);
                }
            } else {
                // Convergence at a contact that has not started subducting is resisted: compression
                // builds there until initiation (Rule IV: subduction is hard to start).
                f = add(f, scale(b.n, -p.k_resist * len));
            }
        } else if b.conv < -CONV_EPS {
            // Ridge push: away from the ridge, into the plate. Secondary to slab pull (Rule II: ~20%).
            f = add(f, scale(b.n, -p.k_ridge * len));
        }
        // Young rift: thermal uplift of the shoulders pushes the halves apart until a ridge exists.
        let key = (pc.plate.min(b.other), pc.plate.max(b.other));
        if let Some(&t_rift) = w.rift_pairs.get(&key) {
            if w.t - t_rift < p.rift_push_myr && b.conv > -CONV_EPS {
                f = add(f, scale(b.n, -p.k_rift * len));
            }
        }
        let occ = occupancy.get(&(pc.plate, b.other, ckey(r))).copied().unwrap_or(1.0);
        let norm = (expected / occ).min(1.0);
        torque[a] = add(torque[a], scale(cross(r, f), norm));
        // coarse boundary tractions for the intraplate stress field
        let e = cells.entry((pc.plate, ckey(r))).or_insert(([0.0; 3], [0.0; 3], 0.0));
        e.0 = add(e.0, scale(f, norm));
        e.1 = add(e.1, r);
        e.2 += 1.0;
    }
    let mut ct: HashMap<u32, Vec<(V3, V3)>> = HashMap::new();
    let mut ckeys: Vec<(u32, i32)> = cells.keys().copied().collect();
    ckeys.sort_unstable();
    for k in ckeys {
        let (fsum, psum, n) = cells[&k];
        if norm(fsum) > 1e-12 { ct.entry(k.0).or_default().push((normalize(psum), scale(fsum, 1.0))); let _ = n; }
    }
    w.cell_tractions = ct;

    for a in 0..np {
        if !w.plates[a].alive { continue; }
        let d = drag[a];
        let tr = (d[0][0] + d[1][1] + d[2][2]).max(1e-12);
        let eps = 1e-3 * tr;
        let m = [
            [d[0][0] + eps, d[0][1], d[0][2]],
            [d[1][0], d[1][1] + eps, d[1][2]],
            [d[2][0], d[2][1], d[2][2] + eps],
        ];
        let mut om = solve3(m, torque[a]);
        let v = norm(om) * R_KM;
        // A detached arc rolls back at mantle-flow-limited speed, not at the slab-suction limit.
        let cap = if w.arc_plates.contains_key(&(a as u32)) { p.rollback_v } else { p.v_max };
        if v > cap { om = scale(om, cap / v); }
        // relaxation toward the quasi-static balance, expressed per Myr so the step size does not matter
        let k = 1.0 - (1.0 - p.omega_relax).powf(p.dt);
        let pl = &mut w.plates[a];
        pl.omega = add(scale(pl.omega, 1.0 - k), scale(om, k));
        pl.slab = slab[a];
        pl.suction = suction[a];
    }
}

fn solve3(m: [[f64; 3]; 3], b: V3) -> V3 {
    let det = |m: [[f64; 3]; 3]| -> f64 {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let d = det(m);
    if d.abs() < 1e-18 { return [0.0; 3]; }
    let mut out = [0.0; 3];
    for c in 0..3 {
        let mut mc = m;
        for r in 0..3 { mc[r][c] = b[r]; }
        out[c] = det(mc) / d;
    }
    out
}
