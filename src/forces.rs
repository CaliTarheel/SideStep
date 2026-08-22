//! Force balance (Scotese Rule I & II): plates move because they are pushed (ridges) or pulled (slabs),
//! resisted by mantle drag; continental keels add drag. No inertia: omega solves D*omega = torque.
use crate::geom::*;
use crate::world::*;

pub fn update_omegas(w: &mut World) {
    let np = w.plates.len();
    let s = w.s;
    let p = &w.p;
    let mut torque = vec![[0.0f64; 3]; np];
    let mut drag = vec![[[0.0f64; 3]; 3]; np];
    let mut slab = vec![0.0f64; np];
    let mut suction = vec![0.0f64; np];
    let area = s * s; // each parcel represents ~s^2 steradians
    // Boundary forces act on every parcel within `band` of the other plate. Parcels are consumed at 0.8 s
    // and created at 0.9-1.8 s, so a narrow band's occupancy flickered with each step; a 3 s band always
    // holds ~2.2 rows, and `len` is normalised so the force per unit boundary length is the rule's value
    // (calibrated to the 0.7-row occupancy of the original 1.5 s band) independent of dt and spacing.
    let band = 3.0 * s;
    let len = s * 0.55;

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
                        slab[a] += m;
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
                    suction[a] += m;
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
        torque[a] = add(torque[a], cross(r, f));
    }

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
