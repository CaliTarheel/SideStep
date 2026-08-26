//! Rasterise the parcel field to equirectangular maps: elevation (hypsometric PNG + raw f32),
//! plate IDs with boundaries, and oceanic crust age.
use crate::geom::*;
use crate::world::*;
use image::{Rgb, RgbImage};
use rayon::prelude::*;
use std::io::Write;

/// Elevation (m) of a parcel relative to current sea level.
pub fn parcel_elev(w: &World, pc: &Parcel) -> f64 { parcel_elev_abs(w, pc) - w.sea_level }

/// Eustatic sea level: the ocean's volume is fixed at its t = 0 value, so sea level rises when ridges
/// are young and voluminous and falls after collisions and ridge subduction (Scotese's last rule).
pub fn update_sea_level(w: &mut World) {
    let es: Vec<f64> = w.parcels.iter().filter(|pc| pc.alive).map(|pc| parcel_elev_abs(w, pc)).collect();
    if es.is_empty() { return; }
    let vol = |z: f64| es.iter().map(|&e| if e < z { z - e } else { 0.0 }).sum::<f64>() / es.len() as f64;
    // Datum: the ocean volume after 100 Myr of spin-up (the random initial sea floor is older and
    // deeper than the steady state the engine settles into).
    let v0 = match w.sea_v0 {
        Some(v) => v,
        None => { if w.t >= 100.0 - 1e-6 { w.sea_v0 = Some(vol(0.0)); } w.sea_level = 0.0; return; }
    };
    let (mut lo, mut hi) = (-4000.0, 4000.0);
    for _ in 0..50 { let mid = 0.5 * (lo + hi); if vol(mid) < v0 { lo = mid; } else { hi = mid; } }
    w.sea_level = 0.5 * (lo + hi);
}

/// Elevation (m) of a parcel from its recorded history, in the fixed datum (sea level at t = 0).
pub fn parcel_elev_abs(w: &World, pc: &Parcel) -> f64 {
    let t = w.t;
    let mut e = match pc.kind {
        Kind::Oceanic => {
            // Age-depth (Parsons & Sclater style), flattening for old lithosphere.
            let age = (t - pc.birth).max(0.0);
            let depth = (2600.0 + 350.0 * age.sqrt()).min(5700.0);
            let mut e = -depth;
            let dtr = t - pc.trench_t;
            if dtr < 10.0 { e -= 3500.0 * pc.trench_w * (1.0 - dtr / 10.0); }
            e
        }
        // Airy isostasy: ~180 m of elevation per km of crust above the ~32.8 km sea-level reference.
        Kind::Continental => 180.0 * (pc.thick - 32.8),
    };
    e += pc.volc * 1000.0;
    let dh = t - pc.hot_t;
    if dh < 40.0 { e += 600.0 * (-dh / 15.0).exp(); }
    e
}

struct Px { elev: f32, plate: u32, kind: u8, age: f32, cls: u8, stress: f32 }

/// Boundary / feature class of a parcel, from its current boundary info and history.
/// 0 none, 1 trench (subducting side), 2 arc (overriding side), 3 continental collision,
/// 4 oceanic ridge, 5 continental rift, 6 transform, 7 recent suture, 8 active hotspot.
pub fn classify(w: &World, i: usize) -> u8 {
    let pc = &w.parcels[i];
    let t = w.t;
    if t - pc.rift_t < 1.5 * w.p.dt { return 5; }
    if let Some(Some(b)) = w.binfo.get(i) {
        if b.dist < 1.5 * w.s && b.other != pc.plate {
            let pj = &w.parcels[b.oidx as usize];
            if b.conv > CONV_EPS {
                if pc.kind == Kind::Continental && pj.kind == Kind::Continental { return 3; }
                let key = (pc.plate.min(b.other), pc.plate.max(b.other));
                return match w.polarity.get(&key) { Some(&s) if s == pc.plate => 1, Some(_) => 2, None => 1 };
            } else if b.conv < -CONV_EPS {
                return if pc.kind == Kind::Continental || pj.kind == Kind::Continental { 5 } else { 4 };
            } else {
                return 6;
            }
        }
    }
    if t - pc.hot_t < 1.5 * w.p.dt { return 8; }
    if t - pc.suture_t < 100.0 { return 7; }
    0
}

pub fn class_color(c: u8) -> Option<[u8; 3]> {
    match c {
        1 => Some([255, 59, 48]),   // trench
        2 => Some([255, 176, 0]),   // arc / upper plate
        3 => Some([208, 80, 255]),  // continental collision
        4 => Some([48, 224, 96]),   // oceanic ridge
        5 => Some([0, 229, 255]),   // continental rift
        6 => Some([235, 235, 235]), // transform
        7 => Some([128, 64, 160]),  // recent suture
        8 => Some([255, 240, 80]),  // active hotspot
        _ => None,
    }
}

pub fn render(w: &World) {
    let wd = w.p.width;
    let ht = wd / 2;
    let s = w.s;
    let r = 1.3 * s;
    let rows: Vec<Vec<Px>> = (0..ht)
        .into_par_iter()
        .map(|y| {
            let lat = (90.0 - (y as f64 + 0.5) / ht as f64 * 180.0).to_radians();
            (0..wd)
                .map(|x| {
                    let lon = (-180.0 + (x as f64 + 0.5) / wd as f64 * 360.0).to_radians();
                    let p = from_latlon(lat, lon);
                    let (mut wsum, mut esum) = (0.0, 0.0);
                    let mut best = (f64::MAX, u32::MAX);
                    // Plate id by weighted vote (robust to a stray interleaved parcel).
                    let mut votes: Vec<(u32, f64, u32)> = Vec::with_capacity(4);
                    w.hash.query(p, r, |j| {
                        let pc = &w.parcels[j as usize];
                        if !pc.alive { return; }
                        let d = dist(pc.pos, p);
                        if d < r {
                            let wgt = (-(d / (0.5 * s)).powi(2)).exp();
                            wsum += wgt;
                            esum += wgt * parcel_elev(w, pc);
                            if d < best.0 { best = (d, j); }
                            match votes.iter_mut().find(|v| v.0 == pc.plate) {
                                Some(v) => { v.1 += wgt; }
                                None => votes.push((pc.plate, wgt, j)),
                            }
                        }
                    });
                    if best.1 == u32::MAX { return Px { elev: 0.0, plate: u32::MAX, kind: 0, age: -1.0, cls: 0, stress: 0.0 }; }
                    let winner = votes.iter().max_by(|a, b| a.1.partial_cmp(&b.1).unwrap()).unwrap().0;
                    let mut bi = best.1;
                    if w.parcels[bi as usize].plate != winner {
                        // nearest parcel of the winning plate
                        let mut bd = f64::MAX;
                        w.hash.query(p, r, |j| {
                            let pc = &w.parcels[j as usize];
                            if pc.alive && pc.plate == winner { let d = dist(pc.pos, p); if d < bd { bd = d; bi = j; } }
                        });
                    }
                    let pc = &w.parcels[bi as usize];
                    let age = if pc.kind == Kind::Oceanic { (w.t - pc.birth) as f32 } else { -1.0 };
                    // Boundary class from the nearest parcel; prefer a classified neighbour if the nearest has none.
                    let mut cls = classify(w, bi as usize);
                    if cls == 0 && best.1 != bi { cls = classify(w, best.1 as usize); }
                    let mut elev = esum / wsum;
                    if w.p.detail {
                        // Sub-parcel detail conditioned on tectonic state: abyssal hills on young sea floor,
                        // rugged relief in young orogens and arcs, smooth cratons, shelves and abyssal plains.
                        let n = 0.55 * w.detail_noise[0].eval(p) + 0.3 * w.detail_noise[1].eval(p) + 0.15 * w.detail_noise[2].eval(p);
                        let amp = if pc.kind == Kind::Oceanic {
                            let a = (w.t - pc.birth).max(0.0);
                            25.0 + 160.0 * (-a / 50.0).exp()
                        } else if elev > 0.0 {
                            let young = (w.t - pc.suture_t < 150.0) || (w.t - pc.arc_t < 20.0) || pc.thick > 40.0;
                            40.0 + 0.12 * elev + if young { 120.0 } else { 0.0 }
                        } else { 12.0 };
                        elev += amp * n;
                    }
                    Px { elev: elev as f32, plate: pc.plate, kind: if pc.kind == Kind::Continental { 1 } else { 0 }, age, cls, stress: pc.stress }
                })
                .collect()
        })
        .collect();

    // Name slices by nominal time so steps that do not divide the slice interval still give round names.
    let nominal = (w.t / w.p.slice_every).round() * w.p.slice_every;
    let dir = format!("{}/t{:05}", w.p.out, nominal.round() as i64);
    std::fs::create_dir_all(&dir).expect("slice dir");

    // Elevation: hypsometric PNG + raw little-endian f32.
    let mut elev_img = RgbImage::new(wd as u32, ht as u32);
    let mut raw: Vec<u8> = Vec::with_capacity(wd * ht * 4);
    for (y, row) in rows.iter().enumerate() {
        for (x, px) in row.iter().enumerate() {
            elev_img.put_pixel(x as u32, y as u32, Rgb(hypso(px.elev)));
            raw.extend_from_slice(&px.elev.to_le_bytes());
        }
    }
    elev_img.save(format!("{}/elev.png", dir)).expect("save elev.png");
    std::fs::write(format!("{}/elev_{}x{}_f32le.raw", dir, wd, ht), &raw).expect("write raw");

    // Plates: colour by id, black at boundaries, darker tint over oceans.
    let mut plate_img = RgbImage::new(wd as u32, ht as u32);
    for y in 0..ht {
        for x in 0..wd {
            let px = &rows[y][x];
            let mut boundary = false;
            for (dx, dy) in [(1i64, 0i64), (-1, 0), (0, 1), (0, -1)] {
                let xx = ((x as i64 + dx).rem_euclid(wd as i64)) as usize;
                let yy = (y as i64 + dy).clamp(0, ht as i64 - 1) as usize;
                if rows[yy][xx].plate != px.plate { boundary = true; }
            }
            let c = if boundary { [0, 0, 0] } else {
                let [r, g, b] = plate_color(px.plate);
                if px.kind == 1 { [r, g, b] } else { [(r as f32 * 0.55) as u8, (g as f32 * 0.55) as u8, (b as f32 * 0.55) as u8] }
            };
            plate_img.put_pixel(x as u32, y as u32, Rgb(c));
        }
    }
    plate_img.save(format!("{}/plates.png", dir)).expect("save plates.png");

    // Oceanic age: bright = young; continents brown.
    let mut age_img = RgbImage::new(wd as u32, ht as u32);
    for y in 0..ht {
        for x in 0..wd {
            let px = &rows[y][x];
            let c = if px.age < 0.0 { [110, 85, 60] } else {
                let a = (px.age / 180.0).clamp(0.0, 1.0);
                let v = (235.0 * (1.0 - a) + 20.0 * a) as u8;
                [v / 2, v, 255 - v / 3]
            };
            age_img.put_pixel(x as u32, y as u32, Rgb(c));
        }
    }
    age_img.save(format!("{}/age.png", dir)).expect("save age.png");

    // Boundaries: continuous boundary lines taken from the plate-ownership map, coloured by the
    // kinematic class of the nearest classified parcel; sutures and hotspots drawn as overlays.
    let mut b_img = RgbImage::new(wd as u32, ht as u32);
    let mut counts = [0usize; 9];
    let line_px = ((wd as f64 / 1024.0).round() as i64).max(1);
    let is_boundary = |x: usize, y: usize| -> bool {
        let p = rows[y][x].plate;
        for (dx, dy) in [(1i64, 0i64), (-1, 0), (0, 1), (0, -1)] {
            let xx = ((x as i64 + dx).rem_euclid(wd as i64)) as usize;
            let yy = (y as i64 + dy).clamp(0, ht as i64 - 1) as usize;
            if rows[yy][xx].plate != p { return true; }
        }
        false
    };
    for y in 0..ht {
        for x in 0..wd {
            let px = &rows[y][x];
            // near a boundary line?
            let mut on_line = false;
            'outer: for dy in -line_px..=line_px {
                for dx in -line_px..=line_px {
                    let xx = ((x as i64 + dx).rem_euclid(wd as i64)) as usize;
                    let yy = (y as i64 + dy).clamp(0, ht as i64 - 1) as usize;
                    if is_boundary(xx, yy) { on_line = true; break 'outer; }
                }
            }
            let mut cls = 0u8;
            if on_line {
                // own class if it is a boundary class, else the most common boundary class nearby
                if (1..=6).contains(&px.cls) { cls = px.cls; } else {
                    let mut tally = [0u32; 9];
                    let r = 3 * line_px;
                    for dy in -r..=r { for dx in -r..=r {
                        let xx = ((x as i64 + dx).rem_euclid(wd as i64)) as usize;
                        let yy = (y as i64 + dy).clamp(0, ht as i64 - 1) as usize;
                        let c = rows[yy][xx].cls;
                        if (1..=6).contains(&c) { tally[c as usize] += 1; }
                    } }
                    let (mut best, mut bn) = (0u8, 0u32);
                    for c in 1..=6 { if tally[c] > bn { bn = tally[c]; best = c as u8; } }
                    cls = best;
                }
                if cls == 0 { cls = 6; } // boundary with no resolved kinematics: draw as transform/neutral
            } else if px.cls == 7 || px.cls == 8 { cls = px.cls; }
            let c = match class_color(cls) {
                Some(c) => { counts[cls as usize] += 1; c }
                None => {
                    let h = hypso(px.elev);
                    let lum = 0.3 * h[0] as f32 + 0.59 * h[1] as f32 + 0.11 * h[2] as f32;
                    [0, 1, 2].map(|i| ((h[i] as f32 * 0.3 + lum * 0.7) * 0.55) as u8)
                }
            };
            b_img.put_pixel(x as u32, y as u32, Rgb(c));
        }
    }
    b_img.save(format!("{}/bounds.png", dir)).expect("save bounds.png");

    // Intraplate tension: hot colours where opposing boundary pulls stretch a plate.
    let mut s_img = RgbImage::new(wd as u32, ht as u32);
    for y in 0..ht {
        for x in 0..wd {
            let px = &rows[y][x];
            let h = hypso(px.elev);
            let lum = (0.3 * h[0] as f32 + 0.59 * h[1] as f32 + 0.11 * h[2] as f32) * 0.45;
            let base = [lum as u8, lum as u8, (lum * 1.15) as u8];
            let v = (px.stress / 1.5).clamp(0.0, 1.0);
            let c = if v <= 0.01 { base } else {
                let warm = if v < 0.5 { lerp([70, 40, 90], [220, 80, 40], v * 2.0) } else { lerp([220, 80, 40], [255, 235, 120], (v - 0.5) * 2.0) };
                lerp(base, warm, (0.25 + 0.75 * v).min(1.0))
            };
            s_img.put_pixel(x as u32, y as u32, Rgb(c));
        }
    }
    s_img.save(format!("{}/stress.png", dir)).expect("save stress.png");

    let mut meta = std::fs::File::create(format!("{}/meta.json", dir)).expect("meta");
    let land = rows.iter().flatten().filter(|p| p.elev > 0.0).count() as f64 / (wd * ht) as f64;
    writeln!(meta, "{{ \"t_myr\": {}, \"width\": {}, \"height\": {}, \"plates\": {}, \"parcels\": {}, \"land_frac_pixels\": {:.4}, \"sea_level_m\": {:.1}, \"elev_raw\": \"elev_{}x{}_f32le.raw\", \"projection\": \"equirectangular, lon -180..180 left->right, lat 90..-90 top->bottom\", \"boundary_pixels\": {{ \"trench\": {}, \"arc\": {}, \"collision\": {}, \"ridge\": {}, \"rift\": {}, \"transform\": {}, \"suture_recent\": {}, \"hotspot_active\": {} }} }}",
        w.t, wd, ht, w.alive_plates(), w.stats.n_parcels, land, w.sea_level, wd, ht, counts[1], counts[2], counts[3], counts[4], counts[5], counts[6], counts[7], counts[8]).unwrap();
}

fn lerp(a: [u8; 3], b: [u8; 3], t: f32) -> [u8; 3] {
    let t = t.clamp(0.0, 1.0);
    [0, 1, 2].map(|i| (a[i] as f32 + (b[i] as f32 - a[i] as f32) * t) as u8)
}

fn hypso(e: f32) -> [u8; 3] {
    if e <= 0.0 {
        let d = (-e / 6000.0).clamp(0.0, 1.0);
        if d < 0.15 { lerp([160, 210, 240], [60, 120, 200], d / 0.15) } else { lerp([60, 120, 200], [6, 18, 60], (d - 0.15) / 0.85) }
    } else if e < 200.0 { lerp([80, 150, 80], [110, 170, 90], e / 200.0) }
    else if e < 1000.0 { lerp([110, 170, 90], [190, 180, 110], (e - 200.0) / 800.0) }
    else if e < 3000.0 { lerp([190, 180, 110], [140, 95, 60], (e - 1000.0) / 2000.0) }
    else { lerp([140, 95, 60], [245, 245, 245], (e - 3000.0) / 2500.0) }
}

fn plate_color(id: u32) -> [u8; 3] {
    if id == u32::MAX { return [0, 0, 0]; }
    let h = ((id as f64 * 0.618_033_988_7) % 1.0) as f32 * 6.0;
    let sat = 0.55 + 0.3 * (((id * 7919) % 3) as f32 / 2.0);
    let val = 0.9;
    let c = val * sat;
    let xx = c * (1.0 - ((h % 2.0) - 1.0).abs());
    let (r, g, b) = match h as u32 { 0 => (c, xx, 0.0), 1 => (xx, c, 0.0), 2 => (0.0, c, xx), 3 => (0.0, xx, c), 4 => (xx, 0.0, c), _ => (c, 0.0, xx) };
    let m = val - c;
    [((r + m) * 255.0) as u8, ((g + m) * 255.0) as u8, ((b + m) * 255.0) as u8]
}
