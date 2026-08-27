//! Identify plates in a rendered slice: per-plate stats on stdout and a labeled plate map.
//! Usage: plateinfo <run_dir> <t_myr>
use image::{Rgb, RgbImage};
use std::collections::HashMap;

// Must match render.rs::plate_color exactly (including truncations).
fn plate_color(id: u32) -> [u8; 3] {
    let h = ((id as f64 * 0.618_033_988_7) % 1.0) as f32 * 6.0;
    let sat = 0.55 + 0.3 * (((id * 7919) % 3) as f32 / 2.0);
    let val = 0.9;
    let c = val * sat;
    let xx = c * (1.0 - ((h % 2.0) - 1.0).abs());
    let (r, g, b) = match h as u32 { 0 => (c, xx, 0.0), 1 => (xx, c, 0.0), 2 => (0.0, c, xx), 3 => (0.0, xx, c), 4 => (xx, 0.0, c), _ => (c, 0.0, xx) };
    let m = val - c;
    [((r + m) * 255.0) as u8, ((g + m) * 255.0) as u8, ((b + m) * 255.0) as u8]
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let dir = args.get(0).cloned().unwrap_or_else(|| "out/run".into());
    let t: i64 = args.get(1).and_then(|v| v.parse().ok()).unwrap_or(1000);
    let path = format!("{}/t{:05}/plates.png", dir, t);
    let img = image::open(&path).expect("open plates.png").to_rgb8();
    let (w, h) = (img.width() as usize, img.height() as usize);

    // colour -> (plate index, continental?)
    let mut table: HashMap<[u8; 3], (u32, bool)> = HashMap::new();
    // truncated colours collide for distant ids: the FIRST (lowest) id keeps the colour
    for id in 0..2048u32 {
        let c = plate_color(id);
        table.entry(c).or_insert((id, true));
        table.entry([(c[0] as f32 * 0.55) as u8, (c[1] as f32 * 0.55) as u8, (c[2] as f32 * 0.55) as u8]).or_insert((id, false));
    }

    // per plate: pixel count, continental count, per-cell counts for label placement
    const GX: usize = 96;
    const GY: usize = 48;
    let mut count: HashMap<u32, (u64, u64)> = HashMap::new();
    let mut cells: HashMap<(u32, usize), (u64, f64, f64)> = HashMap::new();
    for y in 0..h {
        for x in 0..w {
            let p = img.get_pixel(x as u32, y as u32).0;
            let Some(&(id, cont)) = table.get(&p) else { continue };
            let e = count.entry(id).or_insert((0, 0));
            e.0 += 1;
            if cont { e.1 += 1; }
            let cellid = (y * GY / h) * GX + (x * GX / w);
            let ce = cells.entry((id, cellid)).or_insert((0, 0.0, 0.0));
            ce.0 += 1;
            ce.1 += x as f64;
            ce.2 += y as f64;
        }
    }

    // label point per plate: centroid of its best cell (lands inside C/S shapes)
    let mut best_cell: HashMap<u32, (u64, f64, f64)> = HashMap::new();
    for (&(id, _), &(n, sx, sy)) in &cells {
        let e = best_cell.entry(id).or_insert((0, 0.0, 0.0));
        if n > e.0 { *e = (n, sx / n as f64, sy / n as f64); }
    }

    // perimeter (boundary pixels) per plate, for an isoperimetric compactness number
    let mut perim: HashMap<u32, u64> = HashMap::new();
    for y in 0..h {
        for x in 0..w {
            let p0 = img.get_pixel(x as u32, y as u32).0;
            let Some(&(id, _)) = table.get(&p0) else { continue };
            let mut edge = false;
            for (dx, dy) in [(1i64, 0i64), (-1, 0), (0, 1), (0, -1)] {
                let xx = ((x as i64 + dx).rem_euclid(w as i64)) as u32;
                let yy = (y as i64 + dy).clamp(0, h as i64 - 1) as u32;
                let pn = img.get_pixel(xx, yy).0;
                match table.get(&pn) { Some(&(idn, _)) if idn == id => {}, _ => { edge = true; break; } }
            }
            if edge { *perim.entry(id).or_insert(0) += 1; }
        }
    }
    let total: u64 = count.values().map(|v| v.0).sum();
    let mut ids: Vec<(u32, u64, u64)> = count.iter().map(|(&id, &(n, nc))| (id, n, nc)).collect();
    ids.sort_by(|a, b| b.1.cmp(&a.1));

    let mut out = img.clone();
    println!("{:>6} {:>7} {:>7} {:>8} {:>9} {:>9}", "plate", "share%", "cont%", "compact", "label_lat", "label_lon");
    let mut wsum = 0.0;
    let mut csum = 0.0;
    for &(id, n, nc) in &ids {
        if (n as f64) / (total as f64) < 0.0005 { continue; }
        let (_, lx, ly) = best_cell[&id];
        let lat = 90.0 - (ly + 0.5) / h as f64 * 180.0;
        let lon = -180.0 + (lx + 0.5) / w as f64 * 360.0;
        let pm = perim.get(&id).copied().unwrap_or(0) as f64;
        let compact = if pm > 0.0 { (4.0 * std::f64::consts::PI * n as f64 / (pm * pm)).min(1.0) } else { 0.0 };
        wsum += n as f64;
        csum += compact * n as f64;
        println!("{:>6} {:>7.2} {:>7.1} {:>8.3} {:>9.1} {:>9.1}", 100 + id, n as f64 / total as f64 * 100.0, nc as f64 / n as f64 * 100.0, compact, lat, lon);
        draw_label(&mut out, lx as u32, ly as u32, &format!("{}", 100 + id));
    }
    println!("area-weighted mean compactness: {:.3}", csum / wsum.max(1.0));
    let outp = format!("{}/t{:05}/plates_labeled.png", dir, t);
    out.save(&outp).expect("save labeled");
    println!("{}", outp);
}

/// 3x5 digit glyphs drawn at 3x scale with a black outline, centred on (cx, cy).
fn draw_label(img: &mut RgbImage, cx: u32, cy: u32, text: &str) {
    let glyph = |c: char| -> [u8; 5] {
        match c {
            '0' => [0b111, 0b101, 0b101, 0b101, 0b111],
            '1' => [0b010, 0b110, 0b010, 0b010, 0b111],
            '2' => [0b111, 0b001, 0b111, 0b100, 0b111],
            '3' => [0b111, 0b001, 0b111, 0b001, 0b111],
            '4' => [0b101, 0b101, 0b111, 0b001, 0b001],
            '5' => [0b111, 0b100, 0b111, 0b001, 0b111],
            '6' => [0b111, 0b100, 0b111, 0b101, 0b111],
            '7' => [0b111, 0b001, 0b001, 0b001, 0b001],
            '8' => [0b111, 0b101, 0b111, 0b101, 0b111],
            '9' => [0b111, 0b101, 0b111, 0b001, 0b111],
            _ => [0; 5],
        }
    };
    let scale = 3u32;
    let tw = text.len() as u32 * 4 * scale;
    let x0 = cx.saturating_sub(tw / 2);
    let y0 = cy.saturating_sub(5 * scale / 2);
    for pass in 0..2 {
        let col = if pass == 0 { Rgb([0, 0, 0]) } else { Rgb([255, 255, 255]) };
        let offs: &[(i32, i32)] = if pass == 0 { &[(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (1, 1), (-1, 1), (1, -1)] } else { &[(0, 0)] };
        let mut cxx = x0;
        for ch in text.chars() {
            let g = glyph(ch);
            for (row, bits) in g.iter().enumerate() {
                for colb in 0..3 {
                    if bits & (0b100 >> colb) != 0 {
                        for dy in 0..scale {
                            for dx in 0..scale {
                                for &(ox, oy) in offs {
                                    let px = cxx as i32 + (colb * scale + dx) as i32 + ox;
                                    let py = y0 as i32 + (row as u32 * scale + dy) as i32 + oy;
                                    if px >= 0 && py >= 0 && (px as u32) < img.width() && (py as u32) < img.height() {
                                        img.put_pixel(px as u32, py as u32, col);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            cxx += 4 * scale;
        }
    }
}
