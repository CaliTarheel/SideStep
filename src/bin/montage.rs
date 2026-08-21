//! Lay out the elevation slices of a run as a contact sheet: `montage <run_dir> [every_myr] [cols]`.
use image::{imageops, Rgb, RgbImage};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let dir = args.get(0).cloned().unwrap_or_else(|| "out/run".into());
    let every: i64 = args.get(1).and_then(|v| v.parse().ok()).unwrap_or(100);
    let cols: u32 = args.get(2).and_then(|v| v.parse().ok()).unwrap_or(3);
    let layer = args.get(3).cloned().unwrap_or_else(|| "elev".into());

    let mut slices: Vec<(i64, String)> = std::fs::read_dir(&dir)
        .expect("run dir")
        .filter_map(|e| e.ok())
        .filter_map(|e| {
            let name = e.file_name().to_string_lossy().to_string();
            name.strip_prefix('t').and_then(|n| n.parse::<i64>().ok()).map(|t| (t, e.path().to_string_lossy().to_string()))
        })
        .collect();
    // keep, for each multiple of `every`, the slice nearest to it
    slices.sort();
    let mut picked: Vec<(i64, String)> = vec![];
    let mut by_mark: std::collections::BTreeMap<i64, (i64, String)> = std::collections::BTreeMap::new();
    for (t, path) in slices {
        let mark = ((t as f64) / every as f64).round() as i64 * every;
        let d = (t - mark).abs();
        if d > (every / 4).max(1) { continue; }
        match by_mark.get(&mark) { Some((bt, _)) if (bt - mark).abs() <= d => {}, _ => { by_mark.insert(mark, (t, path)); } }
    }
    for (_, v) in by_mark { picked.push(v); }
    let slices = picked;
    if slices.is_empty() { eprintln!("no slices in {}", dir); std::process::exit(1); }

    let first = image::open(format!("{}/{}.png", slices[0].1, layer)).expect("open first slice").to_rgb8();
    let (w, h) = (first.width() / 2, first.height() / 2);
    let rows = (slices.len() as u32 + cols - 1) / cols;
    let pad = 6u32;
    let label_h = 14u32;
    let mut sheet = RgbImage::from_pixel(cols * (w + pad) + pad, rows * (h + pad + label_h) + pad, Rgb([20, 20, 24]));
    for (i, (t, path)) in slices.iter().enumerate() {
        let img = image::open(format!("{}/{}.png", path, layer)).expect("open slice").to_rgb8();
        let small = imageops::resize(&img, w, h, imageops::FilterType::Triangle);
        let x = pad + (i as u32 % cols) * (w + pad);
        let y = pad + (i as u32 / cols) * (h + pad + label_h) + label_h;
        imageops::overlay(&mut sheet, &small, x as i64, y as i64);
        draw_label(&mut sheet, x, y - label_h + 2, &format!("{} MYR", ((*t as f64) / every as f64).round() as i64 * every));
    }
    let out = format!("{}/montage_{}.png", dir, layer);
    sheet.save(&out).expect("save montage");
    println!("{}", out);
}

/// Minimal 3x5 pixel font for digits, space and the letters M Y R.
fn draw_label(img: &mut RgbImage, x0: u32, y0: u32, text: &str) {
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
            'M' => [0b101, 0b111, 0b111, 0b101, 0b101],
            'Y' => [0b101, 0b101, 0b010, 0b010, 0b010],
            'R' => [0b110, 0b101, 0b110, 0b101, 0b101],
            _ => [0; 5],
        }
    };
    let scale = 2u32;
    let mut cx = x0;
    for c in text.chars() {
        let g = glyph(c);
        for (row, bits) in g.iter().enumerate() {
            for col in 0..3 {
                if bits & (0b100 >> col) != 0 {
                    for dy in 0..scale { for dx in 0..scale {
                        let px = cx + col * scale + dx;
                        let py = y0 + row as u32 * scale + dy;
                        if px < img.width() && py < img.height() { img.put_pixel(px, py, Rgb([230, 230, 230])); }
                    } }
                }
            }
        }
        cx += 4 * scale;
    }
}
