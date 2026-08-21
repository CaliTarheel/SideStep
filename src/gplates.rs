//! GPlates export.
//!
//! The simulation frame is the absolute (hotspot) reference frame, so every plate's accumulated
//! finite rotation is directly a rotation relative to GPlates' anchored plate 000. "Present day"
//! is the end of the run; ages are `years - t` in Ma.
//!
//! Written into `<out>/gplates/`:
//! - `rotations.rot`            total reconstruction poles for every plate at every slice time
//! - `plates.csv`               plate id, birth/death (Myr and Ma), final parcel counts
//! - `continents_present.gpml`  present-day continental crust as multipoint features per plate
//! - `boundaries_<age>Ma.gmt`   classified boundary points per slice (OGR-GMT, plate id + class)
//! - `rasters/<layer>/<layer>-<age>.png`  time-named copies of the slice images for raster sequences
use crate::geom::*;
use crate::render::classify;
use crate::world::*;
use std::fmt::Write as _;
use std::io::Write;

/// GPlates plate ids: 000 is the anchor, so plate index i becomes 100 + i.
fn gid(i: usize) -> u32 { 100 + i as u32 }

/// Euler pole (lat, lon, angle in degrees) in GPlates' convention. The engine frame
/// (x, y = pole, z, lon = atan2(z, x)) is left-handed with respect to the geographic frame, so a
/// positive rotation angle here is a negative (clockwise) angle in GPlates: flip the sign.
fn pole(m: M3) -> (f64, f64, f64) {
    let (axis, ang) = mat_to_axis_angle(m);
    let (lat, lon) = to_latlon(axis);
    (lat.to_degrees(), lon.to_degrees(), -ang.to_degrees())
}

pub fn export(w: &World) {
    let dir = format!("{}/gplates", w.p.out);
    std::fs::create_dir_all(format!("{}/rasters", dir)).expect("gplates dir");
    let years = w.p.years;
    let age = |t: f64| (years - t).max(0.0);

    // ---- lifetimes from the rotation samples ----
    let np = w.plates.len();
    let mut first: Vec<Option<usize>> = vec![None; np];
    let mut last: Vec<Option<usize>> = vec![None; np];
    for (si, (_, row)) in w.rot_hist.iter().enumerate() {
        for i in 0..np.min(row.len()) {
            if row[i].is_some() { if first[i].is_none() { first[i] = Some(si); } last[i] = Some(si); }
        }
    }

    // ---- rotations.rot ----
    let mut rot = String::new();
    writeln!(rot, "! tectonic engine export: absolute (hotspot) frame, present day = {} Myr of simulation", years).unwrap();
    writeln!(rot, "! moving  age(Ma)  lat  lon  angle  fixed").unwrap();
    for i in 0..np {
        let (Some(f), Some(l)) = (first[i], last[i]) else { continue };
        // reference frame: the plate's last recorded position (present day for living plates)
        let r_ref = w.rot_hist[l].1[i].unwrap();
        let r_ref_t = mat_t(r_ref);
        // GPlates wants ages increasing within a sequence
        for si in (f..=l).rev() {
            let Some(r) = w.rot_hist[si].1[i] else { continue };
            let q = mat_mul(r, r_ref_t);
            let (lat, lon, ang) = pole(q);
            writeln!(rot, "{:>4} {:8.1} {:9.4} {:10.4} {:10.4}  000 !plate {} born {:.0} Myr", gid(i), age(w.rot_hist[si].0), lat, lon, ang, gid(i), w.plates[i].born).unwrap();
        }
    }
    std::fs::write(format!("{}/rotations.rot", dir), rot).expect("rotations.rot");

    // ---- plates.csv ----
    let mut csv = String::from("gplates_id,index,born_myr,died_myr,begin_ma,end_ma,alive_at_end,parcels,continental_parcels\n");
    for i in 0..np {
        let (Some(f), Some(l)) = (first[i], last[i]) else { continue };
        let t0 = w.rot_hist[f].0;
        let t1 = w.rot_hist[l].0;
        writeln!(csv, "{},{},{},{},{},{},{},{},{}", gid(i), i, w.plates[i].born, if w.plates[i].alive { String::from("") } else { format!("{}", t1) },
            age(t0), age(t1), w.plates[i].alive, w.plates[i].n, w.plates[i].n_cont).unwrap();
    }
    std::fs::write(format!("{}/plates.csv", dir), csv).expect("plates.csv");

    // ---- continents_present.gpml ----
    let mut g = String::new();
    g.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    g.push_str("<gpml:FeatureCollection xmlns:gpml=\"http://www.gplates.org/gplates\" xmlns:gml=\"http://www.opengis.net/gml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" gpml:version=\"1.6.0339\" xsi:schemaLocation=\"http://www.gplates.org/gplates ../xsd/gpml.xsd http://www.opengis.net/gml ../xsd/gml.xsd\">\n");
    for i in 0..np {
        if !w.plates[i].alive { continue; }
        let pts: Vec<&Parcel> = w.parcels.iter().filter(|pc| pc.alive && pc.plate == i as u32 && pc.kind == Kind::Continental).collect();
        if pts.is_empty() { continue; }
        let Some(f) = first[i] else { continue };
        let begin = age(w.rot_hist[f].0);
        writeln!(g, "  <gml:featureMember>\n    <gpml:UnclassifiedFeature>").unwrap();
        writeln!(g, "      <gpml:identity>GPlates-tectonic-plate-{}</gpml:identity>\n      <gpml:revision>GPlates-tectonic-rev-{}</gpml:revision>", gid(i), gid(i)).unwrap();
        writeln!(g, "      <gml:name>continental crust of plate {}</gml:name>", gid(i)).unwrap();
        writeln!(g, "      <gpml:reconstructionPlateId><gpml:ConstantValue><gpml:value>{}</gpml:value><gml:description></gml:description><gpml:valueType xmlns:gpml=\"http://www.gplates.org/gplates\">gpml:plateId</gpml:valueType></gpml:ConstantValue></gpml:reconstructionPlateId>", gid(i)).unwrap();
        writeln!(g, "      <gml:validTime><gml:TimePeriod><gml:begin><gml:TimeInstant><gml:timePosition gml:frame=\"http://gplates.org/TRS/flat\">{}</gml:timePosition></gml:TimeInstant></gml:begin><gml:end><gml:TimeInstant><gml:timePosition gml:frame=\"http://gplates.org/TRS/flat\">0</gml:timePosition></gml:TimeInstant></gml:end></gml:TimePeriod></gml:validTime>", begin).unwrap();
        g.push_str("      <gpml:unclassifiedGeometry><gpml:ConstantValue><gpml:value><gml:MultiPoint>\n");
        for pc in pts {
            let (lat, lon) = to_latlon(pc.pos);
            writeln!(g, "        <gml:pointMember><gml:Point><gml:pos>{:.4} {:.4}</gml:pos></gml:Point></gml:pointMember>", lat.to_degrees(), lon.to_degrees()).unwrap();
        }
        g.push_str("      </gml:MultiPoint></gpml:value><gml:description></gml:description><gpml:valueType xmlns:gpml=\"http://www.gplates.org/gplates\">gml:MultiPoint</gpml:valueType></gpml:ConstantValue></gpml:unclassifiedGeometry>\n");
        g.push_str("    </gpml:UnclassifiedFeature>\n  </gml:featureMember>\n");
    }
    g.push_str("</gpml:FeatureCollection>\n");
    std::fs::write(format!("{}/continents_present.gpml", dir), g).expect("gpml");

    // ---- boundaries at present day (OGR-GMT points with plate id and class) ----
    let mut gmt = String::new();
    gmt.push_str("# @VGMT1.0 @GPOINT\n# @Nplate_id|class|class_name\n# @Tinteger|integer|string\n# FEATURE_DATA\n");
    let names = ["none", "trench", "arc", "collision", "ridge", "rift", "transform", "suture", "hotspot"];
    for (i, pc) in w.parcels.iter().enumerate() {
        if !pc.alive { continue; }
        let c = classify(w, i);
        if c == 0 || c > 6 { continue; }
        let (lat, lon) = to_latlon(pc.pos);
        writeln!(gmt, "# @D{}|{}|{}\n{:.4} {:.4}", gid(pc.plate as usize), c, names[c as usize], lon.to_degrees(), lat.to_degrees()).unwrap();
    }
    std::fs::write(format!("{}/boundaries_0Ma.gmt", dir), gmt).expect("gmt");

    // ---- raster sequences: time-named copies ----
    for layer in ["elev", "plates", "age", "bounds"] {
        let ldir = format!("{}/rasters/{}", dir, layer);
        let _ = std::fs::create_dir_all(&ldir);
        for (t, _) in &w.rot_hist {
            let src = format!("{}/t{:05}/{}.png", w.p.out, t.round() as i64, layer);
            let dst = format!("{}/{}-{}.png", ldir, layer, age(*t).round() as i64);
            let _ = std::fs::copy(src, dst);
        }
    }
    let mut readme = std::fs::File::create(format!("{}/README.txt", dir)).expect("readme");
    writeln!(readme, "GPlates export from the tectonic engine.\n\nrotations.rot          load as a rotation file (File > Open Feature Collection). Anchored plate 000 = absolute/hotspot frame; present day = end of run ({} Myr), ages in Ma.\ncontinents_present.gpml present-day continental crust per plate (multipoints, valid from the plate's birth to 0 Ma). Load with the rotation file and drag the time slider.\nboundaries_0Ma.gmt     present-day boundary parcels with plate id and kinematic class (1 trench, 2 arc, 3 collision, 4 ridge, 5 rift, 6 transform).\nrasters/<layer>/       time-named images; import one folder as a time-dependent raster (File > Import > Import Time-Dependent Raster), global extent -180..180 / -90..90.\nplates.csv             plate lifetimes.\n\nLimitations: parcels that changed plate during the run (accreted blocks) are exported under their final plate only; no topological (continuously closing) plate polygons yet.", years).unwrap();
}
