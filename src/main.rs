mod forces;
mod geom;
mod gplates;
mod render;
mod step;
mod world;

use std::io::Write;
use std::time::Instant;
use world::{Params, World};

fn main() {
    let params = Params::from_args(std::env::args().skip(1).collect());
    std::fs::create_dir_all(&params.out).expect("create output directory");
    std::fs::write(format!("{}/params.json", params.out), params.to_json()).expect("write params.json");
    std::fs::write(format!("{}/params_full.txt", params.out), format!("{:#?}\n", params)).expect("write params_full.txt");
    let t_init = Instant::now();
    let mut w = World::new(params);
    eprintln!(
        "init: {} parcels ({:.0} km spacing), {} plates, cont_frac {:.3}, {:.2}s",
        w.stats.n_parcels, w.s * geom::R_KM, w.alive_plates(), w.stats.cont_frac, t_init.elapsed().as_secs_f64()
    );
    let mut log = std::fs::File::create(format!("{}/log.csv", w.p.out)).expect("log.csv");
    writeln!(log, "t,plates,parcels,cont_frac,mean_v,max_v,subducted,created,accreted,absorbed,cont_lost,cont_grown,rifts,merges,initiations,active_rifts,sea_level,backarcs,arc_plates,weld_static,weld_locked,retired,dissolved,deposited,sediment,stress_p50,stress_p95,stress_max,split_off,enclaves").unwrap();

    let n_steps = (w.p.years / w.p.dt).round() as usize;
    let slice_steps = ((w.p.slice_every / w.p.dt).round() as usize).max(1);
    let t0 = Instant::now();
    for k in 0..=n_steps {
        if k % slice_steps == 0 || k == n_steps { render::update_sea_level(&mut w); render::render(&w); w.record_rotations(); }
        if k == n_steps { break; }
        step::step(&mut w);
        let s = &w.stats;
        writeln!(log, "{},{},{},{:.4},{:.2},{:.2},{},{},{},{},{},{},{},{},{},{},{:.1},{},{},{},{},{},{},{},{:.4},{:.2},{:.2},{:.2},{},{}", w.t, s.n_plates, s.n_parcels, s.cont_frac, s.mean_v, s.max_v,
            s.subducted, s.created, s.accreted, s.absorbed, s.cont_lost, s.cont_grown, s.rifts, s.merges, s.initiations, w.rifts.len(), w.sea_level, s.backarcs, w.arc_plates.len(), s.weld_static, s.weld_locked, s.retired, s.dissolved, s.deposited, w.sediment, s.stress_p50, s.stress_p95, s.stress_max, s.split_off, s.enclaves).unwrap();
        if k % 10 == 0 || s.rifts > 0 || s.merges > 0 {
            eprintln!(
                "t={:6.0} Myr  plates={:3}  parcels={:6}  land={:.3}  v_mean={:5.1} v_max={:5.1} km/Myr  sub={:4} new={:4} dock={:3} absorb={:3}{}{}  [{:.1}s]",
                w.t, s.n_plates, s.n_parcels, s.cont_frac, s.mean_v, s.max_v, s.subducted, s.created, s.accreted, s.absorbed,
                if s.rifts > 0 { format!("  RIFT x{}", s.rifts) } else { String::new() },
                if s.merges > 0 { format!("  SUTURE x{}", s.merges) } else { String::new() },
                t0.elapsed().as_secs_f64()
            );
        }
    }
    gplates::export(&w);
    eprintln!("done: {} Myr in {:.1}s -> {}", w.t, t0.elapsed().as_secs_f64(), w.p.out);
}
