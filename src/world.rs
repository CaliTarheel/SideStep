//! World state: parcels (Lagrangian crust samples), plates (rigid rotations), parameters, initial conditions.
use crate::geom::*;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use std::collections::HashMap;

/// Relative-velocity threshold (km/Myr) separating converging / transform / diverging boundaries.
pub const CONV_EPS: f64 = 2.0;
pub const IDENT: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
pub const NEVER: f64 = -1.0e9;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Kind { Oceanic, Continental }

#[derive(Clone, Copy, Debug)]
pub struct Parcel {
    pub pos: V3,
    pub plate: u32,
    pub kind: Kind,
    /// Time of formation (Myr). Age = t - birth. Continental crust gets a very old birth.
    pub birth: f64,
    /// Crustal thickness (km). Drives continental elevation by Airy isostasy.
    pub thick: f64,
    /// Volcanic pile (km): hotspot seamounts/islands, island arcs, flood basalts.
    pub volc: f64,
    /// Last time this parcel sat next to an active trench on the subducting side (flexural deepening).
    pub trench_t: f64,
    /// Last time this parcel was part of a continental collision (suture = future weakness).
    pub suture_t: f64,
    /// Last time a mantle plume sat under this parcel.
    pub hot_t: f64,
    /// Last time this parcel was on the upper plate right above an active trench (arc volcanism).
    pub arc_t: f64,
    /// Last time a propagating rift passed through this parcel.
    pub rift_t: f64,
    /// Intraplate tension proxy (opposing boundary pulls), updated every `stress_every` Myr.
    pub stress: f32,
    /// Constriction amplification that went into `stress` (1 = full-width plate, >1 = neck).
    pub amp: f32,
    /// Flexural deepening weight (1 at the trench, 0 at the outer edge of the flexural bulge).
    pub trench_w: f64,
    pub alive: bool,
}

#[derive(Clone, Debug)]
pub struct Plate {
    pub omega: V3,
    pub alive: bool,
    /// Accumulated extensional "tension" used by the rifting rule.
    pub tension: f64,
    pub n: usize,
    pub n_cont: usize,
    pub n_weak: usize,
    pub mean_v: f64,
    /// Total slab-pull magnitude acting on this plate last step.
    pub slab: f64,
    /// Total slab-suction magnitude acting on this plate last step.
    pub suction: f64,
    pub born: f64,
}

/// A rift propagating through a continent; when both tips reach the ocean the plate is split along the path.
#[derive(Clone, Debug)]
pub struct ActiveRift {
    pub plate: u32,
    pub nucleus: V3,
    pub normal: V3,
    pub path: Vec<V3>,
    pub tip: [V3; 2],
    pub dir: [V3; 2],
    pub done: [bool; 2],
    pub born: f64,
}

/// Per-parcel boundary information: the nearest parcel that belongs to another plate.
#[derive(Clone, Copy, Debug)]
pub struct BInfo {
    pub other: u32,
    pub oidx: u32,
    /// Unit tangent at this parcel pointing toward the other parcel.
    pub n: V3,
    /// Convergence rate (km/Myr): >0 approaching, <0 separating.
    pub conv: f64,
    pub dist: f64,
}

#[derive(Default, Clone, Debug)]
pub struct Stats {
    pub subducted: usize,
    pub created: usize,
    pub accreted: usize,
    pub absorbed: usize,
    pub cont_lost: usize,
    pub cont_grown: usize,
    pub initiations: usize,
    pub backarcs: usize,
    pub weld_static: usize,
    pub weld_locked: usize,
    pub retired: usize,
    pub dissolved: usize,
    pub deposited: usize,
    pub split_off: usize,
    pub enclaves: usize,
    pub stress_p50: f32,
    pub stress_p95: f32,
    pub stress_max: f32,
    pub rifts: usize,
    pub merges: usize,
    pub mean_v: f64,
    pub max_v: f64,
    pub cont_frac: f64,
    pub n_plates: usize,
    pub n_parcels: usize,
}

#[derive(Clone, Debug)]
pub struct Params {
    pub seed: u64,
    pub n_parcels: usize,
    pub n_plates: usize,
    pub years: f64,
    pub dt: f64,
    pub slice_every: f64,
    pub width: usize,
    pub out: String,
    pub cont_frac: f64,
    pub n_hotspots: usize,
    // --- force balance (Scotese Rule I/II) ---
    pub k_slab: f64,
    pub k_ridge: f64,
    pub k_coll: f64,
    /// Slab suction on the upper plate, as a force per unit trench length (fraction of k_slab).
    pub k_suction: f64,
    /// Push across a young rift (thermal uplift of the rift shoulders), per unit length.
    pub k_rift: f64,
    /// How long after a split the rift push acts (Myr).
    pub rift_push_myr: f64,
    /// Rift tip propagation speed (km/Myr).
    pub rift_prop_v: f64,
    /// How often (Myr) the intraplate stress field is re-evaluated.
    pub stress_every: f64,
    /// Decay length (km) of a boundary pull's influence on interior stress.
    pub stress_l_km: f64,
    /// Fraction of the net (plate-moving) force subtracted from the opposing-pull tension.
    pub stress_beta: f64,
    /// Arc parcels need this much tension (x rift_threshold) to detach for rollback.
    pub backarc_stress: f64,
    /// Oceanic lithosphere is this many times stronger than continental against rifting.
    pub ocean_strength: f64,
    /// Reference load-bearing width (km): stress amplifies by width_ref/width at constrictions.
    pub width_ref_km: f64,
    /// Cap on the constriction amplification.
    pub width_amp_max: f64,
    /// Resistance at a convergent contact that has no subduction zone yet (per unit length).
    pub k_resist: f64,
    /// Oceanic lithosphere this old (Myr) can start subducting after a little compression.
    pub init_age: f64,
    /// Accumulated shortening (km) that forces subduction to start regardless of age.
    pub init_short: f64,
    /// Continental shortening (km) after which a collision belt locks and the trench jumps.
    pub lock_km: f64,
    /// A convergent contact this close (km) to an existing trench starts subducting at once ("virus").
    pub virus_km: f64,
    /// Slab age (Myr) above which the arc can detach and roll back (back-arc opening).
    pub rollback_age: f64,
    /// Probability per Myr that an eligible arc detaches.
    pub rollback_rate: f64,
    /// Width (km) of the arc sliver that detaches.
    pub backarc_km: f64,
    /// Lifetime (Myr) of a detached arc plate before it re-welds (back-arc basins stay small).
    pub backarc_myr: f64,
    /// Speed cap (km/Myr) for a rolling-back arc plate.
    pub rollback_v: f64,
    /// Add tectonically-modulated sub-parcel detail to the rendered elevation.
    pub detail: bool,
    pub drag_ocean: f64,
    pub drag_cont: f64,
    /// Phanerozoic speed limit, km/Myr (Rule II: ~20 cm/yr).
    pub v_max: f64,
    pub omega_relax: f64,
    // --- rifting (Rules III, V, XI): stress-based ---
    /// Tension/strength ratio above which a rift can nucleate.
    pub rift_threshold: f64,
    pub rift_rate: f64,
    /// Overstressed area (steradians) at which the nucleation rate reaches rift_rate.
    pub rift_area_ref: f64,
    /// Plates covering more than this fraction of the surface can fail without a constriction
    /// (broad oceanic breakup / ridge jump), at extra strength.
    pub mega_frac: f64,
    // --- crust evolution ---
    pub arc_rate: f64,
    pub hot_rate: f64,
    pub erosion_tau: f64,
    pub volc_tau: f64,
    pub thick_coeff: f64,
    pub thin_coeff: f64,
    /// Probability per Myr that an active continental-arc parcel converts an adjacent forearc parcel.
    pub accrete_rate: f64,
}

impl Params {
    pub fn default() -> Params {
        Params {
            seed: 42, n_parcels: 40_000, n_plates: 12, years: 1000.0, dt: 1.0, slice_every: 10.0,
            width: 1024, out: "out/run".into(), cont_frac: 0.30, n_hotspots: 12,
            k_slab: 0.015, k_ridge: 0.004, k_coll: 0.2, k_suction: 0.003, k_rift: 0.012, rift_push_myr: 60.0, rift_prop_v: 150.0, stress_every: 2.0, stress_l_km: 4000.0, stress_beta: 0.6, backarc_stress: 0.02, ocean_strength: 2.0, width_ref_km: 1500.0, width_amp_max: 6.0, k_resist: 0.05, init_age: 60.0, init_short: 150.0, lock_km: 500.0, virus_km: 600.0, rollback_age: 70.0, rollback_rate: 0.03, backarc_km: 300.0, backarc_myr: 40.0, rollback_v: 30.0, detail: true, drag_ocean: 1.0, drag_cont: 3.0,
            v_max: 200.0, omega_relax: 0.5,
            rift_threshold: 20.0, rift_rate: 0.025, rift_area_ref: 0.005, mega_frac: 0.25,
            arc_rate: 0.04, hot_rate: 2.0, erosion_tau: 40.0, volc_tau: 60.0,
            thick_coeff: 0.015, thin_coeff: 0.011, accrete_rate: 0.01,
        }
    }
    pub fn from_args(args: Vec<String>) -> Params {
        let mut p = Params::default();
        let mut i = 0;
        while i < args.len() {
            let key = args[i].trim_start_matches("--").to_string();
            if key == "help" || key == "h" {
                println!("{}", HELP);
                std::process::exit(0);
            }
            let val = match args.get(i + 1) {
                Some(v) => v.clone(),
                None => { eprintln!("missing value for --{}", key); std::process::exit(2); }
            };
            let f = || val.parse::<f64>().unwrap_or_else(|_| { eprintln!("bad number for --{}: {}", key, val); std::process::exit(2) });
            match key.as_str() {
                "seed" => p.seed = f() as u64,
                "parcels" => p.n_parcels = f() as usize,
                "plates" => p.n_plates = f() as usize,
                "years" => p.years = f(),
                "dt" => p.dt = f(),
                "slice" => p.slice_every = f(),
                "width" => p.width = f() as usize,
                "out" => p.out = val.clone(),
                "cont-frac" => p.cont_frac = f(),
                "hotspots" => p.n_hotspots = f() as usize,
                "k-slab" => p.k_slab = f(),
                "k-ridge" => p.k_ridge = f(),
                "k-coll" => p.k_coll = f(),
                "k-suction" => p.k_suction = f(),
                "k-rift" => p.k_rift = f(),
                "rift-push-myr" => p.rift_push_myr = f(),
                "rift-prop-v" => p.rift_prop_v = f(),
                "stress-every" => p.stress_every = f(),
                "stress-l-km" => p.stress_l_km = f(),
                "stress-beta" => p.stress_beta = f(),
                "backarc-stress" => p.backarc_stress = f(),
                "ocean-strength" => p.ocean_strength = f(),
                "width-ref-km" => p.width_ref_km = f(),
                "width-amp-max" => p.width_amp_max = f(),
                "k-resist" => p.k_resist = f(),
                "init-age" => p.init_age = f(),
                "init-short" => p.init_short = f(),
                "lock-km" => p.lock_km = f(),
                "virus-km" => p.virus_km = f(),
                "rollback-age" => p.rollback_age = f(),
                "rollback-rate" => p.rollback_rate = f(),
                "backarc-km" => p.backarc_km = f(),
                "backarc-myr" => p.backarc_myr = f(),
                "rollback-v" => p.rollback_v = f(),
                "detail" => p.detail = val != "0" && val != "false",
                "drag-ocean" => p.drag_ocean = f(),
                "drag-cont" => p.drag_cont = f(),
                "v-max" => p.v_max = f(),
                "omega-relax" => p.omega_relax = f(),
                "rift-threshold" => p.rift_threshold = f(),
                "rift-rate" => p.rift_rate = f(),
                "rift-area-ref" => p.rift_area_ref = f(),
                "mega-frac" => p.mega_frac = f(),
                "arc-rate" => p.arc_rate = f(),
                "hot-rate" => p.hot_rate = f(),
                "erosion-tau" => p.erosion_tau = f(),
                "volc-tau" => p.volc_tau = f(),
                "thick-coeff" => p.thick_coeff = f(),
                "thin-coeff" => p.thin_coeff = f(),
                "accrete-rate" => p.accrete_rate = f(),
                _ => { eprintln!("unknown option --{}\n{}", key, HELP); std::process::exit(2); }
            }
            i += 2;
        }
        p
    }
    pub fn to_json(&self) -> String {
        format!(
            "{{\n  \"seed\": {}, \"parcels\": {}, \"plates\": {}, \"years\": {}, \"dt\": {}, \"slice\": {}, \"width\": {},\n  \"cont_frac\": {}, \"hotspots\": {},\n  \"k_slab\": {}, \"k_ridge\": {}, \"k_coll\": {}, \"k_suction\": {}, \"k_rift\": {}, \"rift_push_myr\": {}, \"drag_ocean\": {}, \"drag_cont\": {}, \"v_max\": {}, \"omega_relax\": {},\n  \"rift_threshold\": {}, \"rift_rate\": {}, \"arc_rate\": {}, \"hot_rate\": {}, \"erosion_tau\": {}, \"volc_tau\": {}, \"thick_coeff\": {}, \"thin_coeff\": {}, \"accrete_rate\": {}\n}}\n",
            self.seed, self.n_parcels, self.n_plates, self.years, self.dt, self.slice_every, self.width,
            self.cont_frac, self.n_hotspots, self.k_slab, self.k_ridge, self.k_coll, self.k_suction, self.k_rift, self.rift_push_myr, self.drag_ocean, self.drag_cont,
            self.v_max, self.omega_relax, self.rift_threshold, self.rift_rate, self.arc_rate, self.hot_rate,
            self.erosion_tau, self.volc_tau, self.thick_coeff, self.thin_coeff, self.accrete_rate)
    }
}

const HELP: &str = "tectonic - time-evolved spherical plate tectonics (Scotese rules of thumb)
  --seed N            RNG seed (42)
  --parcels N         crust parcels on the sphere (40000 ~ 113 km spacing)
  --plates N          initial plate count (12)
  --years MYR         simulated time (1000)
  --dt MYR            time step (1)
  --slice MYR         write a time slice every N Myr (10)
  --width PX          equirectangular output width (1024)
  --out DIR           output directory (out/run)
  --cont-frac F       initial continental area fraction (0.30)
  --hotspots N        fixed mantle plumes (24)
  physics: --k-slab --k-ridge --k-coll --k-suction --k-rift --rift-push-myr --rift-prop-v
           --k-resist --init-age --init-short --virus-km --rollback-age --rollback-rate
           --backarc-km --backarc-myr --rollback-v --detail 0|1 --drag-ocean --drag-cont --v-max --omega-relax
           --rift-threshold --rift-rate --arc-rate --hot-rate --erosion-tau --volc-tau
           --thick-coeff --thin-coeff --accrete-rate";

pub struct World {
    pub p: Params,
    pub t: f64,
    /// Mean parcel spacing (radians, ~chord).
    pub s: f64,
    /// Parcel-count scale relative to the 40 000-parcel reference (count thresholds scale with this).
    pub n_scale: f64,
    /// Continent-continent contact parcels per plate pair (from the last census), for shortening bookkeeping.
    pub pair_ncc: HashMap<(u32, u32), usize>,
    /// Cumulative finite rotation of each plate since t = 0 (absolute / hotspot frame).
    pub rot: Vec<[[f64; 3]; 3]>,
    /// Rotation samples at slice times: (t, per-plate rotation if alive).
    pub rot_hist: Vec<(f64, Vec<Option<[[f64; 3]; 3]>>)>,
    pub parcels: Vec<Parcel>,
    pub plates: Vec<Plate>,
    /// Fixed fine grid used to detect gaps (divergent boundaries) where new crust must be created.
    pub grid: Vec<V3>,
    pub hotspots: Vec<V3>,
    /// For each plate pair (min,max): which plate subducts. Fixed once established (Rule IV).
    pub polarity: HashMap<(u32, u32), u32>,
    /// Collisional shortening (km, averaged along the contact) per plate pair since contact.
    pub pair_absorbed: HashMap<(u32, u32), f64>,
    /// Fraction of each plate pair's contact that is continent-continent (from the last census).
    pub pair_ccf: HashMap<(u32, u32), f64>,
    /// Plate pairs created by a rift, with the rift time: the young rift is pushed apart for `rift_push_myr`.
    pub rift_pairs: HashMap<(u32, u32), f64>,
    /// How long (Myr) each plate pair's contact has had essentially no relative motion.
    pub static_myr: HashMap<(u32, u32), f64>,
    /// Shortening (km) taken up at convergent contacts that have not started subducting yet.
    pub pair_compress: HashMap<(u32, u32), f64>,
    /// Boundary tractions coarsened to ~250 km cells, per plate: (cell centroid, net force), rebuilt
    /// each step in forces.rs and consumed by the stress-field evaluation.
    pub cell_tractions: HashMap<u32, Vec<(V3, V3)>>,
    /// Last time the stress field was evaluated.
    pub stress_eval_t: f64,
    /// Last time plate connectivity was enforced.
    pub conn_t: f64,
    /// Rifts currently propagating.
    pub rifts: Vec<ActiveRift>,
    /// Detached arc plates: arc plate -> (parent upper plate, lower plate, detachment time).
    pub arc_plates: HashMap<u32, (u32, u32, f64)>,
    /// Slow drift of each plume (rad/Myr tangent vector), Rule X: "some aren't fixed".
    pub hot_drift: Vec<V3>,
    /// Sub-parcel detail noise octaves for rendering (wavelengths ~1000, 300, 100 km).
    pub detail_noise: Vec<Noise>,
    /// Ocean volume at t = 0 (mean water column per parcel, m) and the current eustatic sea level (m).
    pub sea_v0: Option<f64>,
    pub sea_level: f64,
    /// Sediment reservoir: crustal volume (km thickness x steradian area) eroded off high ground,
    /// waiting to prograde continental margins.
    pub sediment: f64,
    pub rng: ChaCha8Rng,
    pub hash: SpatialHash,
    pub binfo: Vec<Option<BInfo>>,
    pub stats: Stats,
    pub weak_noise: Noise,
}

impl World {
    pub fn new(p: Params) -> World {
        let mut rng = ChaCha8Rng::seed_from_u64(p.seed);
        let n = p.n_parcels;
        let s = (4.0 * std::f64::consts::PI / n as f64).sqrt();
        let positions = fibonacci_sphere(n);
        let grid = fibonacci_sphere(4 * n);

        // --- plate seeds by farthest-point sampling, plate regions = noisy Voronoi ---
        let mut seeds: Vec<V3> = vec![random_unit(&mut rng)];
        while seeds.len() < p.n_plates {
            let mut best: (Option<V3>, f64) = (None, -1.0);
            for _ in 0..256 {
                let c = random_unit(&mut rng);
                let d = seeds.iter().map(|&q| angle(c, q)).fold(f64::MAX, f64::min);
                if d > best.1 { best = (Some(c), d); }
            }
            seeds.push(best.0.unwrap());
        }
        let plate_noise: Vec<Noise> = (0..p.n_plates).map(|_| Noise::new(&mut rng, 5, 2.0, 6.0)).collect();

        // --- continents: warped caps (cratons) covering ~cont_frac of the sphere ---
        let n_crat = ((p.n_plates as f64 * 0.6).round() as usize).max(3);
        let mut cratons: Vec<(V3, f64, Noise)> = vec![];
        let total_area = p.cont_frac * 4.0 * std::f64::consts::PI;
        for _ in 0..n_crat {
            let area = total_area / n_crat as f64 * rng.gen_range(0.5..1.5);
            let r = (1.0 - area / (2.0 * std::f64::consts::PI)).clamp(-1.0, 1.0).acos();
            cratons.push((random_unit(&mut rng), r, Noise::new(&mut rng, 6, 3.0, 9.0)));
        }
        let thick_noise = Noise::new(&mut rng, 8, 6.0, 20.0);
        let age_noise = Noise::new(&mut rng, 6, 2.0, 6.0);

        let mut parcels = Vec::with_capacity(n * 2);
        for &pos in &positions {
            let mut best = (0u32, f64::MAX);
            for (i, &sd) in seeds.iter().enumerate() {
                let c = angle(pos, sd) * (1.0 + 0.3 * plate_noise[i].eval(pos));
                if c < best.1 { best = (i as u32, c); }
            }
            let mut kind = Kind::Oceanic;
            let mut thick = 7.0;
            for (c, r, nz) in &cratons {
                let d = angle(pos, *c);
                let rr = r * (1.0 + 0.4 * nz.eval(pos));
                if d < rr {
                    kind = Kind::Continental;
                    let core = (1.0 - d / rr).clamp(0.0, 1.0);
                    let t = 33.0 + 7.0 * core + 1.5 * thick_noise.eval(pos);
                    if t > thick { thick = t; }
                }
            }
            let birth = if kind == Kind::Continental { -3000.0 } else { -(70.0 + 60.0 * age_noise.eval(pos)).clamp(2.0, 160.0) };
            parcels.push(Parcel { pos, plate: best.0, kind, birth, thick, volc: 0.0, trench_t: NEVER, suture_t: NEVER, hot_t: NEVER, arc_t: NEVER, rift_t: NEVER, stress: 0.0, amp: 1.0, trench_w: 0.0, alive: true });
        }

        let plates = (0..p.n_plates).map(|_| {
            let axis = random_unit(&mut rng);
            let v = rng.gen_range(10.0..40.0);
            Plate { omega: scale(axis, v / R_KM), alive: true, tension: 0.0, n: 0, n_cont: 0, n_weak: 0, mean_v: v, slab: 0.0, suction: 0.0, born: 0.0 }
        }).collect();
        // Plumes: keep them off the poles (equirectangular maps smear a polar island into a band)
        // and give each a slow random drift of ~2 km/Myr.
        let mut hotspots: Vec<V3> = vec![];
        while hotspots.len() < p.n_hotspots {
            let h = random_unit(&mut rng);
            if h[1].abs() < 0.94 { hotspots.push(h); }
        }
        let hot_drift: Vec<V3> = hotspots.iter().map(|&h| { let d = normalize(cross(h, random_unit(&mut rng))); scale(d, 2.0 / R_KM) }).collect();
        let weak_noise = Noise::new(&mut rng, 8, 6.0, 18.0);
        let detail_noise = vec![Noise::new(&mut rng, 10, 15.0, 40.0), Noise::new(&mut rng, 12, 60.0, 150.0), Noise::new(&mut rng, 14, 250.0, 500.0)];

        let n_scale = n as f64 / 40_000.0;
        let n_plates = p.n_plates;
        let mut w = World {
            hash: SpatialHash::new(1.5 * s), p, t: 0.0, s, n_scale, pair_ncc: HashMap::new(),
            rot: vec![IDENT; n_plates], rot_hist: vec![], parcels, plates, grid, hotspots,
            polarity: HashMap::new(), pair_absorbed: HashMap::new(), pair_ccf: HashMap::new(), rift_pairs: HashMap::new(), static_myr: HashMap::new(), pair_compress: HashMap::new(), cell_tractions: HashMap::new(), stress_eval_t: -1.0e9, conn_t: 0.0, rifts: vec![], arc_plates: HashMap::new(), hot_drift, detail_noise, sea_v0: None, sea_level: 0.0, sediment: 0.0, rng, binfo: vec![], stats: Stats::default(), weak_noise,
        };
        w.rebuild_hash();
        crate::step::plate_stats(&mut w);
        w
    }

    pub fn rebuild_hash(&mut self) {
        let it = self.parcels.iter().enumerate().filter(|(_, pc)| pc.alive).map(|(i, pc)| (i as u32, pc.pos));
        self.hash.build(it);
    }

    pub fn alive_plates(&self) -> usize { self.plates.iter().filter(|p| p.alive).count() }

    /// Kilometres on the surface as an angular (chord) distance.
    #[inline] pub fn km(&self, km: f64) -> f64 { km / R_KM }
    /// A length that is at least `mult` spacings and at least `km` kilometres.
    #[inline] pub fn reach(&self, mult: f64, km: f64) -> f64 { (mult * self.s).max(self.km(km)) }
    /// Record the current plate rotations (called at slice times).
    pub fn record_rotations(&mut self) {
        let row: Vec<Option<[[f64; 3]; 3]>> = self.plates.iter().enumerate().map(|(i, pl)| if pl.alive { Some(self.rot[i]) } else { None }).collect();
        self.rot_hist.push((self.t, row));
    }

    /// Cheap deterministic hash-based uniform in [0,1) from (seed, time, index); used where a
    /// per-parcel coin flip is needed inside a read-only pass.
    pub fn rng_f64_hash(&self, i: usize) -> f64 {
        let mut x = (i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ (self.t.to_bits().rotate_left(17)) ^ self.p.seed.wrapping_mul(0xD1B5_4A32_D192_ED03);
        x ^= x >> 31; x = x.wrapping_mul(0x7FB5_D329_728E_A185); x ^= x >> 27; x = x.wrapping_mul(0x81DA_DEF4_BC2D_D44D); x ^= x >> 33;
        (x >> 11) as f64 / (1u64 << 53) as f64
    }
}
