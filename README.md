# SideStep

SideStep is a time-evolved plate-tectonics engine on a sphere. Crust is carried as Lagrangian parcels on
rigid plates that move by a Scotese-style force balance (slab pull, slab suction, ridge push,
collision resistance, mantle drag). Crust is created at gaps, consumed at trenches, thickened
in collisions, accreted, rifted, and volcanised by fixed plumes — and every parcel remembers
its history, so the relief at any time is *derived* from how it got there rather than painted on.

Physics rules follow C. R. Scotese, *Plate Tectonics: Rules of Thumb* — see
[docs/RULES.md](docs/RULES.md) for the rule-by-rule mapping, the crust budget, and the list
of rules that were tried and abandoned.

## Build

Any recent Rust toolchain (developed on stable 1.98, `x86_64-pc-windows-gnu`; no MSVC needed).
Dependencies are `image`, `rand`, `rand_chacha`, `rayon` only.

```bash
cargo build --release
```

## Run

```bash
./target/release/tectonic --years 1000 --slice 10 --out out/run1
```

Defaults: 40 000 parcels (~113 km spacing), 12 plates, 30 % continental area, 12 drifting hotspots,
1 Myr steps, a slice every 10 Myr, 1024×512 equirectangular output. A 1000 Myr run takes
~90 s on 12 cores. `--help` lists every physics knob.

## Globe viewer

Interactive 3D globe of a run: rotate, zoom, scrub or play through time, switch layers.

```bash
./target/release/viewer out/run1 8077
```

Then open <http://127.0.0.1:8077/>. A **Globe | Map** toggle (or `m`) switches between the 3D
globe and a flat equirectangular projection — in map mode, drag pans, the wheel zooms about the
cursor, longitude wraps, a 30° graticule overlays the layer, and the cursor's lat/lon is shown.
On the globe: drag to rotate, wheel to zoom. Either way ←/→ step one slice, space
plays, click the timeline to jump; ▲/▼ marks on the timeline are rift / suture events. The
viewer is a single self-contained page ([viewer/index.html](viewer/index.html), WebGL, no
external libraries) embedded in the `viewer` binary; it serves the run directory read-only on
localhost. On Windows, `run-viewer.cmd [run_dir] [port]` starts it and opens the browser.

Contact sheet of a run's history (every 100 Myr, 3 columns, `elev` / `plates` / `age` layer):

```bash
./target/release/montage out/run1 100 3 elev
```

Other seeds / sizes:

```bash
./target/release/tectonic --seed 7 --plates 16 --cont-frac 0.35 --out out/seed7
```

```bash
./target/release/tectonic --parcels 160000 --width 4096 --slice 50 --out out/hires
```

## Output

`out/<run>/`
- `params.json` — the parameters used
- `log.csv` — per step: plates, parcels, continental fraction, mean/max speed, parcels subducted / created / accreted / absorbed, continental parcels lost / grown, rift and suture events
- `tNNNNN/` per slice:
  - `elev.png` — hypsometric elevation
  - `elev_<W>x<H>_f32le.raw` — raw little-endian float32 elevation in metres, row-major, lon −180→180 left→right, lat 90→−90 top→bottom. This is the hand-off for regridding to a GCM grid.
  - `plates.png` — plate IDs (continents bright, oceans dark, boundaries black)
  - `age.png` — oceanic crust age (bright = young), continents brown
  - `stress.png` — intraplate tension (hot colours = corridors where opposing boundary pulls stretch a plate)
  - `bounds.png` — classified plate boundaries on a muted relief: red trench (subducting side) / amber arc (upper plate), magenta continental collision, green oceanic ridge, cyan continental rift, white transform; purple = suture younger than 100 Myr, yellow = active hotspot. Pixel counts per class are in `meta.json`.
  - `meta.json`
- `montage_<layer>.png` — from the `montage` tool
- `gplates/` — rotation file, plate lifetimes, present-day continental multipoints (GPML), boundary points (GMT), time-named raster sequences. See [docs/GPLATES.md](docs/GPLATES.md).

## Resolution

All physics rules are expressed in kilometres or scaled by the parcel count, so resolution is a
free choice. Double the linear resolution (56 km spacing) with a half-Myr step so plates still move
less than a spacing per step:

```bash
./target/release/tectonic --parcels 160000 --width 2048 --dt 0.5 --out out/hires
```

This takes ~15× longer than the 40K default (about 25 minutes for 1000 Myr on 12 cores).
Three times the default linear resolution (37 km spacing, 360K parcels, 3072-px maps) is
about 90 minutes; four times (28 km, 640K parcels, 4096-px maps) about 9 hours:

```bash
./target/release/tectonic --parcels 360000 --width 3072 --dt 0.33 --out out/hr360
```

```bash
./target/release/tectonic --parcels 640000 --width 4096 --dt 0.25 --out out/hr640
```

Before trusting a rule change, check that event rates do not depend on the step or the parcel
count — run a couple of seeds at `--dt 1` and `--dt 0.25` and compare rifts / initiations /
back-arcs / mean speed in `log.csv` (they should agree to within seed-to-seed scatter):

```bash
for sd in 3 7; do for dt in 1 0.25; do ./target/release/tectonic --seed $sd --dt $dt --slice 50 --out out/sweep_s${sd}_dt${dt}; done; done
```

The physical resolution is the parcel spacing (113 km at 40K, 56 km at 160K), not the pixel
count: the renderer samples each pixel with a Gaussian of half a spacing, so features smaller
than ~2 spacings do not exist. Pick the map width so there are 2–3 pixels per parcel at the
equator (1024 px ≈ 39 km/px, 2048 px ≈ 20 km/px).

## Layout

- `src/geom.rs` — sphere math, Fibonacci sampling, band-limited noise, spatial hash
- `src/world.rs` — parcels, plates, parameters, initial conditions
- `src/step.rs` — one time step: advect → detect boundaries → subduct / collide / accrete → fill gaps → crust evolution → forces → rifting → welding / cleanup
- `src/forces.rs` — torque / drag force balance per plate
- `src/render.rs` — history → elevation, rasterisation, PNG / raw writers
- `src/bin/montage.rs` — contact-sheet tool
- `src/bin/viewer.rs` + `viewer/index.html` — localhost server + WebGL globe viewer
- `src/gplates.rs` — GPlates export (rotations, GPML, GMT, raster sequences)

## Status (v0.3)

Runs are reproducible (same seed ⇒ identical log) and the dynamics are step- and
resolution-independent to within chaotic scatter (see the sweep recipe above). Latest
production run: 640K parcels / 28 km / 4096 px, 1000 Myr in ~16 h — 4–12 plates (ending in a
supercontinent phase), 10 rifts, 122 subduction initiations, 15 back-arc basins, sea level
−270…+100 m, continental fraction 0.29 → 0.27 over the Gyr (the sediment return path covers
about half of collisional loss at this resolution), continents coherent, with collision belts,
arcs, back-arc basins, cratonic highs and textured sea floor.
Subduction initiates only after compression builds (old passive margins, exhausted shortening
budget, or next to an existing trench); rifts nucleate at weaknesses and propagate along them
before splitting a plate; old slabs roll back and open back-arc basins behind detached arcs;
sea level is eustatic; rendered elevation carries tectonically-modulated sub-parcel detail;
rifting is driven by a spatially resolved intraplate stress field (opposing boundary pulls with
relief near open rifts — the `stress` viewer layer shows it) instead of a per-plate clock.
The crust cycle is closed: eroded volume returns as prograding shelf at continental margins,
and continental fraction holds flat over a Gyr at 160K+ parcels (mild upward drift at 40K).
Plate speeds agree across resolutions (boundary forces are normalised by smoothed boundary
length). Known gaps: no fracture zones, arc detachment is probabilistic, and the detail layer
is noise rather than a process model. See
[docs/RULES.md](docs/RULES.md).

## Reference

C. R. Scotese, *Plate Tectonics: Rules of Thumb* (1993, updated 2012). The handout itself is
not redistributed here; the rule-by-rule mapping is in [docs/RULES.md](docs/RULES.md).
