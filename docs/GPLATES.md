# GPlates export

Every run writes `<out>/gplates/`. The mapping from the engine to GPlates is unusually clean
because the engine already works the way GPlates thinks:

| Engine | GPlates |
|---|---|
| Simulation frame = fixed hotspots (Scotese Rule X) | Anchored plate **000** (absolute / hotspot reference frame) |
| Plate = rigid rotation with an Euler vector ω(t), integrated every step into a finite rotation R(t) | **Total reconstruction pole** per plate: `Q(age) = R(t) · R(present)⁻¹`, written as (lat, lon, angle) |
| Plate index `i`, born at rift / dies at weld or consumption | Plate id `100 + i`, rotation sequence covering its lifetime |
| End of run | Present day, 0 Ma; ages are `years − t` |
| Continental parcels at the end | Present-day multipoint features with `reconstructionPlateId` and a valid-time range |
| Boundary parcels with kinematic class | OGR-GMT point file with plate id + class (1 trench, 2 arc, 3 collision, 4 ridge, 5 rift, 6 transform) |
| Slice images | Time-dependent raster sequences (`<layer>-<age>.png`) |

## Files

- `rotations.rot` — one line per plate per slice time: `moving age lat lon angle fixed`. Ages increase within each sequence as GPlates requires. Plates that died before the end are referenced to their last position (identity at their death age).
- `plates.csv` — plate id, index, birth/death in Myr and Ma, alive-at-end flag, final parcel counts.
- `continents_present.gpml` — present-day continental crust, one multipoint feature per living plate, valid from the plate's birth age to 0 Ma. Load together with `rotations.rot` and drag the time slider.
- `boundaries_0Ma.gmt` — present-day boundary parcels as points with attributes `plate_id | class | class_name`.
- `rasters/<layer>/<layer>-<age>.png` — `elev`, `plates`, `age`, `bounds`. In GPlates: *File → Import → Import Time-Dependent Raster*, pick a folder, set the extent to −180…180 / −90…90.
- `README.txt` — the same instructions, next to the data.

## Conventions and checks

- The engine frame is `x, y (pole), z` with `lon = atan2(z, x)`; that is left-handed relative to the geographic frame, so exported rotation **angles are sign-flipped**. Positions (lat = asin y, lon = atan2(z, x)) are unaffected.
- Sanity check in GPlates: load `rotations.rot` + `continents_present.gpml`, set the time to a slice age, and compare with `rasters/plates/plates-<age>.png` loaded as a raster at that age. The multipoints should sit on the coloured continents.

## What is not exported yet (and why)

1. **Terrane histories.** A parcel that changed plate (continental block accreted after a collision, docking) is exported under its *final* plate only; reconstructed back past the accretion it will move with the wrong plate. GPlates' answer is a plate id per block whose rotation sequence changes *fixed plate* at the crossover time (block fixed to A until accretion, then fixed to B, zero relative rotation). The engine has the information (parcel plate ids per slice); the export needs connected-block tracking across slices and crossover lines in the `.rot` file.
2. **Topological (continuously closing) plate polygons.** GPlates' `TopologicalClosedPlateBoundary` features reference boundary line features with plate ids and are resolved at every time step. The engine's boundaries are parcel contacts, not lines; producing gap-free polygons needs boundary tracing (ordering the classified boundary parcels into polylines, splitting at triple junctions). Doable, but it is a separate piece of work; the `bounds` raster and the boundary point file are the interim products.
3. **Oceanic crust as features.** The age raster carries it; isochrons (age contours on the parcel field) would be the GPlates-native form.
4. **netCDF grids.** GPlates prefers netCDF age grids and paleo-DEMs. The raw `elev_*_f32le.raw` slices have everything needed; writing netCDF without a library is possible (CDF-1 is a simple format) but not done yet. For ROCKE-3D boundary conditions this is the same adapter, so it is worth doing once, properly.
