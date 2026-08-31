# World: hr640v4-seed42

The stress-field production world: 640K parcels (28 km), 4096x2048, 1000 Myr, completed
2026-08-28. Full physics: spatially resolved constriction-aware intraplate stress on all
lithosphere, propagating rifts, neck snapping, back-arc rollback gated on measured arc tension,
eustatic sea level, closed crust cycle, resolution-invariant (area-scaled) boundary forces.

## Exact regeneration (deterministic, ~17 h on 12 cores)

```
git checkout world/hr640v4-seed42
cargo build --release
./target/release/tectonic --seed 42 --parcels 640000 --width 4096 --dt 0.25 --years 1000 --slice 10 --out out/hr640v4
```

## Headline numbers

8-11 plates, 31 rifts, 22 back-arc basins, land 0.30 -> 0.27, area-weighted plate
compactness 0.196 at 1000 Myr (hr640v2, the pre-shape-physics world: 0.124).
Validation ladder: 160K gate passed with stress p50/p95 = 8.4/35.3 vs 40K 6.6/26.1.

## Known open issues carried into this world

- Two dominant plates at the end (48% + 31% of the surface): broad oceanic breakup /
  ridge jumps are still unmodelled, so plates without constrictions grow unchecked.
- Back-arc rate shows a residual resolution slope (40K ~50/Gyr, 160K 27, 640K 22).
