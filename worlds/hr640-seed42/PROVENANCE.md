# World: hr640-seed42

The 1000 Myr, 640K-parcel (28 km) world completed 2026-08-26. Rendered products live in
`out/hr640/` (local, ~5.2 GB, not in git); this folder is the permanent record.

## Exact regeneration

The engine is deterministic: the same code and parameters reproduce this world bit-for-bit,
including every intermediate slice and the GPlates export (~16 h on 12 cores):

```
git checkout world/hr640-seed42
cargo build --release
./target/release/tectonic --seed 42 --parcels 640000 --width 4096 --dt 0.25 --years 1000 --slice 10 --out out/hr640
```

- Generating commit: `edb87b4` (tagged `world/hr640-seed42`); no engine source changed
  between that commit and the tag's creation.
- Full parameter set: [params.json](params.json) (identical to the command line above plus defaults).
- Per-step history: [log.csv](log.csv). Final-slice metadata: [meta.json](meta.json).
- Present-day elevation: [final_elev.png](final_elev.png) (hypsometric, 4096x2048; the raw
  float32 field is `out/hr640/t01000/elev_4096x2048_f32le.raw` when regenerated).

## Headline numbers

4-12 plates (supercontinent phase at the end), 10 rifts, 122 subduction initiations,
15 back-arc basins, sea level -270..+100 m, continental fraction 0.29 -> 0.267,
mean plate speed ~35 km/Myr.
