# World: hr640v2-seed42

Second production world (2026-08-26): 1000 Myr, 640K parcels (28 km), 4096 px, with the
spatially resolved intraplate stress field (v0.4 physics). Same seed as hr640-seed42, so the
differences between the two worlds are purely the physics changes.

## Exact regeneration (~10 h on 12 cores)

```
git checkout world/hr640v2-seed42
cargo build --release
./target/release/tectonic --seed 42 --parcels 640000 --width 4096 --dt 0.25 --years 1000 --slice 10 --out out/hr640v2
```

Generating commit: `08df7d2` (tagged `world/hr640v2-seed42`).

## Headline numbers

5-18 plates with full Wilson cycling (supercontinent minimum at 400-500 Myr, recovery to 18 by
800 Myr), 39 rifts, 412 subduction initiations, 9 back-arc basins, sea level -550..+485 m
(post-breakup high stand at 900 Myr), continental fraction flat at 0.29-0.30, mean plate
speed ~36 km/Myr. Compare hr640-seed42: 10 rifts, decay to 4 plates, no recovery.
