# Scotese's "Rules of Thumb" → engine implementation

Source: C. R. Scotese, *Plate Tectonics: Rules of Thumb* (12/19/93, updated 08/13/12) —
`RULESOFTHUMB081312copy.pdf` in the repo root. This records how each rule is realised (or
deliberately not yet) in the engine as of v0. File references are to `src/`.

| Rule | Statement (short) | Implementation | Status |
|---|---|---|---|
| **I** | Plates are pushed or pulled, not dragged; mantle is passive; continental keels slow plates; ridge-surrounded plates are slow; plates with neither ridge nor slab don't move. | `forces.rs`: quasi-static balance `D·ω = τ` every step. Torque from boundary parcels only; drag tensor `D = Σ c·(I − r rᵀ)` over all parcels, `c = drag_cont = 3×` oceanic. No inertia, no mantle-flow term. A plate with no boundary forces stops; encircling ridges cancel by symmetry. | done |
| **II** | Slab pull ≫ ridge push (80/20); fast plates have slabs; speed limit ~20 cm/yr; subduction implies spreading. | `k_slab : k_ridge ≈ 4 : 1`; slab pull ∝ slab age (`age/80 Myr`, clamped 0.2–1). `v_max = 200 km/Myr` hard cap. Spreading follows from area conservation (gaps are filled). Plus **slab suction** on the upper plate at 20 % of slab pull — not in the PDF, but the standard force that lets a trench-ringed supercontinent pull itself apart; without it the model locked into 2–3 plates. | done |
| **III** | Ridges are passive: crust breaks where pulled; continents break first; ridges align parallel to trenches; fracture zones point to trench. | `step.rs::fill_gaps`: zero-age oceanic crust appears wherever parcels separate beyond 0.9 spacing, owned by the plate with the local majority; nothing drives a ridge. Continental breakup is the explicit rifting rule (V). Rift lines are noisy great circles ⟂ plate velocity (⇒ parallel to the pulling trench). Fracture zones are not modelled. | done (fracture zones: no) |
| **IV** | Subduction is forever; only continent–continent collision stops it; hard to start; trench jumps over small colliding terranes. | `polarity[(a,b)]` set on first convergent overlap, never flipped. Continent arriving at a continental upper plate: **absorbed** (parcel removed, its crustal volume thickens the belt across ~2.5 spacings — volume-conserving shortening) for the first 25 parcel-rows; after that the belt is locked and the **whole connected continental block jumps to the upper plate in one event**, subduction resuming behind it. Plates weld only when their shared contact is >60 % continental and static. Initiation is instantaneous on first contact (no passive-margin-age rule yet). | done except initiation difficulty |
| **V** | Pacific vs Tethyan systems; ridge in a one-sided ocean gets subducted; the pulled continent then rifts, India-style tear or Pangea-style breakup; sutures/hotspots make it likelier. | Emergent + rule: continental plates accumulate `tension` at a base rate (÷220 Myr) × (1 if slab-pulled, 0.35 if not) × √size; threshold lowered by the fraction of weak parcels (recent sutures, recent hotspots). A rift splits the plate along a noisy great circle; for the next 60 Myr the young rift is pushed apart (`k_rift`, thermal uplift of the rift shoulders — Rule X's "little push") until a true ridge exists and ordinary ridge push takes over; the pair cannot re-weld for 80 Myr. Without the sustained push most rifts stalled and were inverted into fold belts. Ridge subduction is emergent (a plate consumed below 30 parcels dissolves). | done |
| **VI** | Plates subduct normally; oblique convergence → terranes, highest Andean peaks. | Only the trench-normal relative velocity drives consumption and slab pull. Transform/near-static continental contacts consume nothing. Oblique topographic effects: no. | partial |
| **VII** | Andean vs W-Pacific style set by absolute motions. | Emergent from the upper plate's own motion; no explicit back-arc opening. | emergent |
| **VIII** | Island arcs don't ride across oceans; always a continent nearby. | Thin (<33 km) continental parcels with no normal-thickness continental neighbour are **subducted** at intra-oceanic trenches instead of docking; thicker or attached crust docks (terrane accretion) and obducts the oceanic crust it overrides. | done |
| **IX** | Rollback makes odd trapped intracontinental basins. | Trapped ocean appears only when plates weld across a partly-oceanic contact; no rollback. | partial |
| **X** | Hot spots: deep, roughly fixed, weaken continents, build tracks, not the drivers. | Fixed plumes; parcels within 0.8 spacing get `hot_t`, a volcanic pile (3 km/Myr oceanic, capped 8 km, decays τ = 60 Myr) and a transient 1 km thermal swell. Tracks emerge from plate motion. Recent hotspot = weak parcel for rifting. | done |
| **XI** | Sutures are long-lived weaknesses, future rifts. | `suture_t` recorded at collision contacts and block accretions; 5× weight when choosing the rift nucleus; 500 Myr memory. | done |
| **Last** | Catastrophic, not chaotic; collision and ridge subduction are the instabilities; both lower sea level. | Collisions lock and jump abruptly; ridge subduction removes a plate and re-balances neighbours in one step. Sea level fixed at 0 — no eustasy. | partial |

## Failed rifts (aulacogens)

Two rules added after the 160K-parcel run showed deep slots criss-crossing continents:

- **Static contacts close.** A plate pair whose contact has had < 4 km/Myr relative motion for 40 Myr (and is older than the 80 Myr rift refractory) is welded: no strain, no boundary. This is how a failed rift stops being a plate boundary.
- **Enclosed ocean floor fills.** An intraplate oceanic parcel older than 15 Myr that is enclosed by continental crust of its own plate (≥ 4 continental neighbours on opposing sides — short resultant of the neighbour directions — and at least half of all neighbours) becomes 28 km continental crust: a sediment-filled basin ~900 m below sea level that then thickens by diffusion. The "opposing sides" test matters: with plain adjacency the rule converted every passive margin one row per step and ate the oceans.

## Crust budget (not in the PDF, needed for a 1 Gyr run)

- **Continental growth:** island arcs whose volcanic pile reaches 5 km convert, as contiguous segments of ≥3 parcels, into 30–38 km continental crust; continental arcs grow a forearc parcel seaward at 1 %/Myr. Together ≈ 500–1200 parcels per Gyr.
- **Continental loss:** collisional shortening (≈1700–3400 parcels/Gyr, volume kept as thickness) and subducted thin slivers (≈300–1300). Net continental area drifts −10 to −20 % per Gyr depending on seed. Earth's is ~constant; closing this gap (post-orogenic extension that adds area without making confetti) is open work.
- **Erosion:** diffusion of crustal thickness between continental neighbours (τ = 25 Myr per-neighbour relaxation; volume-conserving) plus slow loss of crust above 35 km (τ = 100 Myr).

## Elevation model (history → relief)

| Crust | Elevation |
|---|---|
| Oceanic | −min(2600 + 350·√age, 5700) m; −3000 m flexural deepening fading over 10 Myr next to a trench on the lower plate; + volcanic pile; + 1000 m thermal swell decaying after a plume passage. |
| Continental | 180 m per km of crust above 32.8 km (Airy). Initial 33–40 km (thicker cratonic cores); collisions add the absorbed parcel's full thickness to the belt; arcs add 0.04 km/Myr; rifted margins thin at `thin_coeff`×divergence (floor 20 km ⇒ shelf); cap 70 km. |

## Rules tried and abandoned (do not reintroduce)

- Parcel-by-parcel transfer of a colliding continent to the upper plate → plates interleave through the continent, later boundaries cut through the mixture, continents shred.
- Converting oceanic parcels next to over-thick crust into thin "shelf" crust (plateau collapse) → manufactures wide submerged confetti.
- Rift tension scaled by plate speed → stagnant supercontinents never break.
- Welding plates after a shortening budget regardless of contact type → the whole planet welds into one plate.
- Basin fill by neighbour *count* alone (no enclosure test) → a conversion front creeps out from every passive margin and turns the oceans continental.
- A one-step velocity "kick" at rifting → decays in ~5 Myr, rifts stall and invert; replaced by a 60 Myr rift push.

## Known gaps to revisit

- Subduction initiation difficulty; back-arc extension/rollback; fracture zones; eustatic sea level.
- Rift lines follow a noisy great circle rather than tracing sutures as curves.
- Continental area budget (see above).
- Relief is coarse (113 km parcels, diffused): an amplification stage is needed for Orogen-level detail.
