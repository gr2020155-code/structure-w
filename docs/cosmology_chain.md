# The cosmology chain: matrix records in English

27 November 2026 · Theory: Hrachya Danielyan

This file translates the matrix records on cosmology: redshift without expansion, the law of link density, the relic (CMB) spots, the share of mass, and the Hubble tension. It is a reading aid. The authoritative text is the Russian matrix [`system/matrix.json`](../system/matrix.json); if the two differ, the matrix is right.

## How to read

- **IDs.** In the matrix, IDs use Cyrillic letters. Here they are written in Latin: А = A (author's postulate), Д = D (data), Р = R (calculation), В = V (statement), С = S (withdrawn), Х = X (contradiction). V-374 here is В-374 in the matrix.
- **Status** of a statement: *derived* (follows from the named postulates without extra assumptions), *computed* (output of a program whose rules were written and pushed to GitHub before it was run), *candidate* (agrees with data, not proven), *open*, *contradiction*. Postulates, data and calculation records have no status.
- **Data rule (author, 26.11).** Only measured data count as a check. A result of official theory may be used only if it is derived in W; otherwise it is marked "bridge, model". So the Planck value H₀ = 67.4 is not taken as data: it is a fit of a model with expansion.
- **Fit rule (author, 26.11).** A choice first seen in the data is not a fit once it is explained strictly within the structure. It is still disclosed.
- **Calculation rule.** For every calculation R-xxx the rules (`alpha/*_rules_*.md`), the pass criterion and the author's expectation are committed before the program runs. The commit hashes are in each R record.
- Words used for W objects: *link* (a stable sequence of acts of change, A-02); *figure-eight* (a local self-closed loop — matter, A-36); *pre-material links* (links below the level of matter, A-106); *tick* and *step* (one change, A-01); *k* — link density; *u = k/k₀*.

## The chain at a glance

| Step | What | Records | Status |
| --- | --- | --- | --- |
| 1 | Space does not expand; link density grows; the number of points is constant | A-90, A-91, A-92, A-93, A-103 | postulates |
| 2 | Unfolding chain 00 → 10 → 20 → 30 → 3W; no global "now"; density = share of points distinguishable in 3D; the start is finite | A-98 … A-102, V-370, V-371 | postulates, derived |
| 3 | How light moves in a medium of growing density: 1 + z = k₀/k; same kinematics as a scale factor | V-374 | candidate |
| 4 | Power laws and shape-of-space ideas tested and rejected | V-375 … V-386, R-122 … R-131 | computed / open |
| 5 | The law of density from the BAO ruler: g² = A·u⁻³ + B | V-387, R-132 | open (observation) |
| 6 | Postulate: the squared rate is the sum of shares; light share u⁻⁴ | V-388, A-105, V-389, R-133 | postulate, candidate |
| 7 | Mass share = figure-eights (0.040) + pre-material links (0.256) | A-106, V-391, R-134 | computed |
| 8 | Official cosmology as a special case of W (same form, other reading) | V-390 | candidate |
| 9 | Exact spot angle vs local rate: the Hubble tension appears inside W | V-392, R-135, X-13 | contradiction (open) |
| 10 | Tests of the tension: local, candles, Cepheids, early share, sky, spectrum tail, medium, masses, law of N | V-393 … V-402, R-136 … R-143 | computed / candidate |
| 11 | Pre-material links must cluster like matter (full CMB spectrum passes at h = 0.687) | V-398, R-141 | computed |

Resolved on the way: X-11 (supernovae vs a concave law) and X-12 (BAO ruler vs power laws), both resolved by A-105 on 27.11. Open: X-13 (Hubble tension).

## Bridges used (official tools, not derived in W)

| Bridge | Where | Why it is needed |
| --- | --- | --- |
| Factor 3/(8πG) between a share and the squared rate | R-133 and later | To turn the measured T₀ into the light share Ω_r |
| Neutrino factor 1 + 0.2271·3.046 for the light share | R-133, R-135 … R-143 | W has not derived a neutrino background; branch ν0 (photons only) is always computed too |
| Three-level hydrogen atom (Peebles / RECFAST, fudge F = 1.00 and 1.14) | R-135 … R-143 | Recombination redshift z* and drag redshift |
| Scaling of atomic physics with m_e and α | R-142, R-143 | Energies ∝ α²m_e, α_B → α²m_e⁻²·α_B(T/α²m_e), 2s decay ∝ α⁸m_e, σ_T ∝ α²/m_e² |
| Primordial nucleosynthesis fits (Fields et al. 2020): Y_p ∝ η^0.039 (N/3)^0.163, D/H ∝ η^−1.597 (N/3)^0.396 | R-139 | Helium for an extra early share |
| CAMB 2.0.5 (perturbation code) | R-141 | Full TT/TE/EE spectrum; the "smooth links" branch as a fluid with c_s² = 1 |

## Data

| ID | Data | Source (level) |
| --- | --- | --- |
| D-62 | Lunar laser ranging: Ġ/G = (7.1 ± 7.6)·10⁻¹⁴ per year | processed |
| D-64 | Relic temperature far away: T = 5.08 ± 0.10 K at z = 0.89; T(z) = T₀(1+z)^(1−β), β = 0.006 ± 0.013 | processed |
| D-65 | Oldest globular clusters t = 13.5 Gyr; distance ladder H₀ = 73.04 ± 1.04 (SH0ES), 70.39 (CCHP) | processed |
| D-66 | Pantheon+ supernovae Ia: 1701 SNe, corrected magnitudes, z_HD, full covariance | processed |
| D-79 | Relic acoustic angle 100θ* = 1.0411 ± 0.0003; T₀ = 2.72548 ± 0.00057 K | processed |
| D-80 | BAO ruler, DESI DR2: D_V, D_M, D_H over r_d at z = 0.295 … 2.33, full covariance | processed |
| D-81 | Constancy of atomic lines: Δα/α = −0.5 ± 0.5 ± 0.4 ppm at z = 0.6–2.4; abs(Δμ/μ) < 5·10⁻⁶ at z = 2–4 | processed |
| D-82 | Primordial deuterium D/H = (2.527 ± 0.030)·10⁻⁵; 100Ω_bh² = 2.166 (b1, computed rate) or 2.235 (b2, measured rate) | processed |
| D-83 | CCHP TRGB (HST + JWST; 25 SNe, 20 galaxies) and JAGB (8 galaxies), zero point NGC 4258 (arXiv:2408.06153v3) | processed |
| D-84 | Primordial helium Y_p: 0.2453 ± 0.0034 (Aver), 0.2436 ± 0.0040 (Hsyu), 0.2370 ± 0.0034 (EMPRESS); the series disagree | processed |
| D-85 | Planck 2018 binned spectra: TT (83 bins, ℓ = 48–2499), TE and EE (66 bins each, ℓ = 48–1988) | processed |

## Author's postulates

| ID | Postulate |
| --- | --- |
| A-90 | The universe does not expand in the sense claimed by official science. |
| A-91 | The total number of points is constant. |
| A-92 | Space does not expand; the local density of links changes. |
| A-93 | A change of link density: 1 tick — 1 change. |
| A-98 | Unfolding chain: 00 → 0W = 10 (a point closes on itself and appears) → 1W = 20 (one-dimensional space closes into a strictly symmetric figure-eight) → 2W = 30 (two-dimensional space reaches its limit and becomes a circle — the primary atom of 3D space) → 3W, a black hole: absolute indistinguishability, no 3D space. Each level adds one degree of freedom. |
| A-99 | A point is a point locally and globally a field of possible distinctions; points can be distinguishable in 1D, 2D or 3D. Until space closes globally it does not exist. Local figure-eights are already matter. Statics is a contradiction: 1 tick — 1 local change. So there is no "now" in the universe: a "now" would be a static slice. |
| A-100 | Link density (A-92) is the share of points distinguishable in 3D. There are many directions in 3D, but not infinitely many: infinity is a contradiction. Anything not closed does not exist; there is no choice; contradictory trajectories do not exist. |
| A-101 | Symmetry of unfolding (1W = 20): from the start, 1D unfolds left-up and right-down, each closing into a circle; they meet at point "zero" in counter-motion, forming a figure-eight, and +1 degree of freedom (2D) must appear. Zero is not static; it is the place of closure. |
| A-102 | In 3D, points, figure-eights and circles all get a degree of freedom. 3W as a shape is a ball with a marked central place where all figure-eights and circles close; apparently a horn torus. |
| A-103 | "Expansion" of space is the addition of links: their density grows. Light goes 1 step per tick, but in the same tick the medium changes locally by one minimal change. So a second photon after the first passes through a different medium. |
| A-104 | The energy released by an explosion (a light source) per tick is split equally among the three degrees of freedom of 3D space. |
| A-105 | The square of the relative growth rate of link density equals the sum of the shares (V-388 item 4). On 27.11 first accepted as an assumption, then confirmed by the author as a postulate. |
| A-106 | There are no dark substances; there are links. Matter is local figure-eights; there is also a pre-material level of links of different density. |
| A-36 | Mass is the density of self-closed links; continuations into the medium are not mass. |

## Statements

### Kinematics and the start

**V-370 (derived).** What A-99 means for earlier records. Time T in the density law is not a global "now" but the local count of ticks from level 30 along the path of the observer or source. Formulas do not change. The formalisation with a common tick is a simplification of counting, not a claim of a global slice.

**V-371 (derived).** The redshift of the start is finite. A law 1 + z = (T₀/T)^p gives infinite z at the start and an infinite sound path; by A-100 infinity is a contradiction, so such a law is excluded at the very start. Link density at 30 is above zero, and 1 + z_max = k(now, at the receiver)/k(30) is finite.

**V-374 (candidate).** How light moves while link density grows. (1) Places do not move (the number of points is constant, A-91). The number of steps between two places is k·χ: steps are points already distinguishable in 3D, so there are more steps without new points. A photon makes one step per tick, so it passes dχ = dT/k, and links are added on the rest of the path in the same tick (A-103). (2) [mathematics] Hence: the gap between two photons stretches by k(reception)/k(emission), so 1 + z = k₀/k and bursts stretch by exactly 1 + z; all wavelengths stretch by one factor, so a thermal spectrum stays thermal and the radiation temperature was higher by 1 + z (D-64); the path in steps is D = ∫(1 + z) dT. (3) [reading] As kinematics of light this is the same mathematics as a metric with scale factor a = k(T) and cosmic time = local tick count. The difference from official science is the cause (places do not fly apart; links are added) and the law k(T), which in W must follow from the growth of the share of points distinguishable in 3D, not from the equations of gravity.

### Laws tested and rejected (steps 4)

**V-375 (open).** If 3D distinguishability spreads one step per tick from circle 30, the covered region is a tube around the circle, volume 2π²R·r², which at r = R becomes exactly a horn torus (the shape 3W). Then the share s ∝ T². None of six laws (s ∝ T, T², T³ times steps ∝ s or s^(1/3)) passes both supernovae and the relic angle (R-123).

**V-376 (candidate).** Author: link density grows faster at the start, then slower; the decline is gradual, not linear, not two regimes. The tube with steps ∝ s^(1/3) gives k ∝ T^(2/3), a concave growth with no fitted numbers; the relic agrees within 7 %.

**V-377 (candidate; was X-11).** With V-374, any concave growth makes distant supernovae brighter than linear growth, but Pantheon+ shows them fainter (Δχ² ≥ +107; +702 for k ∝ T^(2/3)). Possible causes were luminosity at the place of the explosion, spreading of light on the torus, or data processing. Resolved on 27.11 by A-105.

| Tested idea | Result | Records |
| --- | --- | --- |
| Power laws 1 + z = (T₀/T)^p for the relic | infinite sound path for p ≥ 1; the needed p = 0.64 is a fit, supernovae want 1.41 | R-122 |
| Constant curvature (torus as curvature) | supernovae +135 even with fitted radius; not removed | V-378, R-124 |
| Relic as "reverse" light along the figure-eight | open idea of the author; on the torus surface reverse light through zero spreads infinitely | V-379, V-380, R-125 |
| Light on the horn-torus surface, one fitted place | supernovae +8.7 (threshold 10), relic 2.3 % off; in the full torus volume no place fits both | V-380, V-381, R-125, R-126 |
| Pantheon+ processing (selection bias, free α and β, z_CMB) | still +638 … +781; not the cause | V-382, R-127 |
| Supernova luminosity depending on local density, L ∝ k^q | best q = 0.78, +10.2 (threshold 10); not passed | V-383, R-128 |
| Luminosity linear in ticks | +0.2 with one fitted number; observation only | V-384, R-129 |
| Luminosity = mean over three tube axes | supernovae +3.9, but the BAO ruler fails (χ² = 1420/12): the dimming is geometry, not luminosity | V-385, V-386, R-130, R-131 |
| Any power law k ∝ T^p against the BAO ruler | best p = 0.91, χ² = 97/11; not passed (was X-12) | V-386, R-131 |

### The law of density

**V-387 (open; observation, R-132).** With kinematics V-374 the BAO ruler measures directly the relative rate of density growth per tick g = d ln k/dT = 1/D_H. The data choose g² = A·u⁻³ + B (free power 2.94): ruler χ² = 10.3 on 11. With the same function supernovae pass with no luminosity corrections: Δχ² = +3.9 (threshold 11). In W variables u³ ∝ s, so g² = A·s₀/s + B: early the A/s term dominates and gives k ∝ T^(2/3) (the tube law is right as the initial limit); late B dominates, a constant relative rate. The share A/(A + B) today is 0.30. The relic angle with this law misses (0.0124–0.0209 vs 0.0104): another term growing faster than u⁻³ is needed early.

**V-388 (open; reading).** (1) [mathematics] g² = A·s₀/s + B is equivalent to k̈ = B·k − (A/2)·k₀³/k². (2) Term B: if every existing link realises new distinctions each tick (A-93), the number of new links per tick is ∝ the number of links, g is constant — exponential growth. (3) Term A/s: self-closed links (mass, A-36) make no changes "outward" and slow the growth; their number does not change, while the number of points distinguishable in 3D grows ∝ s, so their share ∝ 1/s — exactly the form the data chose. (4) Squares add, not rates (sum of rates: 23 on 11; sum of squares: 10.3 on 11): the two contributions are independent, like perpendicular components. This was the missing postulate, later accepted as A-105. (5) Light is also a share: the number of light changes is conserved and the energy of each ∝ k, so its share ∝ u⁻⁴ — the early term the relic needs.

**V-389 (candidate; R-133).** Under A-105: g² = g₀²[Ω_r u⁻⁴ + Ω_m u⁻³ + Ω_B]. Ω_m = 0.297 from the ruler; Ω_r from the measured relic temperature T₀ = 2.725 K and the own rate H₀ = 73.5 ± 1.1 (Pantheon+ with Cepheid calibrators); the factor 3/8πG is a bridge. With no fitted numbers, θ = 0.0080–0.0143 (photons only) and 0.0071–0.0129 (with neutrinos) — both contain 0.0104. The check was weak: the width comes from the unknown sound speed. The exact check (R-135, V-392) later failed. Ruler, supernovae and relic agree with one law; X-11 and X-12 are resolved, since A-105 was accepted as a postulate (27.11; the earlier wording "if A-105 becomes a postulate" is kept in `was`).

**V-390 (candidate).** The law V-387–V-389 has the same form as the official Friedmann equation; the reading is different. Official cosmology, where it is right, is a special case of W:

| Term | Official reading | W reading |
| --- | --- | --- |
| Redshift | space itself stretches | space does not expand (A-92); 1 + z is a ratio of link densities (V-374) |
| Constant rate (≈ 70 %) | dark energy / cosmological constant (vacuum estimate off by ~120 orders) | every link makes a change each tick (A-93); growth drives itself; no substance needed |
| u⁻³ (≈ 30 %) | gravity of matter, ~85 % of it dark matter | share of self-closed links and pre-material links (V-391) |
| u⁻⁴ | radiation diluted and stretched | number of light changes conserved, energy of each ∝ k |
| Squares | Einstein's equations | postulate A-105 |
| Start | singularity, inflation | infinity is a contradiction; the start is finite (V-371), unfolding from circle 30 (A-98) |

Item 7 (Hubble tension), corrected on 27.11: officially it is a conflict of the local ladder with a model fit to the relic; in W the measurements do not disagree (V-367), but under A-105 the exact relic angle needs the far rate h ≈ 0.685 while the local rate is 0.735 (V-392, X-13), so the tension appears inside W too. The earlier text ("the relic agrees with 73.5") is kept in `was`.

**V-391 (computed; R-134).** The mass share is local figure-eights plus pre-material links whose number does not change. Figure-eights alone (ordinary matter, share 0.040–0.041 from deuterium, D-82) fail: ruler χ² ≈ 4900 on 12, supernovae +440. Pre-material links with share ∝ u⁻³ (their number after 30 does not change, like figure-eights) pass: ruler 10.3 on 11, supernovae +3.9; their share is 0.256–0.257, 6.2–6.4 times the figure-eights. Shares ∝ u⁻², u⁻¹ or constant fail. So official "dark matter" reads in W as the pre-material level of links (A-106). Open: why their number is constant and why there are ~6.3 times more of them.

### The Hubble tension

**V-367 (derived, by the data rule).** Measured: SH0ES 73.04 ± 1.04 and CCHP 70.39 ± 1.94; difference 1.2σ — no tension between measurements. 67.4 (Planck) is not a measurement but a model fit. The question becomes: W must explain the relic angle θ* itself.

**V-392 (contradiction, open; R-135).** Sound speed from the postulates: light gives pressure ρ/3 (A-104), figure-eights give inertia, c_s² = c²/(3(1 + R)). Recombination is atomic physics under our law of rate. With the own h = 0.735, figure-eight share from deuterium and mass share from the ruler, θ is 3.3–3.5 % above the data (4.1–4.5σ) with photons and neutrinos; 16 % above with photons only. An exact hit needs h = 0.683–0.687 — the same as the official fit. So under A-105 W itself arrives at the far rate, and the local rate (z < 0.15) and the rate the relic needs differ by ~7 %. Without the neutrino share of light the relic misses by 16 %: W needs a light share ~1.7 of the photon share (candidate: a neutrino background, not derived).

**X-13 (open).** The relic angle under A-105 with the local rate h = 0.735 is 3.3–3.5 % above the data (4σ); an exact hit needs h ≈ 0.685. Updated on 27.11 after the tests below: the far rate is confirmed three ways, the local rate is 0.735 (Cepheids) or 0.710 (red giants); the causes listed earlier were tested and rejected or not supported; open are a new postulate about N, the full spectrum with other masses (CLASS) and the calibrator sample.

| Test | Result | Records |
| --- | --- | --- |
| Far rate without the relic: BAO ruler + sound until drag + deuterium | h = 0.6835 ± 0.0102 (b1), 0.6891 ± 0.0106 (b2); −3.5σ and −3.1σ from 0.7354 | V-393, R-136 |
| Local share of redshift, 1 + z_obs = (1 + ε)(1 + z) (no global "now", A-99) | free ε = −0.0003; at h = 0.685 Δχ² = +28.8 — rejected | V-393, R-136 |
| Local zone with another rate, z_v = 0.03–0.5 | Δχ² = +11.6 … +29.5 — rejected | V-393, R-136 |
| Supernova scale split by host mass, colour, stretch, survey, calibrator distance, flow z | all h = 0.717–0.752, differences ≤ 1.9σ; needed shift 0.153 mag not seen | V-394, R-137 |
| Cepheid link vs TRGB and JAGB (D-83) | Cepheid − TRGB = −0.024 ± 0.020 mag, Cepheid − JAGB = −0.041 ± 0.027; full TRGB sample h = 0.710 ± 0.011; tension 4.8σ → 2.2σ | V-395, R-138 |
| Early share growing like light from the first minutes (D-84) | at h = 0.735 needs +12 % of light, Y_p = 0.258 (3.6–3.9σ) — rejected; at h = 0.710 Y_p = 0.2524–0.2529 (2.1–2.3σ) — rejected, borderline | V-396, R-139 |
| Author's assumption: the start is a 2D circle, density falls from centre to edge | no gradient (≤ ~13 % per Hubble distance, 95 %); the only direction is a local dipole of ~320 km/s that disappears at z > 0.15 | V-397, R-140 |
| Full CMB spectrum TT/TE/EE (D-85, CAMB bridge) | links cluster: h = 0.6867, χ² = 239.2/213, p = 0.105 — passes; h = 0.7354: 14368; early share at 0.7354: 403.8; at 0.710: 277.4 (p = 0.002); smooth links: 70 000–109 000 | V-398, R-141 |
| Another medium at recombination: m_e or α changed | exact θ* at h = 0.735 needs m_e +4.83 % or α +2.37 %; at 0.710 +2.31 % / +1.14 %; at 0.685 ≈ 0; ruler passes (weak check) | V-399, R-142 |
| What may change in W | α cannot (V-171); mass m = M₀/N is global and not forbidden, but then all masses change together | V-400 |
| All masses ×λ for z > 200 | h = 0.735: +4.05 % (links unchanged) / +3.20 % (links heavier too); h = 0.710: +1.93 % / +1.54 %; ruler passes (weak check) | V-401, R-143 |
| Law of N from the postulates | N is set by the total number of points, which is constant (A-91): N constant; heavier masses not supported | V-402 |

**V-393 (computed; R-136).** The Hubble tension in W is not local and not due to redshift. It is a disagreement of two rulers at the same distances: candles give 0.735, sound gives 0.685. Either the sound length is ~7 % shorter (r_d ≈ 138 Mpc instead of 148) or the candle scale is shifted by ~0.15 mag.

**V-394 (computed; R-137).** The supernova scale is uniform. If the candle scale is to blame, it can only be the link "geometry → Cepheids", which these data do not test.

**V-395 (computed; R-138).** Three different stars (young Cepheids, old red giants, carbon stars) on one geometry give one scale to ~2 %. H₀ by red giants: 72.8 on the same supernovae as the Cepheids; 71.0 ± 1.1 on the full sample (with elliptical galaxies of the Fornax cluster). Part of the shift comes from the calibrator sample, not the method. The disagreement with the far rate 68.5 falls from 4.8σ to 2.2σ but does not vanish.

**V-396 (computed; R-139).** An early share growing like light cannot alone remove the tension: too much helium. Helium allows at most ≈ +5 % of light. Condition (author's question, 27.11): laboratory reaction rates, neutron lifetime, masses and α are taken as in our medium; if the medium at the start was different and changed them, the conclusion is wrong. Quasars (D-81) show constancy to z ≈ 4; there are no direct data at z ≈ 1100 or in the first minutes.

**V-397 (computed; R-140).** The author's assumption was not recorded as a postulate ("This is an assumption. A check will give us the real picture"). Not confirmed within z < 2.3. Not excluded: a gradient on scales beyond z ≈ 2 with us near the centre; the relic is the same in all directions to ~0.1 %.

**V-398 (computed; R-141; bridge CAMB).** Pre-material links cluster around figure-eights like matter; they cannot be smooth. With smooth links the potential wells for sound are shallow and the peak heights are wrong (χ² 300 times larger). With clustering, its own shares (ruler, deuterium) and the far rate h = 0.687, the whole spectrum including the tail passes (p = 0.10); only the initial A_s and n_s are free. Mathematically this is the same as cold dark matter; the difference is the reading: links with a constant number (V-391), not particles. The far rate 0.687 is now confirmed three ways: spot angle, sound ruler, full spectrum. An early share does not save the local rate.

**V-399 (candidate; R-142).** A heavier electron (+4.8 %) or larger α (+2.4 %) at recombination would fit the angle and the ruler at h = 0.735. Quasars require the change to vanish between z ≈ 1100 and z ≈ 4. The full spectrum with another m_e was not computed (needs CLASS).

**V-400 (candidate).** α is set by the local structure (the share of uncancelled exchanges at the free end, cancellation inside the electron size, Poisson placement); by V-171 α does not change while the whole unfolds (confirmed by D-36). So the α branch of V-399 contradicts the postulates and is excluded. Mass is a property of the global structure: m = M₀/N (V-17, V-247); the postulates do not say N is fixed during unfolding, so a different mass at the start is neither forbidden nor derived. Mass ratios are expressed through α, so at fixed α all masses change together.

**V-401 (candidate; R-143).** If all masses were 3–4 % heavier at recombination, the spot angle and the sound ruler agree with h = 0.735. In W this would mean the number of ticks N per closure was 3–4 % smaller then and reached today's value by z ≈ 4. Not derived.

**V-402 (candidate; rests on V-317).** The law of N from the postulates: N is constant. N is not a local quantity; it is set by the order of the single trajectory of the whole (V-247, A-14). By the main path V-317 (b) N is set by the total number of points, constant by A-91 — so N was the same at recombination (in agreement with lunar ranging, D-62). V-401 does not follow from the postulates; it would need a new postulate. The only thing that changes with time in W is the state of the trajectory (saturation, V-154, A-41); if N depended on saturation, repetitions would grow and masses would increase — lighter in the past, the opposite sign. From the present postulates "another medium at recombination" is not derived.

## Where it stands (27.11)

- **Derived or confirmed by data:** no expansion; redshift as a ratio of link densities with the same kinematics as a scale factor; the law g² = g₀²[Ω_r u⁻⁴ + Ω_m u⁻³ + Ω_B] (A-105) fits the BAO ruler, supernovae and the relic spectrum; the mass share is 0.040 figure-eights + 0.256 pre-material links; pre-material links cluster like matter.
- **Far rate:** h = 0.685–0.689 from three independent routes.
- **Local rate:** 0.735 with Cepheids, 0.710 with red giants.
- **Rejected as causes of the tension:** a local share of redshift, a local zone, a centre-to-edge gradient, non-uniform supernova scale, the Cepheid link, an early share (helium and spectrum tail), a varying α (postulates).
- **Not supported:** heavier masses at recombination (N constant, V-402).
- **Open:** X-13. Possible routes: a new postulate about N; the full spectrum with other masses (CLASS); the calibrator sample. Also open: the neutrino share of light in W; why the pre-material links are ~6.3 times the figure-eights; their distribution in galaxies (rotation curves, lensing) was not tested.

## Programs

All in [`alpha/`](../alpha), with rules (`*_rules_*.md`, committed before the run) and output (`*_output.txt`): `relic_scale`, `density_share`, `torus_light`, `reverse_light`, `torus_volume`, `sn_processing`, `sn_luminosity`, `sn_shape`, `sn_axes`, `light_test`, `bao_density`, `relic_light`, `mass_share`, `relic_sound`, `hubble_variants`, `candle_scale`, `cepheid_link`, `early_share`, `sky_gradient` (+ `sky_gradient_split`, an observation outside the rules), `relic_tail`, `recomb_medium`, `all_masses`. Data files: `desi_dr2_bao_data.txt`, `cchp_trgb_jagb_2024.txt`, Planck `COM_PowerSpect_CMB-*-binned_R3.0*.txt`; Pantheon+ is read from its public release.
