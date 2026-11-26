# The α and G chain: matrix records in English

26 November 2026 · Theory: Hrachya Danielyan

This file translates the matrix records on which the main result rests: 1/α = 137.280 ± 0.130 far from the charge (data: 137.036, +1.9σ) and the link α/G = 4π²q². It is a reading aid. The authoritative text is the Russian matrix [`system/matrix.json`](../system/matrix.json); if the two differ, the matrix is right.

## How to read

- **IDs.** In the matrix, IDs use Cyrillic letters. Here they are written in Latin: А = A (author's postulate), Д = D (data), Р = R (calculation), В = V (statement), С = S (withdrawn), Х = X (contradiction). A-31 here is А-31 in the matrix.
- **Status** of a statement: *derived* (follows from the named postulates without extra assumptions), *computed* (output of a program whose rules were written before it was run), *candidate* (agrees with data, not proven), *open*. Postulates, data and calculation records have no status.
- **Rests on** is the `proof` field: the records the statement is built from. **Condition** is the `cond` field: a premise accepted conditionally, so the statement reads "if the condition holds, then …". `check_matrix.py` enforces that a *derived* or *computed* statement never rests on a *candidate*.
- Labels inside a text, such as \[derived\] or \[candidate\], give the status of that single item.
- **Data rule (author, 26.11).** Only measured data count as a check. A result of official theory may be used only if it is derived in W; otherwise it is marked "bridge, model". Earlier comparisons with 1/α(M_Z) = 128.9 (a QED calculation) are therefore not data checks.
- **Fit rule (author, 26.11).** A choice first seen in the data is not a fit once it is explained strictly within the structure: it follows from the postulates and excludes the other branches. It is still disclosed when a branch was noticed before the rules were written. "Fit" means a choice that rests on the data alone.
- **Not included.** The full dependency closure of V-360 has 230 records. This file contains the 76 records that the α result and the α–G link use directly. Records cited below but not translated (for example V-109, V-216, R-52, R-98) are in the Russian matrix.

## The chain at a glance

| Step | What | Records |
| --- | --- | --- |
| 0 | Every particle's pattern moves at exactly c (condition of most of the chain) | V-12 |
| 1 | One mechanism for gravity and charge; α = c_g·q², G = c_g·c³a²/(4π²ħ), α/G = 4π²q² | A-69, A-31, V-245, V-235, V-269 |
| 2 | c_g = 1 from the exchange size; step a = 2π·l_P | A-80, V-271 |
| 3 | G is constant, only α runs with distance | V-272, V-270 |
| 4 | Screening of charge pushes happens in the medium | A-81 … A-88 |
| 5 | Edge of screening at the electron size (5/6 and γ) | V-283, V-285, V-286 |
| 6 | Slope 2/(3π) = 4·(1/6)·(1/(2π))·2 | V-287, V-288, V-289 |
| 7 | Quark charges, colours, contribution N_c·Q² | A-15, V-304, V-305 |
| 8 | Fraction at the source q₀ from alternation; sign of charge | A-70, A-71, V-293, V-295, V-303 |
| 9 | Placement of places is Poisson; q₀ = 0.12289 ± 0.00004 over 44 clouds | A-89, V-307, R-76, R-113, V-359 |
| 10 | The W boson is not a separate species (a pair of parts already counted); contribution 0 | V-294, V-299, V-358, R-114 |
| 11 | Check by data only: with quark parts by masses the chain falls short (lifted in step 12) | D-48, R-91, V-329, V-330 |
| 12 | Hadronic part from R(s) data through the link "screening ↔ pair production" | D-68, D-76, V-337, R-97, V-338 |
| 13 | Excess of production of parts; it decreases as 1/ln(E/E₀), E₀ = mass of a part | V-339, V-341, R-99, V-343 |
| 14 | Amplitude c = π/12 from the completion loop; parts do not screen in the axis channel | V-352, V-354, R-101, R-109, V-355 |
| 15 | Excess tail above 40 GeV and the narrow resonances | R-100, V-342, R-115, V-360 |
| 16 | Assembly: 1/α = 137.280 ± 0.130 against 137.036 (+1.9σ), no fitted numbers | V-360 |

Superseded: steps 9–11 of the earlier version compared the chain with 1/α(M_Z) = 128.9 and gave 136.86 ± 0.14 (V-307, R-72, R-77). Under the data rule that number is a QED calculation, so those comparisons are not data checks. The current result is assembled from `end_to_end.py` (step, slope, edge, charges), `alpha_data_test.py` (R-91), `q0_more_clouds.py` (R-113), `hadron_dispersion.py` (R-97) and `tail_gee.py` (R-115).

Open in this chain (X-01, spin, was resolved on 26.11 by the author's choice A-97, on which V-348 and so V-358 rest): how a W⁺W⁻ pair is produced in W (V-358 item 3); the reading "the pair remainder is q₀" (V-295).

## Author's postulates

| ID | Postulate |
| --- | --- |
| A-01 | One indivisible change — one step — one tick. Length is a number of steps, time is a number of ticks; c is one step per tick. |
| A-02 | A link is a stable sequence of acts of change (moves of points). The placement is dynamic; a state is a place in the sequence; there are no fixed points. |
| A-03 | Local self-closure is an exchange of places between neighbouring points. An exchange is an indivisible act (step, tick). The side of an exchange is a property of the act. A circle is the path of an exchange, one step long. |
| A-06 | A point is indivisible; its size is one step; there are no gaps between points. (Refined by A-89.) |
| A-15 | Charge is the direction of a loop; a quark fragment has one third of a direction. |
| A-17 | The side of an exchange is set by the global figure-eight (level 3); reversing its direction is the transition "world ↔ antiworld", not annihilation. Refinement (14.11): a one-dimensional particle has no plane of its own; the global figure-eight sets the reference of the sign — which sign is "ours" — not a direction in space. |
| A-27 | The placement of places has no lattice: isotropy "gives it, it cannot be otherwise". |
| A-31 | O12: the indistinguishable cannot have different weight. |
| A-39 | E1: the amount of change per step is divided between transfer and internal change by squares; hence E² = (mc²)² + (pc)², E = mc² at rest, time dilation, a common speed limit, inertia. Derived from the structure (register E1). |
| A-42 | There is no choice: a local system that keeps the global system self-consistent has exactly one consistent trajectory. |
| A-43 | Every point exchanges at every step; inside a black hole there is no such step — a phase transition to the next level. |
| A-64 | Each layer of points around a charge passes the lag on to the next (as in the derivation of 1/r for black holes). |
| A-69 | The analogy between gravity and charge is not accidental: the mechanism is one and the same. |
| A-70 | The side of a free end's outward exchange must alternate with the previous one — like the layers in gravity. |
| A-71 | The constant component of the side is by analogy (not literally) with the centre of mass of two bodies in gravity. |
| A-76 | Places are not arbitrary: they are restricted by the same consistency condition as motion; the uniqueness of the trajectory (A-42) refers to the whole — the placement of places and the motion together. |
| A-80 | The size of an exchange is the distance between the partners, i.e. the diameter of the circle of the act (A-03). |
| A-81 | Charge pushes (unscreened outward exchanges) are screened only in the medium; inside an act there is nothing to screen with. |
| A-82 | The screening boundary is the size of the electron, because an electron changing its orientation in space can touch points within the diameter of its size. |
| A-83 | Pushes of the two loops of one electron-type medium form cancel each other: the loops O↔A and O↔B of opposite orientation (A-07) give pushes of opposite side, and a pair of such pushes is screened (V-274). |
| A-84 | The electron's loops alternate by orbits: one orbit — one loop (O↔A), the next orbit — the other (O↔B). |
| A-85 | A medium point holds the "+" side of a push until a "−" push arrives. |
| A-86 | A charge push is screened at a medium point whose side was set by a loop of an electron-type medium form. |
| A-87 | Electron-type medium forms are induced by the charge itself; far from the charge (beyond the electron size) there are none. |
| A-88 | The side is held by the exchange act as a whole: both partner points of one act hold its side. |
| A-89 | Refinement of A-06: the size of a point is one step on average; individual places have different sizes, there are no gaps between points. |
| A-97 | (author, 26.11) The inertia of a particle is split equally between the two ends of the lag exchange: at the core (moves along the orbit) and at the global figure-eight (does not move along the orbit); half of the inertia moves along the orbit. Choice of branch a2 (V-347, V-348); resolves X-01. |

## Data

**D-48.** Change of α with distance, direct measurements. (1) α at large distances: 1/α = 137.036 (D-45). (2) OPAL (2006, small-angle Bhabha scattering): Δα(−6.07 GeV²) − Δα(−1.81 GeV²) = (440 ± 58 ± 43 ± 30)·10⁻⁵; a constant α is excluded at more than 5σ; leptons alone do not explain the growth at 3σ. (3) L3 (2005, Bhabha scattering at 1800–21600 GeV²): the growth of α equals the expected one with a factor C = 1.05 ± 0.07 ± 0.14; a constant α is excluded. In total: α grows at short distances; the size of the growth is confirmed directly to about 15 %. *Data rule (26.11): 1/α(M_Z) = 127.95 (MS-bar scheme) and 128.9 (with a dispersion calculation), the slope 2/(3π)·Q² per species and "the growth starts at the Compton length" are QED calculations, i.e. model; checks of candidates against these precise numbers do not count as data checks. The previous text is kept in `was`.*

**D-68.** e⁺e⁻ → hadrons: R(s) = σ(hadrons)/σ(μ⁺μ⁻), 1400 points, √s 0.3–188.7 GeV; below 2 GeV the sum of exclusive channels (CMD-3, SND, BaBar, KLOE), above that inclusive measurements. Weighted means: R = 3.51 ± 0.17 at 5.0–9.3 GeV and 3.96 ± 0.14 at 12–40 GeV. Narrow resonances (not resolved in the table): Γ_ee — J/ψ 5.53 keV, ψ(2S) 2.33, Υ(1S) 1.340, Υ(2S) 0.612, Υ(3S) 0.443 keV. *Level: processed data with a model part (corrections for initial-state radiation; joining of experiments).*

**D-76.** Γ(e⁺e⁻) and masses of the narrow resonances (PDG 2026, pdgLive): J/ψ 5.53 ± 0.10 keV, 3096.900 MeV; ψ(2S) 2.33 ± 0.04 keV, 3686.097 MeV; Υ(1S) 1.340 ± 0.018 keV, 9460.40 MeV; Υ(2S) 0.612 ± 0.011 keV, 10023.4 MeV; Υ(3S) 0.443 ± 0.008 keV, 10355.1 MeV. They agree with the values entered from memory in R-97.

## Statements

### V-12 · candidate
*Rests on:* A-01, A-02, A-11

The orbital speed of a particle (of its pattern) is exactly c. Postulate P1 sets the speed of change, not of the pattern: this is assumption W3, selected by data (μ_e ≥ μ_B).

*The condition V-12 is listed in each statement that depends on it.*

### V-245 · derived
*Rests on:* A-69, A-64, A-31, A-39, A-65, V-109, V-112, V-113, V-116, V-117, V-128, V-234 · *Condition:* V-12

One mechanism (A-69): the source gives continuations into the medium, each layer passes them to the next whole, equally among its cells, and acts are indistinguishable (A-31) — so the rule is one. Consequences. (1) The same form of the law: 1/r² for the field, 1/r for the potential, superposition of fields, rearrangement at speed c (V-113, V-117, V-234). (2) Differences come only from the source: for gravity the source is closures (repetitions inside the form); they have no outward side, so the field has no sign and is not screened (V-116). For charge the source is unscreened outward exchanges; they have a side (the loop orientation), so the field has a sign and cancels in neutral matter (A-65). (3) The fraction of a particle's own energy that the medium carries beyond its orbit is, for both, the square of the source rate per tick, halved (the energy of a cell is the square of its share, A-39): electric α/2 = 3.65·10⁻³, gravitational Gm_e²/(2ħc) = 8.8·10⁻⁴⁶; the ratio (q·N)² does not depend on the constants of the rule. (4) Hierarchy of forces: gravity is weaker by (q·N)² = 4.2·10⁴² (two electrons; 1.24·10³⁶ for two protons, as in the data) because the electron closes once per N ≈ 10²² ticks, while an unscreened outward exchange happens once per ~12 ticks.

### V-235 · candidate
*Rests on:* A-31, V-234, V-112, V-128, V-17 · *Condition:* V-12

Acts are indistinguishable (A-31), so the rule of transmission through layers (κ in F4) is the same for gravity and charge. The source of gravity of a particle at rest is 1/N closures per tick (V-128); the source of a charge field is q net outward exchanges of the charge loop per tick. Then for two electrons electric / gravitational = (q·N_e)², without κ and N₃; from the data (α/(Gm_e²/ħc) = 4.17·10⁴²) q·N_e = 2.04·10²¹. Consequences: (1) the electron's figure-eight has exactly one free end at every tick, so q ≤ 1 and N_e ≥ 2.04·10²¹, step a ≤ 1.2·10⁻³³ m — about 6·10¹² times tighter than the LEP bound (7·10⁻²¹ m); (2) if the unit of coupling is one full outward exchange per tick, then α = q², q = 1/11.706, N_e = m_P/m_e = 2.39·10²², step a = 2π·l_P = 1.016·10⁻³⁴ m — the Planck length comes out without being put in. The check "q = number of free ends per tick" fails: the electron and the proton (with parallel ends) have 1 each, the muon cross has 3, but the muon charge is 1. Not derived: why the charge loop gives out a net exchange once per 11.7 ticks.

### V-269 · candidate
*Rests on:* V-239, V-245, V-17 · *Condition:* V-12

The link between G and α as a formula. From Gm²/ħc = c_g/N² and m = h/(cNa), the product m·N = h/(ca) is the same for all particles, so **G = c_g·c³a²/(4π²ħ)**: G is set by the step a and the coupling c_g alone, masses drop out. Together with α = c_g·q²: **α/G = 4π²q²·ħ/(c³a²)**, and in W units (a = c = ħ = 1) **α/G = 4π²q²**. The common factor is c_g; G cannot be expressed through α alone until one more number is derived: q (or, equivalently, the step a in units where α is given). At c_g = 1 the formula gives a = 2π·l_P — the same as V-235.

### V-271 · candidate
*Rests on:* A-03, A-31, A-43, A-62, A-80, V-17, V-109, V-112, V-128, V-246 · *Condition:* V-12

c_g = 1 (given A-80). The size of an exchange is the distance between the partners, the diameter of the circle of the act a/π (A-80); the horizon of a particle with N = 1 (A-43) equals it: 2Gm/c² = a/π with m = h/(ca), hence c_g = 1. Independently consistent: the atom (A-62: no medium, β = q·q, by data β = α — R-59; at c_g = 1/2 the atom size would be 68.5 and the binding energy ×4), the lag is carried by both participants (A-31), the energy of a cell is the square of its share without the half (A-39), the mass ladder 1 + 1/2 per figure-eight in units of 1/q² (V-201). Consequences with data: α = q², q = 1/11.706, step a = 2π·l_P = 1.016·10⁻³⁴ m, N_e = m_P/m_e = 2.389·10²². Weak links: β = q·q (V-228, V-250 are candidates); the horizon formula 2Gm/c² is taken from general relativity. The previous text (factor 2 not closed) is kept in `was`.

### V-270 · candidate
*Rests on:* V-244, V-269, A-31, V-208 · *Condition:* V-12

If G and α are of one kind (common coupling c_g, A-31), then G too has a variable part — the same closure factor as α in V-244: α(L) = c_g·q²·Π(L), G(L) = c_g·c³a²/(4π²ħ)·Π_g(L). Two variants. (a) The factor belongs to c_g, i.e. it is common: Π_g = Π. Then α(L)/G(L) = 4π²q² at any distance — the variable parts cancel, the G–α link is exact at any L; G grows at short distances just as α does (at 10⁻¹⁸ m by about 7 %). (b) The factor belongs to q, i.e. to charge only: by V-208 it comes from the alternating orientation of medium layers around the end of a form, and mass has no end with a sign — then G does not change, and the ratio α/G changes with L. The data do not distinguish: α changes only at distances below the electron's Compton length (4·10⁻¹³ m), gravity is measured no closer than ~5·10⁻⁵ m; neither variant contradicts the data. The choice between (a) and (b) belongs to a postulate: whether closure sits on the common coupling or on the side of charge.

*Note in the matrix (25.11): by V-272 variant (a) is excluded — only exchanges with a side are screened, gravity is not screened; (b) remains.*

### V-272 · candidate
*Rests on:* V-208, V-116, V-245, V-216, A-03, A-63 · *Condition:* V-12

Why q changes as α changes with distance (1/137 → 1/128). (1) \[derived\] The electron does not depend on the distance from which it is viewed; only what lies between the charges changes — the medium (as in V-208). So q(L) is the fraction of unscreened exchanges that reach distance L, not a property of a single form: part of the screening is done by the medium. (2) \[derived\] An opposing pair can screen only what has a side; closures (mass) have no side, the gravitational field is not screened (V-116, V-245 item 2). So only charge changes with distance: variant (b) of V-270; variant (a) is excluded — the ratio α/G changes with L. (3) \[derived\] 1/11.706 is the value after all screening in the medium, at distances beyond the electron's Compton length (further out α does not change, D-48). The own q₀ of the exchange itself is larger. (4) \[derived\] The far α is the same for all charges (V-216), although particle sizes differ (the muon is 207 times smaller): screening is set by the medium, not by the size of the source; its edge is the size of the electron, the lightest charged form. (5) \[open\] The natural candidate is q₀ = 1 at the step (each act has its own side, A-63; inside one act there is nothing to screen with). The standard extrapolation of the logarithmic law (not data) gives 1/α ≈ 105 at the Planck distance, not 1; the electron alone between the step and its own length gives 10.9 to 1/α. So either q₀ < 1, or screening below 10⁻¹⁹ m is not logarithmic. Problem 1 changes: derive q₀ and the law of screening in the medium, not only the shape of the electron.

*q₀ < 1 was later derived in V-293.*

### V-283 · derived
*Rests on:* A-82, A-87, V-234, V-245 · *Condition:* V-12

Where the boundary at the electron size and the steep part at the step come from. (1) \[derived\] Neither the lag (∝ 1/r, V-234) nor the pushes (∝ 1/r², V-245) stop at the electron size: both extend to infinity. So, given A-87, the boundary is set not by them but by the form itself: the charge induces forms that reach it, i.e. the charge within their diameter (A-82). (2) \[derived\] In the reading "the orbit of the form passes through the point of the charge", the points of the form lie at distances from 0 to D from the charge, and their share per point of a layer near the charge is ∝ 1/r² (numerical check: f·r² is constant from 10⁻⁴·D to 0.2·D). So the steep part of screening at the step comes from the geometry of the form alone, without a separate rule for pushes. Its size (≈ 61 units) is not derived: it depends on the number of such forms. (3) \[derived\] The logarithmic part (f ∝ 1/r) is given by forms of the same size D passing through the charge, if the form changes orientation in all directions (A-82): its points lie on the surface of a sphere of diameter D through the charge, and the share of such points per point of a layer is ∝ 1/r up to D itself (numerical check: f·r is constant from 0.01·D to 0.75·D). So the steep part (1/r²) is the form in one orientation (a circle), the logarithm (1/r) is the same form over all orientations (a sphere); forms of other sizes are not needed for the logarithm. (4) \[derived\] Unfolding the form from the step to D (sizes equally per scale) does not give a logarithm: for circles of growing size f ∝ ln(D/r)/r², not 1/r (numerical check). (5) \[computed\] A sphere through the charge cuts off sharply at D: this is the earlier sharp edge at the diameter, which disagrees with the data (V-276 item 2: the decrease from leptons is 5.28 against 3.94 in real distance; the data need a soft fall-off with an effective edge D/8.2). (6) \[open\] What makes the edge soft and smaller than the diameter.

### V-285 · candidate
*Rests on:* R-64, A-07, A-12, A-39, A-72, A-74, A-80, V-228, V-283 · *Condition:* V-12

The screening edge is a figure-eight divided between two loops \[candidate, chosen by data\]. (1) \[data, mathematics\] The data profile decomposes exactly into pairs of shares x and 1−x: weight x(1−x), each pair an exponential with decay length R·√(x(1−x)) (R-64). This gives both puzzling numbers of the edge: 5/6 is the mean of ln√(x(1−x)) with weight x(1−x), γ comes from the fall-off being exponential (a constant fraction of loss per step). (2) \[candidate\] Reading in W: the medium form induced by the charge is a figure-eight (A-07, A-12) of size D (A-80); its diameter is divided between the loops in shares x and 1−x; the weight is the product of the shares, the same rule as β = q·q in the atom (V-228); the decay length is the geometric mean of the loop radii (the square of the fraction of loss per step = the product of the loops' turning fractions, "by squares", A-39). (3) \[computed\] Strict symmetry of each form (A-72, A-74: x = 1/2) gives the data tail exactly (length R/2 = ħ/(2mc)), but the edge with an error of +0.09 in 1/α for leptons; the data exclude this. So induced forms are symmetric only in the second sense of A-74, "no ground for a difference": the distribution of shares is symmetric (x ↔ 1−x), each individual form is not. (4) \[candidate, found after the number\] Slope 2/(3π) = (4/π)·∫x(1−x)dx: if the weight is the product of shares without normalisation, the 1/6 in the slope is the same integral; the factor 4/π is not derived. (5) \[open\] Derive from the postulates: (a) why the induced form is divided with weight x(1−x); (b) why the decay length is the geometric mean of the loop radii; (c) the factor 4/π. If (a)–(c) are derived, the edge and slope of α are obtained without fitting.

### V-286 · candidate
*Rests on:* R-65, A-07, A-10, A-12, A-20, A-31, A-39, A-80, A-84, V-285 · *Condition:* V-12

The division weight and decay length from the postulates \[candidate, chosen by data among 4 readings\]. (1) \[candidate\] Weight: an induced form is a choice of partner in each loop (exchanges O↔A and O↔B, A-07); there are n_A·n_B configurations with division n_A : n_B; indistinguishable configurations have equal weight (A-31), so the weight is ∝ x(1−x). (2) \[candidate\] Decay length: the form goes round the loops in turn (A-84) at one step per tick, so the fraction of time spent in a loop is the fraction of its length; per step it turns by a/r; the amount of change per step is divided by squares (A-39), so the measure of loss of direction is the mean square turn, λ = 1/√⟨(1/r)²⟩ = R·√(x(1−x)). (3) \[computed\] Of the four combinations (pairs or weight 1) × (squares or linear mean), only "pairs × squares" agrees with the data: edge −2.1037 and the half at 0.150·D; the longest decay length D/4 = ħ/(2mc) is common to all (R-65). The linear mean gives R/2 for any division and misses by 0.09 in 1/α. (4) Weak points: reading A-31 as "weight = number of partner pairs"; reading A-39 as "loss of direction is the root-mean-square turn"; the exponential fall-off follows from A-20 as a fraction of configurations, not derived. (5) \[open\] Slope 2/(3π): the factor 4/π and the number of induced forms.

### V-287 · candidate
*Rests on:* A-31, A-86, A-87, V-276 · *Condition:* V-12

The slope 2/(3π) as the loss of a push per scale. (1) \[derived, mathematics\] By the data 1/α grows linearly in ln r (D-48). This is the same as saying that the fraction of a push lost per e-fold distance is α/(3π) = q²/(3π): d ln q/d ln r = (κ/2)·q², κ = 2/(3π). So in 1/α the contributions of scales add up like segments of a path (P = 2/α is the path in orbits, V-233). (2) \[candidate\] The factor q² = q·q is the same product rule as β = q·q in the atom (V-228) and the division weight x(1−x) (V-286): the charge push is coupled to the form with its share q, the form to the push with its own. Then one number remains from the structure, G = 1/(3π) per kind of charged form: how many points with the side of a form fall on a scale and what fraction of the push is screened on them. *(Here G is a coefficient of the slope, not the gravitational constant.)* (3) \[open\] G requires a rule that the postulates do not have: how many forms one unscreened charge push induces and how long a form lives. (4) \[computed, discarded\] The combination "one form per tick, lifetime = time to cross the diameter, the transverse share ⟨sin²⟩ = 2/3 is screened" gives exactly 1/(3π) but is inconsistent: a form induced at the point of the charge that has gone D steps is an arc from the charge, its points give 1/r² (the steep part), not 1/r (the logarithm, the sphere of V-283). Found after the number, geometries mixed — not a derivation. (5) \[candidate\] The factor N_c·Q² for quarks: Q² comes from the product rule (A-15: a quark fragment has a third of a direction); N_c = 3 could be the number of axes of space — not checked. *(Later derived in V-304.)*

### V-288 · candidate
*Rests on:* R-67, A-31, A-80, V-286, V-287 · *Condition:* V-12

Slope 2/(3π) = (2/π)·⟨2x(1−x)⟩. (1) \[derived, given the rule of V-286\] The division weight of the figure-eight is the number of partner pairs n_A·n_B (A-31). As a fraction of all pairs of points of the diameter (partners are equal) this is 2x(1−x); averaged over equally possible positions of O it is (n+1)/(3(n−1)) → 1/3 (R-67). The edge (V-286) sets only the shape of the weight, not its normalisation, so the 1/3 in the slope is not fitted to the slope but follows from the same rule as the edge: one rule gives the 5/6 and γ of the edge and the 1/3 of the slope. (2) \[candidate, found after the number\] The factor 2/π remains. Reading: the distance between partners is the diameter (A-80), the path of an exchange is a half-turn, the ratio of the diameter to the half-turn path is 2/π (as in V-276 item 3). It is not a check until it gives some other number. (3) The rival reading 1/3 = ⟨cos²⟩ over directions in 3D (V-276 item 3) is not distinguished by the data; the argument for (1) is that it is the same rule already chosen by the edge. (4) \[open\] Derive 2/π; explain N_c = 3 for quarks.

### V-289 · candidate
*Rests on:* R-68, A-07, A-31, A-82, A-83, A-86, V-283, V-288, A-88 · *Condition:* V-12

The factor 2/π in the slope \[candidate; after A-88 the choice by data is removed\]. (1) \[derived, given V-283\] π is geometric: the held points of the induced form lie on a sphere of diameter D through the charge; a layer at distance r from the charge crosses this sphere along a circle of length 2πr, and the loss per e-fold distance is α·H/(2π·D²) — the 1/(2π) is not fitted (checked by calculation, R-68). (2) \[derived, given V-288\] 1/3 is the mean ⟨n_A·n_B⟩/D² = 1/6 over equally possible positions of O, multiplied by the count below. (3) \[candidate\] The integer count factor c·p = 4: c = 2 because partners are equal (A-31), the same 2 as in the share 2x(1−x) of V-288; p = 2 because an exchange act holds its side with both partner points (A-88). Of the four combinations (c, p) ∈ {1, 2}² only this one agrees with the data (R-68), and it is also the one that follows from the postulates. Then the slope **2/(3π) = 4·(1/6)·(1/(2π))·2** — π from geometry, 1/6 from the division rule, 4 from the count, 2 from α = q² (1/α = q⁻²). (4) The earlier weak point p = 2 (conflict with A-86) is removed by A-88 (V-290). One caveat remains: the set of readings (sphere geometry, division rule, count c·p) was assembled when the number 2/(3π) was already known. (5) Result for α: π, 1/3 and the form of the law come from the structure; the integer factor 4 from A-31 and A-88; screening at the step (≈ 61) and N_c = 3 are open. *(N_c = 3 was later derived in V-304.)*

### V-304 · candidate
*Rests on:* R-75, A-15, A-39, A-31, A-27, A-79, V-268, V-286, V-289

Quarks as fragments of a direction along three axes \[candidate\]. (1) \[derived\] A direction in three-dimensional space is divided by squares into three axes (A-39), none is singled out (A-31), so the share of an axis is 1/3. A fragment is a part of the axes, k of 3: share (charge) |Q| = k/3, number of different fragments C(3, k). (2) \[computed\] One rule gives three sets of numbers of one generation at once (R-75): charges 0, 1/3, 2/3, 1 (neutrino, d, u, electron); numbers of variants 1, 3, 3, 1 — quarks have three ("colour"), leptons one; contributions to the running of α 0, 1/3, 4/3, 1, i.e. N_c·Q² (D-48). "Colour" is the choice of the fragment's axes. (3) \[derived\] Not "Lego": the variants are equal states of the whole, not parts; each screens with its own share squared, total C(3, k)·(k/3)². Adding "parts" (square of the sum) gives 1 and 4 instead of 1/3 and 4/3 — excluded by the data. (4) \[candidate\] Only whole directions exist freely — all axes (electron) or none (neutrino); fragments (1 or 2 axes) exist only inside a whole. This is a reading of quark confinement. (5) \[open\] Why the u fragment has the sign opposite to the electron and d (+2/3 with −1 and −1/3); why there are three generations (perhaps three figure-eights at O, A-79). *(Both answered in V-305.)* The agreements were found when the data were known; the argument is that three independent sets of numbers come out of one rule.

### V-305 · candidate
*Rests on:* V-304, A-23, A-15, A-52, A-07, V-216, V-186, V-300, A-79

The sign of the u quark and three generations. (1) \[derived, given V-304\] The sign of a fragment is the direction of its loop (A-15). The proton has three loops: (+), (−), (+) (A-23). If the "+" loops carry a axes each and the "−" loop b axes, the proton charge (2a − b)/3 must equal the electron charge, 1 (V-216). There are two solutions: a = 2, b = 1 and a = b = 3. The second is three whole charges, not fragments: wholes exist freely (V-304 item 4), and the proton would decay into them. The only remaining one: fragments of two axes on the "+" loops (+2/3), one axis on the "−" loop (−1/3). This is u, d, u: the "+" sign of u is not a separate puzzle but the sign of the proton's loop. (2) \[derived\] The neutron is a proton plus an electron-type figure-eight (A-52): +2/3 − 1/3 + 2/3 − 1 = 0; in the language of fragments one u (+2/3) becomes d (−1/3), the difference −1 is the electron charge. (3) \[candidate\] W⁻ decays are the division of all three axes into complementary parts (k, 3 − k): (3, 0) — an electron and an antineutrino, (1, 2) — d and anti-u; this is what is observed (consistent with V-300). (4) \[derived, V-186\] There are three generations: point O in three-dimensional space holds at most three figure-eights (e — one, μ — a cross, τ — three), there is no fourth charged lepton. The number 3 for "colour" (choice of the fragment's axes, V-304) and for generations (number of figure-eights at O) is the same: the dimension of space, in different roles. (5) \[open\] Why the quarks of the other generations (c, s; t, b) repeat the leptons: they are fragments of the cross and the tau form — not checked.

### V-293 · candidate
*Rests on:* A-65, A-70, A-71, V-216, R-70, V-291 · *Condition:* V-12

Variant 2: q₀ < 1 at the source — the postulates themselves require it. (1) \[derived\] By A-70 the side of a free end's outward exchange must alternate with the previous one; so at the source successive outward pushes come in pairs (+, −, +, …), and only the constant component of the side remains unscreened (A-71). Its share is q₀, and it is less than 1. The earlier record "q₀ = 1 per tick" (V-273 item 3) did not take A-70 into account — fixed at the source. (2) \[derived\] q₀ is the same for all charged forms: the charges of e, μ, p are equal to 10⁻²¹ (V-216). (3) \[computed\] With the derived edge and slope and the known fermions, the data require 1/α at the step ≈ 65–67, i.e. q₀ ≈ 1/8.1–1/8.2; if the W boson screens with the standard sign — ≈ 108–110, q₀ ≈ 1/10.4 (R-70). The target is not firm (light quarks, W boson). The closeness to 1/8 was found after the number and is not a check. (4) \[open\] Derive the constant component of the side at the source from A-70, A-71 (the old problem V-260, now with a target at the step rather than far away) and decide how the W boson behaves in W.

### V-295 · candidate
*Rests on:* R-52, R-71, A-70, A-39, A-63, A-03, V-269, V-293, V-286, V-289, V-294 · *Condition:* V-12

q₀ at the source from the alternation of the side (A-70) and the assembly of 1/α \[candidate\]. (1) \[computed, rule of 18.11\] A free end changes its outward exchange side to the most opposite available one (A-70); the pair remainder by squares (A-39) is 0.1226–0.1234 on Poisson clouds and 0.1252 on a lattice with shifts (R-52, R-71). (2) \[candidate\] By V-269 α = q², where q is the fraction of unscreened outward exchanges; the pair remainder by squares is a fraction, so it is q₀, not α (on 18.11 it was wrongly compared with α, V-248 item 1). Then 1/α at the step = 1/q₀² ≈ 65.7–66.5. (3) \[computed\] Together with the derived edge and slope (V-286, V-289), without the W-boson contribution (V-294) and with the measured masses of charged fermions: 1/α far away ≈ 135.7–138.2 against 137.036 by the data. The number falls inside a band whose width is set by the uncertainty of the hadronic part (light quarks) and by the dependence on the placement of places (2 %). (4) Weak points: reading (2) was made when the target was known; the direction of the remainders is random (V-248 item 2), so the sign of the charge must come from the direction of the loop (A-15), not from the remainder; the placement of places in W has no lattice (A-27), but which one exactly is not set (A-76), hence the 2 %; masses are calibration (A-14); the edge, slope and the role of the W boson are candidates. (5) Testable: with this arrangement α grows faster than standard above ~160 GeV (V-294), and below 5·10⁻²⁰ m without new particles. *(The sign question was resolved in V-303, the placement in V-307.)*

### V-303 · candidate
*Rests on:* R-74, A-15, A-17, A-07, A-39, A-70, V-295 · *Condition:* V-12

The sign of charge in the q₀ rule. (1) \[derived\] Charge is the direction of a loop (A-15); a one-dimensional particle has no plane of its own, and the global figure-eight sets only the reference of the sign, not a direction in space (A-17, A-07). So the sign of charge is not tied to the direction of the pair remainder in space: the unscreened fraction of a pair is the modulus of the remainder by squares (A-39), |(n₁ + n₂)/2|², always with the sign of the loop (reading Z1). (2) \[computed\] Then q₀ = 0.1232 — the same as in R-71, and the charge does not wander: the sign is constant, only the spatial direction of the remainder is random (|⟨r⟩| = 0.002 — this is what goes to zero; V-248 item 2 referred to it, not to the charge). (3) \[computed\] The reading "the sign is carried by the projection of the remainder onto the direction to the observer" (Z3) gives q₀ = 0.041 and 1/α at the Z distance = 661 against 129.72 — excluded; the reading "projection onto the loop axis" (Z2) is undefined for the electron (A-17). (4) This agrees with the data: the field of a charge is spherically symmetric, a charge has no direction. The weak point of V-295 item 4 about the sign is removed; what remains is the reading "the remainder is q₀" (made when the target was known) and the choice of placement by data. *(The latter was removed in V-307.)*

### V-307 · candidate
*Rests on:* R-76, A-02, A-43, A-31, A-27, A-89, A-76, A-42, V-295, V-296, V-303 · *Condition:* V-12

The placement of places from the postulates and the final 1/α \[candidate\]. (1) \[derived\] By A-89 the size of a place is not fixed, so the postulates contain no rule "points no closer than d". By A-02 and A-43 places are constantly mixed by exchanges, and an exchange does not change the density. By A-31 indistinguishable placements have equal weight, and equal weight of all placements with a given density means independent uniform points, i.e. a Poisson cloud; isotropy (A-27) holds. A-76 and A-42 require uniqueness of the whole but give no local rule ordering the places. So the Poisson placement is singled out by the postulates, and choosing it by data (V-296, V-297) is no longer needed. (2) \[computed\] With it q₀ = 0.12282 ± 0.00013 (8 clouds, spread 0.29 %), and 1/α at the Z-boson distance = 129.54 ± 0.14 against 129.72 ± 0.05 by the data — a difference of −0.18, 1.3σ; far away **1/α = 136.86 ± 0.14 against 137.036 (−0.13 %)**. (3) What remains outside the calculation: the energy → distance conversion for the hadronic part (the γ correction), the role of the W boson (V-294, V-299), the reading "the remainder is q₀" (V-295). Each of them can shift the result by a few tenths.

*Note (26.11, data rule): the comparison here with 129.72 and 137.036 rests on the model value 1/α(M_Z) = 128.9 and does not count as a data check. The data check is V-329 → V-338 → V-360. q₀ over 44 clouds: V-359.*

### V-294 · candidate
*Rests on:* A-08, V-27, A-83, A-86, A-88, A-64, A-28, V-293 · *Condition:* V-12

The W boson in screening. (1) \[candidate\] Spin 1 is a form whose orientation returns after one cycle, i.e. one loop, not a figure-eight (a figure-eight has ħ/2 by the 2:1 rule, A-08, V-27). (2) \[derived, given (1)\] A charge push is screened by a pair of loops of opposite sides of one form (A-83, A-86, A-88); a form of one loop has no such pair — it has nothing to screen with. (3) \[derived\] There is no anti-screening (amplification of charge by the medium) in W: a layer passes on what it received (A-64); if the side matches, it passes it unchanged (V-274 item 4); a closed system as a whole does not change (A-28) — the medium cannot add pushes of the charge's side. (4) \[candidate\] So the W boson does not enter the running of α (contribution 0, not the −42.8 of the standard calculation), and the target at the step is 1/α ≈ 65–67, q₀ ≈ 1/8.1–1/8.2 (V-293). (5) \[testable consequence\] Above ~160 GeV α in W grows only from fermions; in the standard calculation the W boson slows the growth by about 1.1 units of 1/α per e-fold distance, by 4 TeV about 4.4 units. There is no direct measurement of α above 160 GeV with this precision: the sensitivity estimate (arXiv:1410.6810) is a forecast for the 14 TeV LHC. (6) \[counterexample to item 3\] The strong coupling α_s decreases at short distances (standard data; not recorded in the matrix as data) — this is anti-screening. If the strong interaction in W works by the same screening mechanism, item 3 is wrong; the strong interaction is not described in W, so item 3 applies to charge only. (7) \[open\] Describe the W and Z bosons as forms (V-177).

*Note (26.11): item 1 (one loop) is wrong for W by V-348; the conclusion of item 4 (contribution 0) now rests on another argument — V-358.*

### V-299 · open
*Rests on:* V-294, D-52, D-53, V-177, D-54, D-55, D-56

Checking V-294 ("the W boson does not enter the running of α") against data. (1) \[data\] The W boson exists (D-52); the question is not its existence but whether it screens charge in the medium. (2) \[data\] There are no direct measurements of α above 160 GeV with the required precision (V-294 item 5). (3) \[data\] The precise electroweak measurements of the LEP era do not unambiguously require boson loops: the need for them depends on the choice of input quantity (G_F or the W width), and with the W width fermion loops alone suffice (D-53). So the data found do not refute V-294 — but do not confirm it either. (4) Caveats: the data are from 1996, modern global fits were not checked; "boson loops" is broader than screening of charge by the W boson; the weak interaction and the W and Z bosons are not described in W (V-177), so W does not yet give an exact prediction for electroweak quantities. (5) \[data, 25.11\] Modern: the W mass is measured to 7.7 MeV (D-54); the full Standard Model with boson loops describes 17 precise electroweak observables with p = 0.68 and predicts the W mass within 0.7σ (D-55). The 1996 argument (use the W width instead of G_F) is now weaker: the W width is known to 2.4 %, G_F to ~10⁻⁶. I did not find an analysis "fermion loops only" on modern data. (6) \[derived\] So V-294 is not directly refuted by data, but it is under threat: a structure in which the W boson does not screen charge must reproduce the same precise electroweak numbers without this contribution, and W has no description of the weak interaction (V-177). Until it is built, the result 1/α ≈ 137 rests on an unverified link. (7) \[open\] Describe the W and Z bosons as forms and compute at least one precise quantity (the W mass or the mixing angle). (8) \[data, 25.11\] A direct test with high-mass lepton pairs: without the W-boson contribution α above 160 GeV is larger, and more pairs are produced than in the standard calculation, by about 1.6 % at 200 GeV, 3.2 % at 500 GeV, 4.4 % at 1 TeV, 5.7 % at 2 TeV (estimate from the W coefficient of R-60). The CMS data found (D-56) have uncertainties of 12–160 % in these bins — an order of magnitude larger than the effect, so they do not decide. A decision needs 1–2 % precision at masses around 1 TeV and parton distribution functions obtained without these same data. (9) \[search, 25.11\] I did not find more precise measurements: there is no direct measurement of the running of α above 160 GeV at the LHC (the proposal arXiv:1410.6810 is a calculation and forecast); the measurement of high-mass lepton pairs on the full 2015–2018 data set (ATLAS, 140 fb⁻¹) was presented at EPS-HEP 2025, but I did not find a paper with uncertainties; searches for non-resonant deviations on 139 fb⁻¹ see none, but normalise the background to data and depend on parton distribution functions at 5–10 %. Conclusion: the 3–6 % effect is currently hidden in the uncertainty of the quark distributions in the proton. A test in which this uncertainty almost cancels: the ratio of pair production through the photon and Z (neutral current) to production through the W (charged current) at the same mass — in W the former grows by a few percent, the latter does not. I did not find this ratio in published form.

### V-329 · computed
*Rests on:* R-91 · *Condition:* V-307, V-286, V-304, V-294, V-313, V-271

The α chain with quark parts by masses disagrees with the data \[computed; the contradiction X-10 was lifted on 26.11 — V-338\]. (1) \[computed, R-91\] Without relying on QED the chain gives 1/α far away 136.10 ± 0.14 (edge S3), 136.34 (S1), 135.81 (S2); by atom recoil it is 137.036 (D-45). The discrepancy is 0.5–0.9 %, 5–9σ, for every edge branch. (2) The earlier agreement (136.86, −0.13 %, R-77) rested on the model value 1/α(M_Z) = 128.9; without it the discrepancy is five times larger. (3) \[open\] Where the error lies is not established. Suspects: q₀ (placement, V-307), the edge (V-285, V-286), the masses of quark parts (Rb, V-313), the form size D_f = 2ħ/(mc), the step a = 2π·l_P (V-271). +0.94 in 1/α is missing. Fitting is not allowed: a rule from the postulates is needed.

*Note (26.11): X-10 is lifted. With the hadronic part taken from R(s) data (V-337, V-338), 1/α far away is 136.77 ± 0.18. The calculation with parts by masses stays correct; what was wrong was replacing hadrons by quark parts.*

### V-330 · computed
*Rests on:* R-91 · *Condition:* V-286, V-304, V-294, V-313

The growth of α at short distances agrees with the direct OPAL measurement \[computed, under conditions\]. (1) \[computed, R-91\] The difference Δα(−6.07) − Δα(−1.81 GeV²) in W is 306–400·10⁻⁵ over six branches (edge S1/S2/S3 × momentum reading M1/M2); OPAL gives 440 ± 78·10⁻⁵; all within 1.7σ. Without the W boson and with quark parts from the Rb reading. (2) The data do not choose a branch: the OPAL precision is 18 %.

### V-337 · candidate
*Rests on:* A-87, A-39, A-28, A-31, A-64, V-75 · *Condition:* V-304, V-289

Screening and pair production are the same forms \[candidate\]. (1) \[candidate\] An induced form (A-87) exists within its size; in a collision with a reserve of change E ≥ 2Mc² (A-39) the same form becomes free, and it appears as a pair because the total charge does not change (A-28). (2) \[candidate\] Indistinguishable configurations have equal weight (A-31), so the weight of a species in production (R = σ(forms)/σ(μ⁺μ⁻)) is the same as in screening — N_c·Q² (V-304): the slope of 1/α per e-fold distance at the scale s is κ·R(s). (3) \[derived, given (1), (2)\] Influence is causal (V-75) and transfer is linear (A-64); a linear causal response is recovered from its absorptive part (the Kramers–Kronig relation, mathematics): 1/α(0) − 1/α(Q) = (1/(3π))·∫R(s)·Q²/(s(s + Q²)) ds. (4) \[data, R-97\] Check of the weights: R by the data is above Σ N_c·Q² by 5.4 % (5–9 GeV, 1.0σ) and 8.0 % (12–40 GeV, 2.1σ). In standard physics this is the strong-interaction correction; in W the strong interaction is not described, so the excess is not explained. (5) Weak points: the production weight is not derived directly (it is taken from A-31); above 40 GeV, where there are no data free of the Z boson, the weight is taken from W (11/3) — if the 5–8 % excess persists up to the step, the result for 1/α changes by units.

*Note (26.11): the dispersion relation of official science is used here as a consequence of W, as the data rule allows. The excess of production over the weights: V-339; its explanation: V-343, V-355.*

### V-338 · computed
*Rests on:* R-97, R-91 · *Condition:* V-337, V-286, V-304, V-294, V-307, V-333

1/α far away from the postulates and the data on hadron production: 136.77 ± 0.18 against 137.036 \[computed, under conditions\]. (1) \[computed, R-97\] The step, q₀, edge, slope, leptons and t are from W (R-91); the contribution of quark parts below 40 GeV is replaced by R(s) data (D-68) through the link V-337; narrow resonances are added by Γ_ee. (2) The deviation −0.26 (−0.19 %), 1.5σ, is within the 2σ criterion written before the calculation; the expectation before the calculation (136.8–136.9) is disclosed in the rules. (3) The contradiction X-10 is lifted given the candidate V-337: the shortfall was in the quark parts, not in the step, q₀ or the screening length (V-333, V-334). (4) The error is set by the systematics of the R data (±0.11) and by q₀ (±0.14). What remains: the R excess of 5–8 % (V-337 item 4) and the weight above 40 GeV.

*Note (26.11): with the decreasing excess above 40 GeV — V-342, V-355, V-360.*

### V-339 · computed
*Rests on:* R-98 · *Condition:* V-337, V-304

Hadrons are produced 8–9 % more often than the part weights N_c·Q² give \[computed; by data alone\]. (1) \[computed, R-98\] ε = R/ΣN_c·Q² − 1: +0.116 ± 0.014 (3 GeV), +0.053 ± 0.029 (8 GeV), +0.078 ± 0.021 (14 GeV), +0.080 ± 0.009 (31 GeV); on average +0.087 ± 0.007 — the excess is real (8σ). A decrease with energy is not proven by these data (−1.9σ). (2) So the production weight in V-337 item 2 is incomplete: something of about 9 % is added to N_c·Q² that W does not have. In standard physics this is the strong interaction (and near 31 GeV the Z-boson contribution); in W neither is described (V-294 item 6, V-177). (3) For α: the result of V-338 (136.77) is right only if the excess disappears above 40 GeV; if it persists up to the step, 1/α far away ≈ 136.77 + 2.6 — too much. (4) \[open\] Needed: a W rule for the excess (what adds ~9 % to the production of parts) and its behaviour at short distances.

*Note (26.11): answered by V-341, V-343, V-355.*

### V-341 · computed
*Rests on:* R-99 · *Condition:* V-339, V-337

The excess of production of parts decreases with energy — the data choose the direction "attraction of parts" \[computed; the form is a fit\]. (1) \[computed, R-99\] The Z boson alone does not explain the excess (χ² 219/7): direction (c) on its own is excluded; the Z contribution is visible only near 58 GeV. (2) A constant excess (a: completion without energy dependence) is allowed (probability 0.24), but one decreasing as 1/ln(E/E₀) (b) is better: Δχ² = 5.4 and 4.9 at a threshold of 4. (3) The amplitude c = 0.25 (E₀ = mass of a part, 312.8 MeV) or 0.33 (E₀ = m_π) is fitted; the postulates have no rule for it. (4) \[open\] Derive from the postulates why the parts attract each other and why this fraction decreases with energy as 1/ln — that will be the strong interaction in W.

*Note (26.11): the 1/ln form and E₀ = mass of a part follow from the postulates — V-343; the amplitude — V-355.*

### V-342 · candidate
*Rests on:* R-100, R-97 · *Condition:* V-341, V-338, V-337

1/α far away with a decreasing excess: 137.33 ± 0.19 against 137.036 (+1.6σ) \[candidate, chosen by data; the amplitude of the excess is fitted\]. (1) \[computed, R-100\] If the excess decreases as c/ln(E/E₀) also above 40 GeV, its contribution up to the step is finite: with E₀ = mass of a part (312.8 MeV) it is +0.56, result 137.33 ± 0.19 (+1.6σ, passes); with E₀ = m_π it is +0.70, 137.47 (+2.3σ, does not pass). (2) The data choose E₀ = mass of a part — the same scale where the induced form of a part ends (D = 2ħ/(mc), Rb). The margin of the choice is small (+1.6σ against +2.3σ). (3) If the excess does not decrease (a), it adds +2.6, too much (V-339). So the α result now rests on the form of the strong interaction in W (V-341 item 4).

*Note (26.11): with V-354 and V-355 the amplitude is no longer fitted: c = π/12, 1/α ≈ 137.35. The tail branches are checked in R-115 and V-360.*

### V-343 · candidate
*Rests on:* V-245, A-69, A-25, A-28, A-87, V-283 · *Condition:* V-287, V-304, V-341, V-294

Attraction of parts and its 1/ln decrease — from postulates already recorded \[candidate; item 1 derived under conditions\]. (1) \[derived, given V-287 item 1 and V-245\] Why 1/ln. The mechanism is one (A-69, V-245): a layer passes on what it received, and the rule is the same at all scales. Then the contributions of the scales add up like segments of a path, one per e-fold distance (V-287 item 1, mathematics), and the inverse of the coupling is linear in ln r: 1/g(r) = b·ln(r₀/r), i.e. g = 1/(b·ln(E/E₀)). This is the same arithmetic as for α (1/α is linear in ln r); with one mechanism there is no other form. So form (b) in R-99 is not a fitted form but a consequence; only the amplitude is fitted. (2) \[candidate\] Why E₀ is the mass of a part. g becomes infinite where 1/g = 0, at the distance beyond which parts never exist separately. By V-304 item 4 parts exist only inside a whole, and the induced form of a part ends at its size D = 2ħ/(mc) (A-87, V-283). So r₀ is the size of the form of a part and E₀ is its mass (Rb: 312.8 MeV). The data chose exactly this E₀, not m_π (R-99, R-100). (3) \[candidate\] Why they attract. A part is an incomplete direction (k axes out of 3, V-304); only the closed is stable (A-25), and the whole does not change (A-28). So the medium around a part completes it to a whole: induced forms near the part carry the missing axes. They do not screen but add, the effect grows with distance (anti-screening), and at r₀ it becomes complete: the part cannot be separated. This does not contradict V-294 item 3: there a layer cannot add pushes of the same side of charge; here what is added is not the side of the part but the missing axes (V-294 item 6 already restricted item 3 to charge). (4) \[open\] The amplitude — how many completion forms come per e-fold distance (b) and how g enters the production excess (ε ∝ g?) — is not recorded in the postulates. By the data c = 0.251 ± 0.017, i.e. 1/c = 3.98 ± 0.27; the closeness to 4 was noticed after the fit and is not an argument.

*Note (26.11): item 4 — V-352, V-354, V-355.*

### V-352 · candidate
*Rests on:* R-108, R-101, A-07, A-83, A-31 · *Condition:* V-286, V-288, V-289, V-343, V-294, V-345, V-304

The completion form is one loop, and its weight is three times that of a figure-eight; this does not yet give the amplitude c \[candidate; item 1 derived under conditions, items 2–3 computed\]. (1) \[derived, given V-286, V-288, V-343 item 3\] A completion form has no pair of the opposite side (V-294 item 2): it is one loop, and it adds rather than screens. The side in it is held by all pairs of points (n²/2), while in a figure-eight only by the pairs separated by the point O (n²/6): the weight of the loop is exactly 3 times larger (R-108). (2) \[computed\] The data on the growth of c with the number of part species (R-101) limit screening by parts: P/A ≤ 0.09, i.e. parts screen no more than 9 % of what the loops add. With one measure of the fraction for the loop and the part (8 loop variants, parts equally) P/A = 1/8, slightly above the limit: c grows with n more noticeably than in the data. (3) \[computed\] None of 108 branches passes. Two rules are missing: (a) the fraction of a part in the axis channel must be smaller than that of the loop (s_part²/s_loop² ≤ 0.73 with 8 variants); (b) how the coupling g enters the production excess (k): with the loop ×3, 8 variants and measure 1, k ≈ 1.3 is needed, and there is no natural value for it. (4) Noticed: in standard theory the gluon weight in these units is 11/4 of the charge form, here it is 3; P/A there is 0.06. The coincidence was noticed during the calculation and is not an argument.

*Note (26.11): item 3 — V-353, V-354, V-355.*

### V-354 · derived
*Rests on:* A-83, A-86, A-07, A-15, R-101, R-108 · *Condition:* V-304, V-343, V-274

Parts do not screen the attraction in the axis channel (P = 0) \[derived, given V-304, V-343\]. (1) Screening in W is a pair of pushes of opposite side from the two loops of one form (A-83, A-86): in an induced figure-eight the loops O↔A and O↔B have opposite orientation (A-07). (2) The orientation of a loop is the charge (A-15), so in the charge channel the induced forms of parts do screen (weights N_c·Q²). (3) The axes, however, are the choice of axes of the form itself (V-304): both loops of one induced form of a part carry the same axes. There is no opposite "axis side" inside the form — nothing to screen with; a matching side is passed on unchanged (V-274 item 4). So P = 0 exactly. (4) Check by data: c does not grow with the number of part species — P/A = 0.004, limit 0.09 (R-101, R-108). (5) This corrects the premise of the rules of R-101 and of branches p1–p3 in R-102 and R-108 ("parts screen by the same mechanism as charge"): for charge yes, for axes no.

### V-355 · candidate
*Rests on:* R-109, A-07, A-28, A-15, V-354 · *Condition:* V-352, V-287, V-343, V-342, V-304

Link between the excess and the attraction: ε = g, c = π/12 = 0.262 \[candidate; the number of loop variants and P from the postulates, the reading of k chosen by data among three\]. (1) \[derived\] The completion loop is oriented (A-07), so the transfers of axis i → j and j → i are different variants: 3·2 = 6. "Diagonal" combinations (an axis into itself) transfer nothing and change nothing (A-28), so they are not completions. So V = 6, not 8 or 3. (2) \[derived\] The measure in the axis channel is a whole axis: the loop transfers an axis entirely, fraction 1 (for charge the measure is the fraction of direction k/3, A-15). (3) P = 0 (V-354). (4) \[computed, R-109\] With V = 6, measure 1 and P = 0, one of three readings of k passes: K1, the excess equals the coupling (product rule: a pair of parts, each with fraction √g, together g — as α = q·q, V-287 item 2). c = 1/(κ·3·6) = π/12 = 0.2618 against 0.255 ± 0.018 (+0.4σ) and 0.257 ± 0.037 (+0.1σ). K2 (g/2) gives 0.131 (−6.9σ); K3 (increment per e-fold distance) gives 1/3 and 1/2 and does not pass. (5) Then there are no fitted numbers in the α chain: 1/α ≈ 137.35 against 137.036 (+1.7σ). (6) Disclosed: the branch K1·V6 was noticed mentally before the rules of R-109 were written. By the author's refinement (fit rule) this is not a fit, since V = 6, the measure and P = 0 follow from the postulates; the choice between K1, K2 and K3 is made by the data, by the fork rule. (7) \[testable consequence\] The same mechanism for charge (A-69) gives an excess in the production of charged-lepton pairs ε = α·(pair fraction) ≈ 0.7 %, while official theory gives 0.17 %. Not tested: e⁺e⁻ → μ⁺μ⁻ data with a precision better than 0.5 % are needed.

*Note (26.11): item 7 — V-356, V-357: the best data (CMD-3, 0.63 %) do not decide; 0.28 % is needed. With q₀ over 44 clouds and the tail check the result is 137.280 ± 0.130 (V-359, V-360).*

### V-358 · candidate
*Rests on:* R-114, V-304, V-337, A-31 · *Condition:* V-348, V-349, V-12

The W boson does not enter the running of α because it is not a separate species but a pair of parts already counted \[candidate; item 1 derived given V-348, V-349; item 2 computed\]. (1) \[derived\] Screening weight is counted by species of parts: a species is one division of the three axes k, weight N_c·Q² (V-304, V-337). The W consists of two parts (V-348 item 3) and is a pair of divisions already counted, (k, 3 − k): e⁻ν̄, dū, … (V-349). It is not a new division of the axes, so not a separate species; adding it would count the same parts twice. Hadrons work the same way: π, ρ, J/ψ are pairs of parts and are not added on top of the quarks. The same holds for the spin-0 partner of the W (V-351). This replaces the earlier argument of V-294 item 1 ("one loop, nothing to screen with"). (2) \[computed, R-114\] Of four branches only this one passes: contribution 0, 1/α = 137.354 (+1.7σ). A separate species with w = 1 gives 145.06 (+43σ), weight 1/4 gives 139.28 (+12σ), the anti-screening of the standard theory gives 96.90 (−217σ). Both the postulates and the data settle the fork. (3) \[open\] W⁺W⁻ pairs are produced often at LEP2 (about 17 pb at 200 GeV; figure from memory, not checked). How a W pair is produced in W is not derived; if it went through the called form (V-337), R above 161 GeV would grow by about 8 and 1/α by tens of units. So in W the W pair must be produced by exchange of parts, not through the called form; this has to be derived.

### V-359 · computed
*Rests on:* R-113 · *Condition:* V-307, V-355, V-12, A-26, V-08

q₀ over 44 clouds and the 1/α result \[computed\]. (1) q₀ = 0.12289 ± 0.00004 over 44 Poisson clouds (R-113); the earlier 0.12282 ± 0.00013 over 8 clouds (R-76) is confirmed (+0.4σ). (2) The spread between clouds equals the calculation error inside one cloud: different placements give the same q₀, and the spread is only the randomness of the calculation. The property of the placement does not affect q₀ at the 0.2 % level. (3) The result of the chain V-355 with this q₀: 1/α far away = 137.281 ± 0.130 against 137.036, +0.245 (+1.9σ). The earlier 137.354 moved closer by 0.07, and the error from q₀ fell from ±0.14 to ±0.05. (4) The error of the result is now set by the rest of the chain (±0.12: the hadronic part from R(s) data and the excess tail), not by q₀. The link "clouds for q₀" is closed.

### V-360 · computed
*Rests on:* R-115, D-76 · *Condition:* V-337, V-338, V-342, V-343, V-354, V-355, V-359, V-12

The excess tail and Γ_ee are not weak links \[computed, under conditions\]. (1) \[data, D-76\] Γ_ee and the masses of the five narrow resonances, entered in R-97 from memory, agree with PDG 2026. The narrow part, 0.150, does not change. (2) \[computed, R-115\] It is unclear to which α the width refers: α(0) or α at the resonance mass. This changes the narrow part by −0.008. It is a systematic, not a choice. (3) \[computed, R-115\] The excess tail above 40 GeV with c = π/12 and one E₀ = 312.8 MeV is +0.581. Inaccuracies of the same calculation (where t starts, the upper end) give +0.01 each. The result: 1/α far away = 137.280 ± 0.130 (±0.012 systematic) against 137.036, +1.9σ. (4) The branches excluded by the postulates do not pass by the data either: an own E₀ per species gives +3.4σ, a frozen excess +12.9σ, no excess above 58 GeV at all −2.5σ. So the α result requires the excess to decrease up to the step, as follows from V-343 item 1. (5) Neither link shifts the result by more than 0.02. The ±0.13 error is set by the R(s) data below 40 GeV (±0.11) and by q₀ (±0.05).


## Calculations

### R-72
Narrowing the band (rules written before the calculation). (a) Hadrons from data: at r_Z = ħ/(M_Z·c) = 2.164·10⁻¹⁸ m, 1/α by the data in real distance = 128.9 + κγ·20/3 = 129.72 ± 0.05; the W running from the step to r_Z without q₀ = κ·(20/3)·ln(r_Z/a) + κ·(4/3)·(ln(D_t/a) − 2.1037) = 63.25 (no light-quark masses). 1/α_W(r_Z) = 1/q₀² + 63.25: Poisson clouds (q₀ = 0.12262–0.12336) — 129.76…128.96, far away 137.08…136.28; a lattice with shifts (0.12517) — 127.07. An exact match needs q₀ = 0.12266. (b) Matérn II clouds (exclusion of close approach, A-06): q₀ = 0.12390; 0.12477 (share 1/2, two seeds) and 0.12455 (share 1/4); 1/α_W(r_Z) = 128.39; 127.48; 127.72 against 129.72.

*Note (26.11, data rule): 129.72 is obtained from the model value 128.9; not a data check.*

### R-76
`placement_poisson.py` (rules written before the calculation). Rule R-52 (A-70) on 8 new Poisson clouds (seeds 101–108, 120 000 points, 3000 ends × 60 ticks): q₀ = 0.12245; 0.12310; 0.12266; 0.12302; 0.12227; 0.12271; 0.12323; 0.12316. Mean 0.12282, spread between clouds 0.29 %, error of the mean 0.00013. 1/α_W(r_Z) = 1/q₀² + 63.25 = 129.54 ± 0.14 against 129.72 ± 0.05 by the data (−0.18, −1.3σ); far away 136.86 ± 0.14 against 137.036. An exact match needs q₀ = 0.12266.

### R-77
`end_to_end.py` — an end-to-end check by one program. Derived from the postulates: step a = 2π·l_P (c_g = 1); slope κ = 2·2·(1/6)·(1/(2π))·2 = 0.212207 = 2/(3π); edge ln(r_eff/D) = −2.1037; Σ N_c·Q² at the Z distance = 20/3, for t 4/3 (three axes); q₀ = 0.12282 ± 0.00013 (the same 8 clouds as R-76, reproduced to the last digit). From data: G, ħ, c, M_Z, m_t, 1/α(M_Z) = 128.9, 1/α far away. Result: 1/α at the Z distance 129.54 ± 0.14 against 129.72 ± 0.05; far away 136.86 ± 0.14 against 137.036 (−0.13 %). An inaccuracy was found and fixed: at n = 2000 steps ⟨n_A·n_B⟩/D² = 0.16675 instead of the limit 1/6 (κ larger by 0.05 %); a real form has D/a ~ 10²², so the limit is taken.

*Note (26.11, data rule): the comparison with 129.72 and 137.036 rests on the model value 1/α(M_Z) = 128.9 and is not a data check; without it the result is 136.10 (R-91, V-329).*

### R-91

`alpha_data_test.py` (rules written and pushed to GitHub before the calculation, commit 455b68f). The α chain without relying on QED: 1/α(0) = 1/q₀² + κ·Σ N_c·Q²·(ln(D_f/a) + E) over 9 species, D_f = 2ħ/(m_f c), masses of quark parts from the Rb reading (u = d = 312.8; s = 490.2; c = 1660.9; b = 4994.1 MeV), t directly, without the W boson. Target T1 (atom recoil, 137.036): S1 136.34, S2 135.81, S3 136.10 ± 0.14 — deviations −0.51 %, −0.89 %, −0.69 % (5.0σ, 8.7σ, 6.7σ); no branch passed. Target T2 (OPAL, 440 ± 78·10⁻⁵): S1 340/400, S2 306/364, S3 324/383 (point-like M1 / via the amplitude M2) — all six passed (−0.5…−1.7σ), the data do not choose. Contributions at S3: e 10.25; μ 9.11; τ 8.52; u 11.85; d 2.96; s 2.93; c 11.37; b 2.77; t 10.06. For information (not a fit): to agree, +0.94 is missing — that is 1/q₀² = 67.23 (q₀ = 0.12196) or an edge higher by 0.55.

### R-97

`hadron_dispersion.py` (rules written and pushed to GitHub before the calculation, commit 56f898f). 1/α(0) = 1/α_W,parts(0) + (1/(3π))·∫[R_data − R_parts] ds/s up to (40 GeV)² + narrow resonances. ∫R_data ds/s = 27.829 (stat. 0.316, syst. 1.029); ∫R_parts ds/s = 22.866; correction +0.527 ± 0.114; narrow +0.150 ± 0.003 (J/ψ 0.101, ψ(2S) 0.036, Υ 0.014). Result: 1/α far away = 136.097 + 0.527 + 0.150 = 136.774 ± 0.181 against 137.036: −0.26, −1.5σ — passes (criterion 2σ). For information: with s_max = (30 GeV)² it is 136.76, with (20 GeV)² 136.74. Check of the link: R = 3.51 ± 0.17 (5.0–9.3 GeV) against 10/3 — +5.4 %, 1.0σ; R = 3.96 ± 0.14 (12–40 GeV) against 11/3 — +8.0 %, 2.1σ; the 3σ threshold is not exceeded.

*Note (26.11): Γ_ee and the masses of the narrow resonances were checked against PDG 2026 (D-76) and agree; the systematic from the choice of α in the width is −0.008 (R-115).*

### R-99

`r_excess_shapes.py` (rules before the calculation, commit 4384db2). 9 regions (D-68), ε: 2.3 GeV +0.138 ± 0.047; 2.8 +0.145 ± 0.045; 3.4 +0.105 ± 0.008; 5.8 +0.033 ± 0.048; 8.2 +0.056 ± 0.031; 13.7 +0.078 ± 0.021; 25.8 +0.067 ± 0.013; 34.7 +0.084 ± 0.013; 57.8 +0.345 ± 0.030. Fits (all except (c) include a Z part z₁·V₁ + z₂·V₂, V₁ = s/(M_Z² − s)): (a) ε = c — χ² 7.95/6 (probability 0.24); (b1) c/ln(E/312.8 MeV), c = 0.251 ± 0.017 — χ² 2.51/6 (0.87); (b2) c/ln(E/139.57 MeV), c = 0.334 ± 0.023 — χ² 3.06/6 (0.80); (c) Z only — χ² 219/7, excluded. Δχ²(a − b1) = 5.4, (a − b2) = 4.9 at a threshold of 4: the data prefer (b), the excess decreases with energy.

### R-100

`alpha_excess_tail.py` (rules before the calculation, commit 75876e5). Excess (b) above 40 GeV up to the step energy (1.94·10¹⁸ GeV), weight 11/3, from the t threshold 5, without Z: (b1) tail +0.557 ± 0.038 — 1/α far away 137.331 ± 0.185 (+0.30, +1.6σ), passes; (b2) +0.698 ± 0.047 — 137.472 ± 0.187 (+0.44, +2.3σ), does not pass. A mental estimate made before the rules (137.32 and 137.46) is disclosed in the rules.

### R-101

`c_by_flavours.py` (rules before the calculation, commit 501998e). Amplitude of the attraction of parts c = (ε − Z)·ln(E/312.8 MeV) by groups: n = 3 (u, d, s) — 0.255 ± 0.018; n = 4 (+ c) — 0.151 ± 0.082; n = 5 (+ b) — 0.257 ± 0.037. Slope dc/dn = −0.003 ± 0.020 (−0.1σ): the data do not decide; the precision is about 8 % per species, as large as the expected growth. The expectation before the calculation ("they do not decide") is disclosed in the rules.

### R-109

`coupling_k.py` (rules before the calculation, commit 6ccf300). 36 branches (K1 ε = g and K2 ε = g/2, with V = 3/6/8, measure 1 or 1/3, P = 0 or P/A = 1/V; K3 — the increment of the pair's push per e-fold distance, c = f/2, where f is the fraction of loop variants acting on the pair). Three passed: K3·f8a (4 transfers out of 8 variants), P = 0 — c = 0.2500 (−0.3σ/−0.2σ); K1·V6·S1, P = 0 — c = π/12 = 0.2618 (+0.4σ/+0.1σ); K2·V3·S1, P = 0 — the same π/12. All that passed require P = 0 (parts do not screen in the axis channel). Two of the three were found mentally before the rules were written — disclosed in the rules; the third gives the same number. For α (tail ∝ c, R-100): c = 1/4 gives 1/α ≈ 137.33, c = π/12 gives ≈ 137.35 (data 137.036; +1.6σ and +1.7σ).

### R-113

`q0_more_clouds.py` (rules before the calculation, commit 2fd6b23). The same rule as R-76 on 32 new Poisson clouds (seeds 109–140): from 0.12216 to 0.12348, mean 0.12288 ± 0.00005, difference from R-76 +0.4σ. The spread between clouds, 0.00028, equals the error inside one cloud, 0.00027 (ratio 1.01). Main result over 44 clouds (101–108, 11, 21, 31, 41, 109–140): q₀ = 0.12289 ± 0.00004 (spread 0.24 %). 1/α far away = 137.281 ± 0.130 (the part from q₀ ±0.048) against 137.036: +0.245, +1.9σ. The expectation was disclosed in the rules and was met.

### R-114
`w_screening.py` (rules written before the calculation). Contribution of a form of mass M_W from its size: κ·w·(ln(D_W/a) + E) = 7.706·w (ln(D_W/a) = 38.417, E = −2.1037). B0 (w = 0, parts already counted): 137.354 (+1.7σ), passes; B1 (separate species, w = 1): 145.06 (+43σ); B2 (weight 1/4, bridge): 139.28 (+12σ); B3 (anti-screening −21/4, bridge): 96.90 (−217σ).

### R-115

`tail_gee.py` (rules before the calculation, commit ead765c). (1) Narrow resonances with Γ_ee from PDG 2026 (D-76): N1 = 0.1500, as in R-97. If Γ_ee is converted to α(0) (1/α(M) by the link V-337: 133.6 at J/ψ, 132.1 at Υ), N2 = 0.1423, a difference of −0.008. (2) Excess tail, base without the tail 136.700: H1 (one E₀ = 312.8 MeV) +0.581 — 137.280 (+1.9σ), passes; H3 (t from m_t) 137.288; H4 (up to m_P) 137.292 (+2.0σ); H2 (own E₀ per species) 137.484 (+3.4σ); H5 (frozen excess) 138.71 (+12.9σ); H6 (no excess above 58 GeV) 136.715 (−2.5σ). Systematic of H1 ±0.012. The expectation was disclosed in the rules and was met.
