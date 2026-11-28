# The hadrons chain: parts, transfers, annihilation and the nucleon

28 November 2026 · Theory: Hrachya Danielyan

This file translates the matrix records from 27–28 November on hadrons: how a whole closes, annihilation, the colour factors, neutrinos on nucleons, anti-parts and the size of the nucleon. It is a reading aid. The authoritative text is the Russian matrix [`system/matrix.json`](../system/matrix.json). If the two differ, the matrix is right. How to read IDs, statuses and the data and fit rules: see [`particles_chain.md`](particles_chain.md).

**Words for W objects.**
- *Part*: a fragment of a direction with k of the three axes (V-304). d, s, b have k = 1; u, c, t have k = 2. Parts exist only inside a whole.
- *Transfer*: a complement loop that moves an axis i → j. There are 6 oriented transfers (V-355); it is the W analogue of a gluon.
- *Complement*: the medium near a part carries the axes the part lacks (V-343 item 3).
- *Meeting region*: where the acts of two forms rebuild into new forms (V-209).

## The chain at a glance

| Step | What | Records | Status |
| --- | --- | --- | --- |
| 1 | The closure type of a whole follows from the dimension of its parts; a loop closure exists only in 1D | V-409, V-410, S-39 | derived (item 1), candidate |
| 2 | Electron + positron: a rebuild, not destruction; the photon count follows the spin | V-411, D-89 | candidate |
| 3 | p p̄ at rest: equal weight of configurations, meeting region = diameter of the outgoing form's circuit | V-412, V-413, R-147, R-148, D-90 | computed, candidate |
| 4 | p̄ n at rest: a prediction made before the data | V-414, R-149, D-91 | computed |
| 5 | Colour factors: transfer-on-transfer / transfer-on-part = 2 | V-415, R-150, D-92 | candidate |
| 6 | Neutrino on a nucleon: growth with energy; meeting region of the W-form size | V-420, R-151, D-93 | candidate |
| 7 | ν̄/ν on a nucleon of parts only is 1/3; the data need anti-parts | V-421, R-152 | computed |
| 8 | Anti-parts are the taking ends of transfers; the complement stops at the nucleon boundary; nucleon size 0.841 fm | V-422, R-153, R-154 | candidate |

## 1. How a whole closes

**V-409 (candidate).** The closure type of a whole follows from the dimension of its parts. For Ω⁻ the closure is the point O of a cross, X = 2⁴ = 16.
1. **Derived** (from A-07, A-10, A-53, A-55). The closure is the element common to the parts (A-10).
   - Parts of a one-dimensional form (electron, A-07) can only line up as a chain on a line; the common element is the middle loop. This is a loop closure (proton).
   - Parts of a two-dimensional form (a cross) cannot line up on a line. They meet in a point, like the figure-eights of the cross at O. This is a point closure.
2. **Calculation** (R-146). X = 16 gives m = 1681.1 MeV against 1672.43 (+0.52 %).
   - The number 16 was noticed before the rules; this is disclosed.
   - The tau, which has the same 2⁴ term, also deviates by +0.52 %.
3. **Open.**
   - Δ (ddd, uuu) also has spin 3/2, but a chain (1D) with X = 9 gives 946 MeV against 1232 (−23 %).
   - Ω_ccc is not described by this reading.

**V-410 (derived; the author's argument).** A loop closure exists only for a one-dimensional form; a two-dimensional cross with a loop closure is impossible.
- For two more loops to attach across the common loop in the same plane, they need points on it.
- *New points*: in one discrete step the loop would gain two points at once. That is a contradiction (A-99: one tick, one local change; A-100).
- *Old points P₁, P₂*: the loops then meet in a point, and the closure is a point, not a loop.
- The same argument excludes a layout over two planes. The form of V-193 (X = 81, 8.51 GeV) was withdrawn as S-39. Data agree: there are no ground-state baryons at 8.51 GeV or 1891 MeV.

## 2. Annihilation as a rebuild

**V-411 (candidate).**
1. **Derived.** The loops of e⁻ and e⁺ have opposite orientation, and their pushes cancel (A-83). Acts do not disappear (V-209, A-28): closures (mass) become transfer, that is, photon loops (V-65).
2. **Derived.** One photon is impossible. A photon has no closures of its own, so all its act is transfer and its momentum is not zero, while a pair at rest has zero transfer. So there are at least two photons, back to back, with equal energy (A-31): 511 keV each.
3. **Candidate** (at V-348). Spins add with sign: a pair has 0 or 1, a photon ±1. So spin 0 goes into an even number of photons and spin 1 into an odd number. Data D-89 agree (para-positronium → 2γ, 125 ps; ortho-positronium → 3γ, 142 ns). In official physics this is C-parity and the Landau–Yang theorem; in W it is a count of loops. The answer was known before.
4. **Open.** Why ortho-positronium lives about 1100 times longer.

**V-412 (computed; R-147).** In p p̄ at rest each outcome is weighted by the number of distinguishable configurations of forms (A-31). With a meeting volume of circuit radius R = ħ/(mc), no branch fits the data:
- the whole p p̄ and the proton give about 2 pions;
- the pion gives 3.5;
- the data have 5 (D-90).

A diagnostic after the count showed the data want R ≈ 2ħ/(m_π c). Accepting that without a postulate would be a fit.

**V-413 (candidate; R-148).** The meeting region is a ball whose radius is the diameter of the outgoing form's circuit, D = 2ħ/(mc) (A-80, A-82, V-283, V-275). The outgoing forms are born from acts that transferred at the meeting place, so each one's circuit passes through it.
- Result: ⟨n±⟩ = 3.16 and ⟨n₀⟩ = 1.87, against 3.0 ± 0.2 and 2.0 ± 0.2.
- Fractions of events with 0/2/4/6 charged: 2.8/41.2/51.0/4.9 % against 4.1/43.2/48.6/4.1 %. These were not used in the choice.
- Disclosed: the factor 2 was sought after the diagnostic. By the fit rule this is not a fit, but it is not an independent check either.

**V-414 (computed; R-149 before the data, D-91 after).** p̄ n at rest, the same rule with no new numbers (charge −1).

*Prediction*, written and pushed before the data were found: 1.11 π⁺, 2.11 π⁻, 1.83 π⁰ (total 5.05); events with 1/3/5/7 charged at 13.8/61.8/23.9/0.5 %.

*Check*, against Klempt–Batty–Richard 2005, Table 10:
- ⟨n±⟩ is 3.22 against 3.156, that is, +2 %. This is about 4σ by the error computed from the fractions, and 1.6σ by an unverified ±0.04.
- Fractions: 3, 5 and 7 charged are within 2σ; 1 charged is 13.8 against 16.4 ± 0.5 (−5σ).
- The split by π⁰ within one charged topology fails, for example 2π⁻π⁺2π⁰: 29 % against 7.0.

So equal weight of configurations gives the total number of pions to about 2 %, but not their composition. A selection rule is missing; officially this is G-parity. In W the similar rule for photons exists (V-411); for pions it is still to be derived.

## 3. Colour factors

**V-415 (candidate; R-150, D-92).** The ratio "transfer acting on a transfer / transfer acting on a part" in W is 2.
- **Derived, by enumeration.** Two of the six transfers act on any part (V-355, V-304).
- **Computed.** Four act on a transfer itself if it is counted by its two ends, as two places of action. The ratio 2.0 is at −1.3σ (ALEPH 2.20 ± 0.16) and −1.0σ (LEP 2.22 ± 0.22). Counting by variants without ends gives 1.0 and 1.5, both excluded.
- **Data choose:** a transfer acts by its ends.
- **Comparison.** SU(3) Yang–Mills gives 2.25. Its 0.25 above W is two diagonal gluons, which W does not have (A-28). Today's data do not distinguish W from SU(3); a precision of about 0.08 is needed.

## 4. Neutrinos on nucleons

**D-93 (data).** World average on an isoscalar nucleon (Formaggio–Zeller 2012):
- σ_ν/E = (0.677 ± 0.014)·10⁻³⁸ cm²/GeV;
- σ_ν̄/E = (0.334 ± 0.008)·10⁻³⁸ cm²/GeV;
- flat at about 30–200 GeV;
- ν̄/ν = 0.493 ± 0.016.

**V-420 (candidate; R-151).** Growth with energy is the number of outgoing configurations in a region of the W-form size.
1. **Computed.** If the weight of an outcome is the number of final configurations (A-31, V-413), it grows as s = 2E·m_part, so σ/E is flat for E far above the masses. The ratio between 200 and 30 GeV is 1.0005. Speed does not enter.
   - The normalisation "final over initial configurations" gives σ/E ∝ 1/E and is excluded.
   - A region set by the neutrino wavelength is also excluded.
2. **Computed.** The size of the meeting region is the diameter of the W-form's circuit, 2ħ/(M_W c). At the meeting the W-form, a pair of parts (V-349, V-348 item 3), is born first, and the lepton and the part come out of it.
   - With this size the W → eν and μν shares are equal: 1.000 against 1.008 ± 0.021 (D-69).
   - A region set by the outgoing lepton, the literal V-413, would give eν/μν ≈ 9·10⁶ and is excluded.
3. **Comparison.** At MeV on a proton the same rule gives σ ∝ E_e·p_e, the official shape, to 1.5 %.
4. **Disclosed.** All branch outcomes were known before the run.
5. **Open.** The absolute size needs the size of a point (A-89). The slowdown above about 10 TeV needs the exchange law (V-366).

**V-421 (computed; R-152).** ν̄/ν on a nucleon of parts only is 1/3; the data give 0.493 ± 0.016.
- The helicities of ν and ν̄ and the participation of left-handed parts only are taken from data.
- Both branches that agree with the postulates give 1/3: halves by the 2:1 rule (A-08), and the whole by A-107 with the full angle. This is 10.3σ below the data.
- The branch giving 1/2 conflicts with the postulates (a spin-1 whole has no 2:1 rule) and was excluded before the run.
- **Requirement.** With an anti-part share f the ratio is (1 + 2f)/(3 − 2f). The data need f = 0.161 ± 0.017, close to the official antiquarks (about 15–17 % of the momentum).

**V-422 (candidate; R-153, R-154; the main branch chosen by the author).** Anti-parts of a W nucleon are the taking ends of transfers.
1. **Reading.** A transfer acts by its ends (V-415): the putting end acts like a part, the taking end like an anti-part. So anti-parts exist only near parts, always paired with an extra part, and "parts minus anti-parts" = 3.
2. **Computed** (R-153). A whole nucleon always has 9 axes. The proton has 5 in parts and 4 transfers; the neutron has 4 and 5. With the full complement the anti-part share is 0.25 (measure "axis"; the ends split a transfer equally, as in A-97). That is +5.2σ from 0.161: it fails.
3. **Why the complement stops at the nucleon.** The complement exists because a part is incomplete (A-25, A-28, V-343 item 3). Outside the nucleon only the whole is seen (V-304 item 4, A-107), and it is closed: there is nothing to complement. The part size r₀ = 2ħ/(m_part c) = 1.26 fm is larger than the nucleon, so inside the nucleon the medium's complement is partial; the other parts of the whole complete the rest.
4. **Computed** (R-154). Strength g(r) = (π/12)/ln(r₀/r) (V-355, V-343).
   - "Up to the part size, full" fails in all six combinations of conversion and measure.
   - "By the neutrino's resolution" also fails in all six, and it gives a share falling with energy.
   - "At the boundary of the whole" passes in some combinations at every size tried: 0.42, 0.84 and 0.80–0.90 fm.
5. **Author's choice** (after the count, disclosed). Boundary at 2D_N = 4ħ/(m_N c) = 0.841 fm, the diameter of a ball of radius D around a point of the circuit (V-283 item 3, as in V-413). Measure: axis. Conversion: φ = g/(1 + g).
   - Result: φ = 0.392, f = 0.141 (−1.2σ), ν̄/ν = 0.472 (−1.3σ).
   - The share does not depend on energy, as in the data.
6. **Consequence.** The nucleon size in W is 0.841 fm. Electron scattering gives 0.80–0.90 fm (D-03), which agrees. The precise 0.8409 fm from muonic hydrogen comes through QED and cannot serve as a check (D-67).
7. **Open.**
   - The asymmetry d̄ > ū in the proton (SeaQuest), since the type of anti-part at a taking end is not defined.
   - In W the share should not grow with energy; the official gluon share grows slightly. This needs a check.

## Predictions from this chain

| Claim | W | Data today | Record |
| --- | --- | --- | --- |
| p̄ n at rest (made before the data) | ⟨n±⟩ = 3.22, ⟨n₀⟩ = 1.83 | ⟨n±⟩ = 3.156: total within 2 %; split by π⁰ fails | V-414 |
| p̄ with nuclei, p p̄ in flight | same rule, no new numbers | not computed yet | V-413 |
| Colour factor ratio C_A/C_F | 2.0 | 2.20 ± 0.16; SU(3) 2.25 | V-415 |
| Anti-part share of the nucleon | does not depend on energy | ν̄/ν flat at 30–200 GeV | V-422 |
| Nucleon size | 0.841 fm | 0.80–0.90 fm (scattering) | V-422 |

## Programs

All in [`alpha/`](../alpha), with rules (`*_rules_*.md`, committed before the run) and output (`*_output.txt`): `ppbar`, `ppbar2`, `pbarn` (+ `pbarn_compare`), `colour_ratio`, `nu_xsec`, `antipart_share`, `complement_boundary`. Pantheon+ and other external data are read from their public releases.

## Sources

- C. Amsler, Rev. Mod. Phys. 70 (1998) 1293, hep-ex/9708025; C. Amsler, arXiv:1908.08455 (LEAR review)
- E. Klempt, C. Batty, J.-M. Richard, Phys. Rep. 413 (2005) 197, hep-ex/0501020, Table 10
- G. Dissertori, hep-ex/9705016 (ALEPH); S. Bethke, hep-ex/0406058 (LEP summary)
- J. A. Formaggio, G. P. Zeller, Rev. Mod. Phys. 84 (2012) 1307, arXiv:1305.7513
- Positronium: arXiv:1805.06384 (review of measurements)
