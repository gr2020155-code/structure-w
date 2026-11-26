# Structure W: what can be tested

26 November 2026 · Theory: Hrachya Danielyan

Record IDs (V-, R-, D-, A-) refer to `system/matrix.json`, where they are written in Cyrillic: В, Р, Д, А.

## Differences from standard physics

There are two differences in the measurable range. (1) In the structure the W boson does not screen charge: it is a pair of parts already counted (V-358), so above ~160 GeV α grows faster than in the standard calculation (V-294). (2) Charged-lepton pairs get the same production excess as quark parts: ε = α = 0.73 % for e⁺e⁻ → μ⁺μ⁻, against 0.17 % in QED (V-355, V-356). For now both are hidden in the uncertainty of the data.

| What to measure | Structure W | Standard physics | Data today |
| --- | --- | --- | --- |
| Production of high-mass lepton pairs at the LHC | higher by 1.6 % (200 GeV), 3.2 % (500 GeV), 4.4 % (1 TeV), 5.7 % (2 TeV) | baseline | CMS 13 TeV: uncertainties 12–160 %, not decisive |
| Ratio of neutral current (γ, Z) to charged current (W) at the same pair mass | grows with mass a few percent faster | baseline | not published in this form |
| e⁺e⁻ → μ⁺μ⁻ relative to the point-like rate (normalised to Bhabha or γγ) | +0.73 % | +0.17 % | CMD-3 1.0017 ± 0.0063, BaBar 1.0040 ± 0.0111: not decisive; 0.28 % needed (V-356, V-357) |

How to test: a precision of 1–2 % is needed at masses around 1 TeV. In the neutral-to-charged current ratio the uncertainty of the quark distributions in the proton almost cancels, so this is a clean test (V-299 item 9). Risk for the structure: the Standard Model with the W contribution already describes 17 precise electroweak observables (p = 0.68), and the structure must reproduce the same numbers without that contribution.

## Agreement with data

The structure reproduces measured numbers from its postulates. These are explanations, not predictions: the numbers were known in advance, so all of them have the status "candidate".

| What | Structure W | Data | Record |
| --- | --- | --- | --- |
| 1/α at large distances | 137.280 ± 0.130 | 137.036 (+1.9σ) | V-360 |
| Growth of α between 1.81 and 6.07 GeV² | (306–400)·10⁻⁵ | (440 ± 78)·10⁻⁵ (OPAL) | V-330 |
| Shape of the running of α near the electron size (the numbers 5/6 and γ) | identical | Uehling potential | V-286 |
| Slope of the running of α | 2/(3π) | 2/(3π) | V-289 |
| Charges within one generation | 0, 1/3, 2/3, 1 | neutrino, d, u, electron | V-304 |
| Number of "colours" | 1, 3, 3, 1 | leptons colourless, quarks 3 | V-304 |
| Signs of quarks in the proton and neutron | u = +2/3, d = −1/3, neutron charge 0 | uud, udd | V-305 |
| Lepton generations | three, no fourth | no fourth found | V-186 |
| Baryons of three same-sign fragments | exist (fragments on different axes) | Ω⁻, Ξcc⁺⁺ decay only weakly | V-306 |
| α in strong gravity | unchanged, if the density does not change the unevenness of places | white dwarf G191-B2B: Δα/α no more than ~10⁻⁴ | V-298 |
| Muon and tau decays | one per 3.5·10¹⁷ and 4.4·10¹² proper cycles | same | V-301 |
| Electric/gravitational force ratio of two electrons | (q·N_e)² with α/G = 4π²q² in W units | 4.17·10⁴² | V-235, V-269 |

## Not testable yet

These claims lie beyond current measurements, or coincide with standard physics where measurement is already possible.

| Claim | Number | Why it cannot be tested | Record |
| --- | --- | --- | --- |
| Step of space | a = 2π·l_P = 1.016·10⁻³⁴ m | measurements reach ~10⁻¹⁹ m | V-271 |
| Bound on the step without conditions | a ≤ 74·l_P | same | V-252 |
| Fraction of the charge source at the step itself | q₀ = 0.12289, 1/α at the step ≈ 66.2 | same | V-359 |
| New charged particles heavier than ~8 TeV | not needed for the running of α | absence cannot be proven, only not found | V-292, V-307 |
| Gravity at short distances | G constant, α/G changes with distance | gravity is measured no closer than ~5·10⁻⁵ m; agrees with standard physics | V-272 |
| Maximum particle mass | exactly the Planck mass | unreachable energies | V-246 |

Open, and would give new numbers if derived: the weak length λ_F (and with it the W mass and the lifetimes of μ, τ and the neutron) and absolute masses, which are calibration for now (A-14, V-301, V-302).

## Weak points

The α result rests on three unverified links; each could shift it by a few tenths.

- **The W boson does not screen charge** (V-294, V-299): not yet tested, either through the W mass or with LHC data.
- **"The pair remainder is the fraction q₀, not α"** (V-295): the alternation rule was written on 18 November before the calculation, but this reading was made on 25 November, when the target was known.
- **The energy → distance conversion for the hadronic part** (R-72): the γ correction is derived for logarithmic running; for hadrons it is an approximation.

Particle masses in all calculations are measured, not derived (A-14). The rules of all calculations were written before running (folder `alpha/`, files `*_rules_2026-11-25.md`).

### Sources

- [CMS, Drell-Yan at 13 TeV, JHEP 12 (2019) 059](https://arxiv.org/abs/1812.10529)
- [Gfitter, global electroweak fit](https://project-gfitter.web.cern.ch/Standard_Model/)
- [PDG 2026, W mass and width](https://pdgweb.lbl.gov/2026/reviews/rpp2026-rev-w-mass.pdf)
- [Berengut et al., α at a white dwarf, PRL 111, 010801 (2013)](https://arxiv.org/abs/1305.1337)
- [PDG 2026, Ω⁻](https://pdglive.lbl.gov/Particle.action?node=S024) · [PDG 2026, Δ(1232)](https://pdglive.lbl.gov/Particle.action?node=B033)
- [LHCb, Ξcc⁺⁺ lifetime, PRL 121, 052002 (2018)](https://arxiv.org/abs/1806.02744)
- [Alves et al., running electroweak couplings, JHEP 02 (2015) 007](https://arxiv.org/abs/1410.6810)
