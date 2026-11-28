# Structure W

An attempt to derive the fine-structure constant α, its link to the gravitational constant G, and the pattern of particle charges from a small set of postulates about space. Only these postulates and experimental data are used; ready-made postulates of standard physics are not.

**Author of the postulates:** Hrachya Danielyan
**Calculations and checks:** done with Claude (Anthropic) under the author's rules.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24544483.svg)](https://doi.org/10.5281/zenodo.24544483)

**Cite as:** Danielyan, H. *Structure W: the fine-structure constant α and its link to G from postulates about space.* Zenodo. https://doi.org/10.5281/zenodo.24544483 (all versions; v1.0: https://doi.org/10.5281/zenodo.24544484).

*Русская версия: [README_ru.md](README_ru.md).*

## Main result

| Quantity | Structure W | Measured |
|---|---|---|
| 1/α far from the charge | **137.280 ± 0.130** | 137.036 (+0.18 %, +1.9σ) |
| Growth of α between 1.81 and 6.07 GeV² (OPAL) | (306–400)·10⁻⁵ | (440 ± 78)·10⁻⁵ |
| Slope of the running of α | 2/(3π) | 2/(3π) is a QED calculation; the measured growth agrees to ~15 % |
| Charges within one generation | 0, 1/3, 2/3, 1 | ν, d, u, e |
| Number of "colours" | 1, 3, 3, 1 | leptons 1, quarks 3 |

There are no fitted parameters. One choice, which of three readings of the coupling enters the excess of hadron production, was made by the data; the rest of that step follows from the postulates, so by the author's rule it is not a fit (V-355). Below 40 GeV the hadronic part is taken from measured e⁺e⁻ → hadrons data (V-337, V-338). Only measured data are used as a check; the earlier result 136.86 relied on the QED value 1/α(M_Z) = 128.9 and is superseded. The chain is translated record by record in [`docs/alpha_chain.md`](docs/alpha_chain.md).

## α and G: one coupling

In Structure W gravity and charge are carried by the same exchange mechanism. Exchange acts are indistinguishable (postulate A-31), so both forces share one coupling constant c_g. The only difference is the source: for mass it is the closures of the particle's loop, for charge it is the unscreened fraction q of the charge loop's exchanges.

With a the step of space and N the number of steps in one cycle of a particle (m = h/(c·N·a)):

- α = c_g · q²
- G·m²/(ħc) = c_g / N²  →  **G = c_g · c³a² / (4π²ħ)**: masses drop out, G is set by the step a alone
- hence **α / G = 4π² q² · ħ / (c³ a²)**, and in W units (a = c = ħ = 1) simply **α / G = 4π² q²**

The exchange size equals the horizon of a particle with N = 1 (postulate A-80), which gives c_g = 1. Then:

- α = q², with q = 1/11.706;
- the step is a = 2π·l_P = 1.016·10⁻³⁴ m: **the Planck length comes out without being put in**;
- the electric-to-gravitational force ratio of two electrons is (q·N_e)² = 4.17·10⁴², with N_e = m_P/m_e = 2.389·10²².

This is how G enters the α calculation: G fixes the step a, the starting point from which α runs up to the Z distance.

**Limits.** G cannot be expressed through α alone: one more number is needed, namely the absolute mass scale, which is calibration so far (A-14). The horizon formula 2Gm/c² is taken from general relativity. Status: candidate (matrix records V-269, V-271, V-235, V-270, V-272).

**Does G run with distance like α?** No (V-272). Screening in the medium cancels only exchange pushes that carry a side; closures (mass) have no side, so gravity is not screened. G stays constant, while α grows below the electron's Compton length (4·10⁻¹³ m). So α/G = 4π²q² holds at large distances, and at short distances the ratio grows together with α. This agrees with standard physics: gravity is measured no closer than ~5·10⁻⁵ m.

## Honest status

- **All agreements were found with the answer already known.** So they are **candidates**, not proofs.
- **One prediction was made before the data (V-414, 28.11).** The pion count for p̄ n annihilation at rest was written and pushed before the data were found. The total number of charged pions agrees to 2 % (3.22 against 3.156), but the fraction of 1-prong events is off by 5σ and the split by π⁰ fails. So this is a partial confirmation, not a full one: the count works, a selection rule is missing (see [`docs/hadrons_chain.md`](docs/hadrons_chain.md)).
- The α calculation uses **measured** particle masses (leptons, m_t) and the measured rate of e⁺e⁻ → hadrons. In the structure, masses and the weak force are still calibration (postulate A-14), not derivation.
- The result rests on unverified links (see [`docs/alpha_chain.md`](docs/alpha_chain.md)):
  1. the W boson is not a separate species but a pair of parts already counted (V-358); how a W⁺W⁻ pair is produced in W is not derived;
  2. the W as two parts (V-348) rests on the author's choice A-97 (half of the inertia sits with the global figure-eight), which resolved X-01 (spin);
  3. "the pair remainder is q₀" (V-295), a reading made when the target was known;
  4. the excess of hadron production: the reading ε = g was chosen by the data among three (V-355).
- Two testable differences from standard physics: above ~160 GeV α grows faster (V-294, V-358); e⁺e⁻ → μ⁺μ⁻ pairs are produced 0.73 % more often than the point-like rate, against 0.17 % in QED (V-355, V-356). Current data are not precise enough for either: the second needs 0.28 %, the best today (CMD-3) is 0.63 %.
- In the α chain no contradiction with data is open now. X-01 (spin: the orbit would carry ħ, the measured spin is ħ/2) was resolved on 26.11 by the author's choice A-97: the inertia is split equally between the two ends of the lag exchange, so the orbit carries ħ/2. X-07 (hydrogen) was withdrawn under the data rule, and X-10 (a shortfall in 1/α) was lifted by the hadron data. All are kept in `meta.resolved`. Still open: the electron's magnetic anomaly 1.16·10⁻³. In cosmology one contradiction is open, X-13 (the Hubble tension): with the law of link density A-105 the relic angle needs a far rate h ≈ 0.685, while local candles give 0.735 (Cepheids) or 0.710 (red giants); see [`docs/cosmology_chain.md`](docs/cosmology_chain.md).

## How it is organised

- [`system/matrix.json`](system/matrix.json) — the matrix: author's postulates (А-), data (Д-), calculations (Р-), statements (В-), withdrawn statements (С-), contradictions (Х-). Every statement has a status: `вывод` (derived), `расчёт` (computed), `кандидат` (candidate), `данные` (data), `открыто` (open) or `противоречие` (contradiction). Errors are fixed at the source, and the previous text is kept in the `was` field.
- [`system/check_matrix.py`](system/check_matrix.py) checks the matrix. The main rule: derived and computed statements may never rest on a candidate.
- [`alpha/`](alpha) — calculations. The rules of each calculation (`*_rules_*.md`) are written down **before** the program is run; the program output is kept in `*_output.txt`.
- [`docs/article.md`](docs/article.md) — overview article; [`docs/predictions.md`](docs/predictions.md) — what can be tested.
- [`docs/alpha_chain.md`](docs/alpha_chain.md) — the 76 matrix records behind the α result and the α–G link, translated into English.
- [`docs/cosmology_chain.md`](docs/cosmology_chain.md) — the cosmology chain in English: redshift without expansion, the law of link density, the relic spots, the share of mass, the Hubble tension, antimatter, the 10¹²⁰ problem and the age test (A-90 … A-109, V-367 … V-402, V-423, V-424, X-13).
- [`docs/particles_chain.md`](docs/particles_chain.md) — points, the photon and the neutrino in English: what a photon is, why it has no rest mass, neutrino oscillations and mass ordering, why a neutrino passes through matter (A-18, A-19, A-109, V-46, V-47, V-416 … V-419).
- [`docs/hadrons_chain.md`](docs/hadrons_chain.md) — the hadrons chain in English: closure types, annihilation, the p̄ n prediction made before the data, colour factors, neutrinos on nucleons, anti-parts and the nucleon size (V-409 … V-415, V-420 … V-422).

**Language.** The matrix, calculation rules and program comments are in Russian, the author's working language; the records of the α and G chain are translated in `docs/alpha_chain.md`. Record IDs use Cyrillic letters: А = A (postulate), Д = D (data), Р = R (calculation), В = V (statement), С = S (withdrawn), Х = X (contradiction). The English texts write them in Latin letters (A-31 = А-31). Russian originals of the documents are in `docs/*_ru.md`.

## Running

Requires Python 3.10+ with `numpy` and `scipy`.

```
pip install -r requirements.txt
python system/check_matrix.py
python alpha/end_to_end.py                 # step, slope, edge, charges; fast q0 (2 clouds)
python alpha/q0_more_clouds.py             # q0 over 32 new clouds (R-113)
python alpha/hadron_dispersion.py <file>   # hadronic part from R(s) data (R-97)
python alpha/tail_gee.py <file>            # excess tail and narrow resonances (R-115)
```

`<file>` is the PDG compilation of R(s) for e⁺e⁻ → hadrons, `rpp2022-hadronicrpp_page1001.dat` (https://pdg.lbl.gov/2024/hadronic-xsections/). The result in the table (137.280) is assembled from these programs, as described in [`docs/alpha_chain.md`](docs/alpha_chain.md).

On Windows, add `-X utf8` if the output encoding is wrong.

## Related work

- H. Danielyan. *The Electron in the W-Structure: Spin, Magnetic Moment, Charge and Motion.* Version 1.1. Zenodo, 2026. https://doi.org/10.5281/zenodo.24221214 — code: https://github.com/gr2020155-code/w-structure-spin-1.0. Several entries of the matrix come from it (R-09, R-25, V-74, X-01, X-07); its positions P1–P18 are axioms A-… here under other numbers.
- H. Danielyan. *Black Holes in the W-Structure.* Zenodo, 2026. https://doi.org/10.5281/zenodo.24016340

## License

Code (`*.py`) is under the MIT license ([`LICENSE`](LICENSE)). Texts, the matrix and program outputs are under CC BY 4.0 ([`LICENSE-TEXT.md`](LICENSE-TEXT.md)): they may be reused provided the author is credited.
