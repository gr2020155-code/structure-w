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
| 1/α far from the charge | **136.86 ± 0.14** | 137.036 (−0.13 %) |
| 1/α at the Z-boson distance | 129.54 ± 0.14 | 129.72 ± 0.05 |
| Slope of the running of α | 2/(3π) | 2/(3π) |
| Charges within one generation | 0, 1/3, 2/3, 1 | ν, d, u, e |
| Number of "colours" | 1, 3, 3, 1 | leptons 1, quarks 3 |

There are no fitted parameters. One program, [`alpha/end_to_end.py`](alpha/end_to_end.py), goes all the way from the postulates to the number.

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

- **All agreements were found with the answer already known.** So they are **candidates**, not proofs. No new prediction has been confirmed by experiment yet.
- The α calculation uses **measured** particle masses, including M_Z and m_t. In the structure, masses and the weak force are still calibration (postulate A-14), not derivation.
- The result rests on three unverified links, listed in [`docs/predictions.md`](docs/predictions.md) and at the end of the program output:
  1. the W boson does not screen charge (V-294);
  2. "the pair remainder is q₀" (V-295);
  3. the energy → distance conversion through γ.
- The only testable difference from standard physics: above ~160 GeV α grows faster. Current data are not precise enough to test it.
- Two unresolved contradictions with data are recorded openly in the matrix: X-01 (spin) and X-07 (hydrogen: the 21 cm line and the Darwin term).

## How it is organised

- [`system/matrix.json`](system/matrix.json) — the matrix: author's postulates (А-), data (Д-), calculations (Р-), statements (В-), withdrawn statements (С-), contradictions (Х-). Every statement has a status: `вывод` (derived), `расчёт` (computed), `кандидат` (candidate), `данные` (data), `открыто` (open) or `противоречие` (contradiction). Errors are fixed at the source, and the previous text is kept in the `was` field.
- [`system/check_matrix.py`](system/check_matrix.py) checks the matrix. The main rule: derived and computed statements may never rest on a candidate.
- [`alpha/`](alpha) — calculations. The rules of each calculation (`*_rules_*.md`) are written down **before** the program is run; the program output is kept in `*_output.txt`.
- [`docs/article.md`](docs/article.md) — overview article; [`docs/predictions.md`](docs/predictions.md) — what can be tested.
- [`docs/alpha_chain.md`](docs/alpha_chain.md) — the 51 matrix records behind the α result and the α–G link, translated into English.

**Language.** The matrix, calculation rules and program comments are in Russian, the author's working language; the records of the α and G chain are translated in `docs/alpha_chain.md`. Record IDs use Cyrillic letters: А = A (postulate), Д = D (data), Р = R (calculation), В = V (statement), С = S (withdrawn), Х = X (contradiction). The English texts write them in Latin letters (A-31 = А-31). Russian originals of the documents are in `docs/*_ru.md`.

## Running

Requires Python 3.10+ with `numpy` and `scipy`.

```
pip install -r requirements.txt
python system/check_matrix.py
python alpha/end_to_end.py          # fast: 2 clouds, precision ±0.35
python alpha/end_to_end.py --full   # 8 clouds, ±0.14 (the result in the table)
```

On Windows, add `-X utf8` if the output encoding is wrong.

## Related work

- H. Danielyan. *The Electron in the W-Structure: Spin, Magnetic Moment, Charge and Motion.* Version 1.1. Zenodo, 2026. https://doi.org/10.5281/zenodo.24221214 — code: https://github.com/gr2020155-code/w-structure-spin-1.0. Several entries of the matrix come from it (R-09, R-25, V-74, X-01, X-07); its positions P1–P18 are axioms A-… here under other numbers.
- H. Danielyan. *Black Holes in the W-Structure.* Zenodo, 2026. https://doi.org/10.5281/zenodo.24016340

## License

Code (`*.py`) is under the MIT license ([`LICENSE`](LICENSE)). Texts, the matrix and program outputs are under CC BY 4.0 ([`LICENSE-TEXT.md`](LICENSE-TEXT.md)): they may be reused provided the author is credited.
