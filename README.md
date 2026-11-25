# Структура W

Попытка вывести постоянную тонкой структуры α и устройство частиц из небольшого набора положений о пространстве. Используются только положения и экспериментальные данные, без готовых постулатов стандартной физики.

**Автор положений:** Hrachya Danielyan
**Расчёты и проверки:** выполнены с помощью Claude (Anthropic) по правилам автора.

*English summary below.*

## Главный результат

| Величина | Структура W | Измерено |
|---|---|---|
| 1/α вдали от заряда | **136,86 ± 0,14** | 137,036 (−0,13 %) |
| 1/α на расстоянии Z-бозона | 129,54 ± 0,14 | 129,72 ± 0,05 |
| Наклон роста α | 2/(3π) | 2/(3π) |
| Заряды частиц поколения | 0, 1/3, 2/3, 1 | ν, d, u, e |
| Число «цветов» | 1, 3, 3, 1 | лептоны 1, кварки 3 |

Свободных подгоночных параметров нет. Весь путь от положений до числа считает одна программа: [`alpha/end_to_end.py`](alpha/end_to_end.py).

## Честный статус

- **Все совпадения получены, когда ответ был уже известен.** Поэтому они имеют статус **кандидат**, а не доказательства. Нового предсказания, подтверждённого опытом, пока нет.
- В расчёте α используются **измеренные** массы частиц, M_Z и m_t. Массы и слабая сила в структуре пока калибровка (положение А-14), а не вывод.
- Результат держится на трёх непроверенных звеньях. Они перечислены в [`docs/predictions_ru.md`](docs/predictions_ru.md) и в конце вывода программы:
  1. W-бозон не гасит заряд (В-294);
  2. «остаток пары — это q₀» (В-295);
  3. перевод «энергия → расстояние» через γ.
- Единственное проверяемое отличие от стандартной физики: выше ~160 ГэВ α растёт быстрее. Сегодняшние данные недостаточно точны, чтобы это проверить.
- В матрице открыто записаны два нерешённых противоречия с данными: Х-01 (спин) и Х-07 (водород: линия 21 см и член Дарвина).

## Как устроено

- [`system/matrix.json`](system/matrix.json) — матрица: положения автора (А-), данные (Д-), расчёты (Р-), утверждения (В-), снятые утверждения (С-), противоречия (Х-). У каждого утверждения есть статус: `вывод`, `расчёт`, `кандидат`, `данные`, `открыто` или `противоречие`. Ошибки исправляются у источника, прежний текст сохраняется в поле `was`.
- [`system/check_matrix.py`](system/check_matrix.py) проверяет матрицу. Главное правило: `вывод` и `расчёт` не могут опираться на `кандидат`.
- [`alpha/`](alpha) — расчёты. Правила счёта (`*_rules_*.md`) записываются **до** запуска программы, вывод программы хранится в `*_output.txt`.
- [`docs/article_ru.md`](docs/article_ru.md) — обзорная статья, [`docs/predictions_ru.md`](docs/predictions_ru.md) — что можно проверить.

## Запуск

Нужен Python 3.10+ и пакеты `numpy`, `scipy`.

```
pip install -r requirements.txt
python system/check_matrix.py
python alpha/end_to_end.py          # быстро: 2 облака, точность ±0,35
python alpha/end_to_end.py --full   # 8 облаков, ±0,14 (результат в таблице)
```

В Windows при неверной кодировке добавьте `-X utf8`.

## Лицензия

Код (`*.py`) распространяется по лицензии MIT ([`LICENSE`](LICENSE)). Тексты, матрица и выводы — по лицензии CC BY 4.0 ([`LICENSE-TEXT.md`](LICENSE-TEXT.md)): их можно использовать при условии, что указан автор.

---

## English summary

**Structure W** is an attempt to derive the fine-structure constant α and the pattern of particle charges from a small set of postulates about space, using only those postulates and experimental data.

- Result: 1/α far from the charge = **136.86 ± 0.14** vs measured 137.036 (−0.13 %); at the Z distance 129.54 ± 0.14 vs 129.72 ± 0.05. No fitted parameters. One script, `alpha/end_to_end.py`, goes from the postulates to the number.
- The same rule gives quark charges 1/3, 2/3, colour multiplicities 1, 3, 3, 1, and the running slope 2/(3π).
- **Status: candidate.** All agreements were found with the answer already known; no new prediction has yet been confirmed. Measured particle masses are used as input. Three links are unverified (the W boson does not screen; "pair remainder = q₀"; the energy-to-distance conversion). Two open contradictions with data are recorded in the matrix (spin, hydrogen hyperfine).
- Method: every claim in `system/matrix.json` carries a status, and `check_matrix.py` enforces that derivations never rest on candidates. Calculation rules are written down **before** each run.
- Documents are in Russian.
