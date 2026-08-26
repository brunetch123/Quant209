# Validation

Quant 209 is compared against the analysts' own manual integration and quantification. Their
workbooks are the reference. They cover the same peaks, integrated by hand over four years of
routine work.

## What was compared

The laboratory's whole archive holds **64 run folders**. Of those, **49 hold analyst
quantification that joins to the acquired files and can be scored**. Five hold no analyst
quantification at all, being either an instrument sequence log or a workbook superseded by a
re-run. One holds no calibration injection, so nothing can be derived from it. The remaining
nine process cleanly and do have analyst workbooks, but none of their sample columns match the
sample identifiers derived from the acquisition filenames, so no row can be joined. All
fifteen unscored folders are reported and skipped rather than silently dropped. The figures
below say nothing about them.

That gives **41,240 detectable peak comparisons** across several analysts, two instruments and
four years of method drift.

Every run goes through the production entry point with production defaults. No parameter is
fitted per run.

## Headline agreement

`RPD = (auto − manual) / mean`, so it is symmetric and bounded at ±200%. Figures cover native
analytes whose *manual* value clears the reporting threshold (area ≥ 1000 counts·s,
mass ≥ 1 ng).

| | rows | median &#124;RPD&#124; | within ±25% | total ratio |
|---|---|---|---|---|
| Peak areas | 41,240 | **2.13%** | **87.1%** | 0.975 |
| Congener masses | 30,549 | **4.61%** | **78.7%** | 0.991 |
| Surrogate recoveries | 6,442 | 0.030 mean absolute error | | |

Areas meet the acceptance criteria (median ≤ 10%, ≥ 80% within ±25%). Masses meet the median
criterion. They fall 1.3 points short of the coverage one. A mass multiplies three measured
quantities, so it inherits the area disagreement and adds to it.

**No run disagrees by more than 7.4% at the median.** The best is 0.69%.

## Agreement against signal

| manual area (counts·s) | rows | median &#124;RPD&#124; | within ±25% |
|---|---|---|---|
| > 10⁶ | 2,391 | 0.7% | 88.1% |
| 10⁵ – 10⁶ | 6,316 | 1.0% | 91.2% |
| 10⁴ – 10⁵ | 13,017 | 1.7% | 89.4% |
| 10³ – 10⁴ | 19,516 | 3.3% | 84.0% |
| 10² – 10³ | 25,752 | 12.4% | 66.4% |
| < 10² | 51,333 | 109% | 17.7% |

The bottom two bands are below the reporting threshold. In that range, manual integration is
imprecise too, and results are censored against the limit of quantification. They are shown
for completeness, not as a claim.

## Agreement by homolog

| homolog | rows | median &#124;RPD&#124; | within ±25% |
|---|---|---|---|
| mono | 751 | 1.10% | 98.0% |
| di | 3,551 | 1.63% | 97.4% |
| tri | 5,717 | 1.98% | 95.7% |
| tetra | 7,712 | 2.21% | 87.4% |
| penta | 9,029 | 3.01% | 78.8% |
| hexa | 8,930 | 1.97% | 84.8% |
| hepta | 4,416 | 1.77% | 89.2% |
| octa | 975 | 3.52% | 74.6% |
| nona | 147 | 5.19% | 87.8% |
| deca | 12 | 9.17% | 91.7% |

Penta is the weakest well-populated homolog. Its transition carries 37 analytes, the most
crowded window in the method. Crowding is where automated and manual integration differ most.

## Where the disagreement is

1,673 rows (4.1% of comparisons but 12.1% of detectable area) are reported at 1% or less of the
analyst's value. Five congeners (83, 117, 115, 88, 143) supply 675 of them. Each elutes on the
leading edge of a much larger neighbour, where the separation is genuinely ambiguous. Nearly
all carry the `absent` flag. Consequently, they appear in the QA/QC table rather than passing
silently. [`../METHODS.md`](../METHODS.md) section 7 explains the mechanism and the other known
limitations, including a probable sample swap in one run's ground truth.

## Per-run reports

`quant209 validate <run folder>` writes a full report for a single run. It contains
`summary.md`, per-analyte agreement for areas, masses, recoveries and response factors,
`bias_diagnosis.csv` and three figures. `bias_diagnosis.csv` lists every detectable native mass
row outside the band, with its integration flags, boundary method and probable cause. The
Export step in the application writes the same report.

No example is committed here. The reports name samples, and the laboratory's sample
identifiers are not published.

## Reproducing this

```
.venv\Scripts\python.exe tools\master_bench.py  --out docs\validation\master
.venv\Scripts\python.exe tools\master_report.py --data docs\validation\master
```

The first runs the pipeline over every run folder present and writes the pooled comparison.
The second turns it into this report. The pooled CSVs are not committed. They are large,
regenerable and keyed by sample.
