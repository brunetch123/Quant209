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

That gives **41,376 detectable peak comparisons** across several analysts, two instruments and
four years of method drift.

Every run goes through the production entry point with production defaults. No parameter is
fitted per run.

## Headline agreement

`RPD = (auto − manual) / mean`, so it is symmetric and bounded at ±200%. Figures cover native
analytes whose *manual* value clears the reporting threshold (area ≥ 1000 counts·s,
mass ≥ 1 ng).

| | rows | median &#124;RPD&#124; | within ±25% | total ratio |
|---|---|---|---|---|
| Peak areas | 41,376 | **1.86%** | **91.8%** | 1.006 |
| Congener masses | 30,712 | **3.94%** | **87.7%** | 1.005 |
| Surrogate recoveries | 6,711 | 0.021 mean absolute error | | |

Areas and masses both meet the acceptance criteria (median ≤ 10%, ≥ 80% within ±25%). A
mass multiplies three measured quantities, so it inherits the area disagreement and adds to
it, and its coverage clears the gate by less.

**No run disagrees by more than 3.6% at the median.** The best is 0.61%. **No sample's
total native PCB mass differs from the analyst's by 2× or more**, apart from two trace-level
samples where the difference is a few nanograms.

**Every comparison is checked to be the same data file.** The labelled standards spiked into
every file act as a fingerprint: for the same file our areas and the analyst's agree to a
median 0.4%. A sheet is scored against a different file only on a decisive, mutual match.
Every verdict is written to `pairing.csv` (see
[`../METHODS.md`](../METHODS.md) section 6.1).

**The analyst's values are taken from their final workbook, on one basis.** Re-issued
workbooks are read in the order they were issued, not by how their names sort. And where a
workbook enters one member's calibration mass for a co-eluting group instead of the members'
sum, its group masses are put on the summed basis Quant 209 reports on, so every comparison is
like for like.

## Changes since 1.0.0

| | 1.0.0 | 1.1.0 | 1.2.0 | 1.3.0 | 1.4.0 | 1.4.1 |
|---|---|---|---|---|---|---|
| Peak areas: median &#124;RPD&#124;, within ±25% | 2.13%, 87.1% | 2.12%, 87.4% | 2.04%, 88.8% | 1.99%, 89.8% | 1.94%, 90.5% | 1.86%, 91.8% |
| Congener masses: median &#124;RPD&#124;, within ±25% | 4.61%, 78.7% | 4.60%, 79.0% | 4.24%, 84.1% | 4.14%, 85.5% | 4.03%, 86.5% | 3.94%, 87.7% |
| Total absolute mass error, detectable rows | 3,214k ng | 3,004k ng | 1,565k ng | 1,474k ng | 1,434k ng | 1,295k ng |
| Detectable rows reported at 1% or less of the analyst's value | 1,673 | 1,340 | 1,174 | 922 | 935 | 868 |

1.1.0's changes come from the named rules for congeners 83, 117, 90+101+113 and 118
([`../METHODS.md`](../METHODS.md) section 2.7), whose peaks carry the `forced` flag. 1.2.0
adds the displaced-peak rule (section 2.8). In a concentrated sample a dominant peak can drift
into a minor neighbour's retention-time window, and the neighbour was credited with all of it.
The rule gives 108 such peaks back, flagged `displaced`. 1.2.0 also reads a sample whose
internal standards are all diluted in step as a dilution rather than refusing it (section
7.3), which quantifies six more samples, and splits envelopes two comparable congeners merged
into by a two-peak fit, flagged `apportioned`. It no longer lets a calibration injection whose
retention times do not fit anchor other files, searches a file whose labelled standards moved
where they went (flagged `realigned`), reports a missing internal standard as an error with no
mass, and scores every analyst sheet against the file it came from. 1.3.0 splits the envelope
two neighbouring congeners share where a run's calibration does not separate them, flagged
`unresolved` (section 2.9), and settles near ties in the retention-time alignment by the
standard's composition (section 3.1). 1.4.0 keeps ghost peaks, the copy of a congener's peak
that losing HCl or Cl2 in the ion source leaves on a lower homolog's channel, out of
integrations, flagged `ghost` (section 2.10), and follows a file's drift along a curve
(section 3.3). 1.4.1 keeps calibration anchors off ghost copies and in line with the drift the
run's other transitions show (section 3.1), and lets a peak's window reach 0.25 min in every
part of a trace (section 2.3). Together with the ground-truth corrections above, the total
absolute mass error is 57% below 1.1.0's. The area median barely moves, because these act on few rows. The named
rules, and two details of the displaced-peak rule, were chosen by reading this same
comparison, so that share of the improvement is in-sample.

## Agreement against signal

| manual area (counts·s) | rows | median &#124;RPD&#124; | within ±25% |
|---|---|---|---|
| > 10⁶ | 2,394 | 0.5% | 95.4% |
| 10⁵ – 10⁶ | 6,321 | 0.9% | 96.7% |
| 10⁴ – 10⁵ | 13,076 | 1.5% | 94.3% |
| 10³ – 10⁴ | 19,585 | 3.0% | 88.1% |
| 10² – 10³ | 25,903 | 11.4% | 69.2% |
| < 10² | 52,084 | 107% | 17.9% |

The bottom two bands are below the reporting threshold. In that range, manual integration is
imprecise too, and results are censored against the limit of quantification. They are shown
for completeness, not as a claim.

## Agreement by homolog

| homolog | rows | median &#124;RPD&#124; | within ±25% |
|---|---|---|---|
| mono | 757 | 1.05% | 99.3% |
| di | 3,575 | 1.60% | 98.5% |
| tri | 5,740 | 1.96% | 97.3% |
| tetra | 7,730 | 2.20% | 89.6% |
| penta | 9,067 | 2.41% | 85.1% |
| hexa | 8,962 | 1.51% | 91.8% |
| hepta | 4,411 | 1.35% | 95.6% |
| octa | 975 | 2.61% | 92.4% |
| nona | 147 | 5.03% | 89.8% |
| deca | 12 | 9.17% | 91.7% |

Penta is the weakest well-populated homolog. Its transition carries 37 analytes, the most
crowded window in the method. Crowding is where automated and manual integration differ most.

## Where the disagreement is

868 detectable rows (2.1% of comparisons) are reported at 1% or less of the analyst's
value. They are two different failures.

- **Riders on a larger neighbour's leading edge.** The 822 rows below 10⁶ counts·s are many
  but small: together they carry 0.2% of the analyst's detectable area. Four congeners (115,
  88, 98 and 143) supply 428 of them. The separation there is genuinely ambiguous.
- **Dominant peaks in concentrated samples.** The 46 rows above 10⁶ counts·s carry 3.6% of
  the analyst's detectable area, down from 8.7% before the displaced-peak rule. The peak's
  top moves off its retention-time anchor and the area goes to the wrong congener, while each
  file's total area still agrees with the analyst's. The rule gives back the clear cases.
  The ones it leaves are described in [`../METHODS.md`](../METHODS.md) section 7.2.

Nearly all of them carry the `absent` flag. Consequently, they appear in the QA/QC table
rather than passing silently. [`../METHODS.md`](../METHODS.md) section 7 explains both
mechanisms and the other known limitations, including the one sample refused because its
injection has no internal standard, and how analyst sheets are paired with data files.

## Figures

`tools/validation_figures.py` draws the validation figures from the pooled comparison, as SVG:
agreement by integration method and flag, the area and total-mass scatters, Aroclor
fingerprint agreement, detection agreement, and breakdowns of where and why the two
integrations differ. Per-sample totals sum native congeners only, never the spiked labelled
standards. The Aroclor figures need a congener-level reference table, given with `--aroclor`. The figures name run folders and plot per-sample results, so none is committed.

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
.venv\Scripts\python.exe tools\validation_figures.py --data docs\validation\master --aroclor <reference table>
```

The first runs the pipeline over every run folder present and writes the pooled comparison.
The second turns it into a full report, `docs/validation/MASTER-REPORT.md`, from which this
page is taken. The third draws the figures into `docs/images/validation/`. None of the three
outputs is committed. They are large, regenerable and keyed by sample.
