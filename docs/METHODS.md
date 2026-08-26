<picture>
  <source media="(prefers-color-scheme: dark)" srcset="images/logo-dark.png">
  <img src="images/logo.png" alt="Quant 209" height="44">
</picture>

# Methods

**What Quant 209 does to a chromatogram, and why.**

Written for an analyst deciding whether to trust a number this software produced, and for a
reviewer who wants to know what was automated and what was assumed.

Every figure in this document traces to a named source, either the comparison report under
`docs/validation/` or the module that implements the rule.

---

## 1. The method under automation

One analysis run is a single GC-MS/MS sequence: a calibration standard injected several times,
one or more reference (recovery) standards, laboratory and batch blanks, and the field samples.
All of it is acquired in one batch on one column build.

**21 SRM transitions, 185 canonical analytes.** The acquisition method monitors 21
precursor→product transitions, one per chlorination level for the natives plus the labelled
standards' isotope-shifted channels. Each transition carries an ordered list of expected
analytes. The full table (`quant209/config/canonical_peaks.csv`) holds **185 rows: 173 native
congeners and inseparable co-elution groups, and 12 labelled standards**. Each row carries its
transition, its homolog, its order within the transition, its role, which internal standard
normalises it, its certified mass in the calibration standard, and its canonical retention
time.

A "native analyte" is not always one congener. Where the method cannot separate two or three
congeners (`90+101+113`, `85+116`, `40+71`, `87+125`), the canonical table carries the group as
one row. It is quantified as one quantity, exactly as the analyst reports it. This is a
property of the chromatography, not a simplification: the peaks are one peak.

**Standards.** Two internal standards are spiked into every file immediately before injection,
at **10.4125 ng** each. d-PCB30 normalises the low-chlorination homologs and PCB204 the high
ones. The assignment is recorded per analyte in the canonical table's `is_norm` column. Ten
`13C`-labelled surrogates are spiked **pre-extraction** at **9.84 ng** each, one per homolog.
Each homolog's congeners are corrected for extraction and cleanup losses using its own
surrogate (`quant209/core/standards.py`, `HOMOLOG_SURROGATE`). Both spike amounts were
reconstructed from laboratory records rather than read from a certificate. The distinction
between the two kinds of standard matters. An
internal standard measures *injection*, and a surrogate measures *recovery*. Only the latter is
a correction.

---

## 2. Integration

Implemented in `quant209/core/peaks.py` (the primitives) and `quant209/core/integrate.py` (the
per-file decision tree). The tree is calibrated against the analyst's manual integration.

### 2.1 Retention-time anchoring

Retention times drift between column builds, maintenance events and instruments. The canonical
table's retention times are therefore a **starting point, not an address**. The calibration
standard contains every congener. It is therefore the one injection where every analyte's peak
is guaranteed present. Each analyte is located near its canonical retention time,
refined to the local apex, and that observed retention time becomes the anchor used for
**every other file in the run**. Every sample is therefore integrated against retention
times measured on its own sequence.

Anchors come from the **first** calibration injection only, never averaged across
replicates. An anchor is a position. Averaging positions across injections separated by
a whole sequence would blur the drift the anchors exist to track.

### 2.2 Territories, apex refinement, and the presence test

Within a transition, each analyte owns a **territory** bounded by the midpoints between
its anchor and its neighbours' anchors. The apex is refined inside a ±`rt_window` =
0.10 min search around the anchor, capped by `max_reach` = 0.20 min.

The **presence test** then asks whether that apex is a resolved peak at all. An analyte is
accepted only when its apex prominence clears **both** an absolute floor (`min_prominence`
= 150 counts) **and** a noise-scaled floor (`presence_sn` = 4.0 × the trace's global
noise). An absolute floor alone misclassifies quiet and noisy channels differently. A
noise-scaled floor alone accepts a bump on a flat baseline. This two-sided form was arrived
at by comparison against analyst areas, not from first principles.

An analyte that fails carries the **`absent`** flag. Its area is still reported, as a bounded
estimate of whatever occupies its territory. This is what the analyst does with a shoulder they
chose to split off a dominant neighbour. Suppressing it would silently change
what the report means. However, `absent` is the flag every downstream consumer gates on.

### 2.3 Boundaries

Two methods, recorded per peak in `PeakResult.method`:

* **`valley_drop`**. At least one neighbouring analyte was itself a resolved peak, so
  the shared boundary is the valley between them (a perpendicular drop). Adjacent peaks'
  boundaries tile exactly: each peak's end is the next one's start.
* **`rt_window`**. No resolved neighbour bounded this peak. An isolated peak is bounded by
  the local minima within `max_halfwidth` = 0.25 min of its apex. A peak beside an *absent*
  neighbour extends through it to its own baseline edge. A shoulder or absent analyte is
  bounded valley-to-valley inside its own midpoint territory.

**Co-elution: split versus summed.** Where two congeners are chromatographically resolved
they are integrated separately by perpendicular drop, and each is reported. Where the method
does not resolve them at all, the canonical table carries them as one summed row and one
peak is integrated. The decision is made once, in the table, from the chromatography. It
is never made per file or per run. The same congener pair cannot be split in one sample and
summed in the next.

### 2.4 Cluster baselines on crowded windows: the measured justification

Implemented in `_cluster_spans`, in `quant209/core/integrate.py`.

**The problem.** Perpendicular-drop boundaries tile edge to edge. Drawing each slice's
baseline as the chord between *its own* two boundaries rests that baseline on the elevated
valleys either side of it. Where several congeners sit on one unresolved bed, the summed
area of the slices misses the bed entirely.

**The measurement.** On one calibration injection's crowded penta window,
congener by congener against the analyst's own `cal` sheet:

| analyte | Quant 209 (before) | analyst area | ratio | after the cluster baseline |
|---|---|---|---|---|
| `93+100` | 163,141 | 261,319 | 0.62 | 248,448 (−4.9%) |
| `102` | 33,282 | 158,301 | 0.21 | 153,342 (−3.1%) |
| `98` | 11,115 | 105,688 | 0.11 | 109,377 (+3.5%) |
| `88` | 17,662 | 105,095 | 0.17 | 99,006 (−5.8%) |
| `91` | 46,117 | 138,007 | 0.33 | 144,254 (+4.5%) |
| **sum** | **271,316** | **768,411** | **0.35** | **754,427 (−1.8%)** |

Two things fall out of this comparison. First, **the retention times and the boundaries were
already right**: every apex was within 0.03 min of the analyst's, and the boundaries tiled
exactly. This was not a peak-placement problem. Second, **the window as a whole lost 65% of its
signal**. It was not a redistribution between neighbours either. It was the baseline.

That also explains an asymmetry: sample areas agree to about 2.4% while calibration areas
for the same congeners were out by a factor of ten. In a *sample* these penta congeners are
mostly small or absent. No run of consecutive resolved peaks forms, and each slice's chord
baseline *is* the true baseline. The error is specific to the crowded standard. A calibration
injection is exactly that. It contains all 209 congeners at comparable concentration.

**What shipped.** Consecutive analytes are grouped into a cluster. Every member is then
integrated above **one** baseline spanning the whole cluster, keeping its own valley-drop
boundaries. Two neighbours join a cluster only when both:

1. are **present**. An analyte the presence test rejected breaks the run.
2. **do not return to base** between their apexes. The valley between them stands more than
   `overlap_frac` = 0.05 of the smaller apex's height above the pair's outer boundaries.

**Condition 1 is load-bearing.** Without it, a cluster baseline over-assigns a dominant peak's
underlying area to a small shoulder, whereas the analyst skims shoulders. The shoulder case
*is* the absent-analyte case. The present-only restriction leaves it untouched. A test pins
that a cluster never jumps an absent analyte.

The rule is self-limiting. Two peaks whose valley reaches the chromatographic floor score ~0 on
condition 2 and are left alone whatever `overlap_frac` is. The parameter is therefore only a
tolerance on "did it reach the floor". Clusters are short in practice. On one run's `Cal01`
there are 26 clusters, 21 of two members and 5 of three, and no cluster longer than three
members occurs anywhere in either run's calibration injections. That bounds the span of a
single baseline chord, keeping the straight-line assumption valid.

**`overlap_frac` = 0.05 is not a value the data pins.** A sweep over 0.02 to 0.15 moves
the >±25% outlier count on five of six (injection, workbook) columns not at all, while
the merge set shrinks monotonically (37 → 31 → 22 → 20 merged pairs on that `Cal01`)
and 42 analytes' areas differ between the endpoints by as much as −22.2%. These are
genuinely different integrations that happen to agree on how many analytes land outside
a 25% band. 0.05 was chosen over 0.02 (which scores *better* on median |RPD| everywhere)
because taking the endpoint that most flatters a metric would be fitting the benchmark.

**Result.** Response-factor agreement against the analyst improved from 2.6% to **1.9%** median
|RPD| on the first run and 4.0% to **2.1%** on the second, against an analyst reproducibility
floor of 1.1% and 1.3%. That is within a factor of two of the best agreement the ground truth
can support.

### 2.5 Units, and what a flag means

Areas are **counts·seconds**: trapezoidal integration of the raw (unsmoothed) signal above the
straight baseline joining the two boundary points, with the retention-time axis in minutes
converted to seconds. Smoothing (Savitzky-Golay, window 5, order 2) is used for *detection*
only. The raw signal is what is integrated. A reported area is therefore a measurement of the
acquired trace, not of a filtered version of it.

Four flags. Each says something different:

| flag | meaning |
|---|---|
| `absent` | the presence test rejected this analyte. The area is a bounded territory estimate, not a peak. **This is the detection verdict.** |
| `low_sn` | a detected peak whose baseline-corrected apex height is below `min_sn` = 3.0 × trace noise: real, but near the detection limit, so the area is imprecise. |
| `narrow` | fewer than 3 samples between the boundaries, so the area is dominated by where the boundaries landed rather than by the peak shape. |
| `clustered` | this area was measured above a cluster baseline rather than its own chord. **Not a confidence verdict.** It tells anyone recomputing the area by hand which baseline to use. |

---

## 3. Per-run retention-time alignment and the run-fit guard

Implemented in `quant209/core/calibrate.py`.

### 3.1 Why alignment exists

Anchoring searches ±0.10 min around each canonical retention time. That window suits column
drift within one GC method. In contrast, it is far too narrow when a run was acquired under a
different oven program. On one such run, congener 1 elutes at 8.34 min against a canonical
19.56, an 11-minute displacement. Without alignment, every analyte anchors onto whatever noise
sits in its window. The run then produces a plausible, wrong report.

Before anchoring, Quant 209 does what an analyst does by eye when the method changes. It lays
the calibration injection's own peak pattern against the canonical retention-time sequence.
Then, channel by channel, it fits the straight line that maps one onto the other. The line is a
correspondence-finder, not a correction. Once it says which observed peak is which congener,
that peak's own retention time is the search centre.

Three constraints bound the fit:

* Channels with fewer than 3 canonical analytes do not get their own map. Instead, they
  inherit one, because a single labelled standard on a quiet channel could be matched to any
  peak at all.
* The slope is bounded to 0.5 to 2.0. A method change re-times a chromatogram. It does not
  compress or stretch it by more than a factor of two. The extremes measured on this dataset
  are 0.791 and 1.494. Bounding the slope stops the search from "explaining" a channel with a
  degenerate near-flat map that piles every congener onto one peak.
* A fitted map displacing by less than `IDENTITY_TOL` = 0.05 min is reported as the identity
  and not applied. Fitting a line to apexes always returns something slightly off identity.
  On the run the canonical table was extracted from, the fitted maps displace by at most
  0.008 min, which is apex noise, not a method difference. Declaring identity there keeps a
  same-method run's anchors exactly what they would be without alignment.

### 3.2 The run-fit guard

Both failures above are silent by construction. Integration always returns a number, so a run
scored against the wrong retention times produces a plausible-looking report rather than an
error. This is the most dangerous failure mode in the pipeline, because it does not look like a
failure.

The guard (`rt_fit_quality`) checks the calibration injection against the table it was
anchored with. On a run whose retention-time table fits, nearly every analyte in the
calibration standard is a resolved peak. A low found-fraction there is therefore direct
evidence that the table does not describe this run. Measured: the canonical-method run
passes at native 170/173 (98%) and labelled 12/12 (100%). The different-oven-program run
passes after alignment at native 167/173 (97%) and labelled 12/12 (100%).

A failing guard does not stop the run. Every number is still computed and the run's status
becomes `FAILED-QC`. The guard's reason is attached to the result and carried verbatim into the
GUI's red banner. An analyst looking at a flagged run is the person who can tell a GC
method change from a bad injection, and needs the numbers to do it.

---

## 4. Quantification

Implemented in `quant209/core/quantify.py`. Orchestration is in `quant209/core/pipeline.py`.

### 4.1 Response factors

For each analyte, from the calibration injection:

    RRF = (cal_mass / cal_area) / (IS_cal_mass / IS_cal_area)

where the internal standard is the one the canonical table assigns to that analyte. Response
factors are averaged over every usable calibration injection (four on the run below). An RRF is
a ratio of responses, so averaging several measurements of it is more reliable than taking one.
A congener mis-integrated in one injection is still quantified from the others. This is the
opposite of the anchor rule (section 2.1), and deliberately so. An anchor is a position. An RRF
is a magnitude.

A calibration guard refuses a response factor rather than issuing a wrong one. An analyte is
refused when it has no calibration peak, no certified mass, no usable internal-standard
response or zero area. It is also refused when the presence test rejected its calibration peak,
or when its calibration area falls below a floor. That floor is the larger of an absolute
minimum and 2% of the median positive native calibration area. A refused analyte gets no mass
in any file of the run. The reason is carried on the result and shown per analyte in the GUI's
Calibration tab. An empty cell and an absent cell mean different things.

### 4.2 Mass

    Mass (ng) = RRF × Area × (IS_spike / IS_area)

evaluated per file, using that analyte's own internal standard (d-PCB30 or PCB204, as the
canonical table assigns). Both are spiked at 10.4125 ng in this method. That is why one
constant serves. A method that spiked them differently would need a per-standard lookup.

The presence verdict is deliberately not applied here. The analyst integrates and quantifies a
shoulder they chose to split off a dominant neighbour, and so does this chain. The `detected`
verdict travels beside the mass, in its own frame, instead of silently zeroing it.

### 4.3 Per-batch reference bracketing

Recovery is measured against a reference injection. Each file is paired with the reference
injection that most recently preceded it in acquisition order. This is what a bracketed
sequence means, and what the analyst actually did.

The evidence comes from a three-batch run in the archive:

* The instrument sequence log runs `Cal01, Ref04, [8 files], Cal02, Ref05, [7 files], Cal03,
  Ref06, [7 files], Cal04`, three analytical batches, each opened by its own calibration and
  reference injection.
* The analyst kept one workbook per batch, and each workbook's `ref` sheet holds that batch's
  own reference. The first workbook's surrogate masses are `Ref04`'s, the second's are
  `Ref05`'s, and the third's are `Ref06`'s.
* Nearest-preceding bracketing reproduces that pairing for all 22 quantified files. It also
  moves the run's recovery error against the analyst from a mean absolute 0.066 to 0.023.

Nearest-preceding is used rather than nearest-in-time, because one lab blank was acquired 14 h
after `Ref04` and 7 h before `Ref05`. Nearest-in-time would pair it with `Ref05`, where the
analyst paired it with `Ref04`. A bracketing standard opens a batch. It does not retroactively
claim the batch before it. A file with no preceding reference falls back to the earliest one
after it, the only reference that can be said to bracket it.

Recovery is measured against one reference injection, never an average. The reference defines
100% recovery, so it is a denominator rather than a measurement to be improved. Averaging three
would make each sample's recovery depend on injections it was never compared to.

Per homolog:

    recovery = surrogate mass in this file / surrogate mass in its bracketing reference

### 4.4 Recovery correction (EPA-1668)

    corrected mass = mass / recovery,  applied only when 0 < recovery < 1

One-sided by design, following EPA-1668. Correcting for extraction losses is standard practice.
Deflating a mass because a surrogate over-recovered is not. An apparent gain is far more often
an integration or spiking artefact than a real one. A recovery of exactly 1, above 1, at or
below 0, or missing altogether leaves the mass untouched. Labelled standards carry no homolog
and are never corrected. This is why an internal standard still reads exactly its own spike
amount after this step. That provides a simple check on the mass chain.

The three stages are kept as three separate frames (`masses`, `masses_corrected`,
`masses_censored`), so an analyst asking "what did the recovery correction do to this
congener?" subtracts two columns instead of re-running the pipeline.

---

## 5. QA/QC

Implemented in `quant209/core/qaqc.py`.

### 5.1 Limit of quantification

    LOQ = mean + k·SD  of the recovery-corrected blank masses,  k = 3.0

This is computed per analyte, over the run's usable blanks. `k = 3` is roughly a one-sided 99%
interval. The limit is computed from corrected masses because it has to live in the same units
as the numbers it censors. A blank that reported no mass for an analyte contributes 0, not a
gap. The analyte was looked for and nothing was quantifiable. That is evidence of a low blank
level, not an absence of evidence. Dropping the cell would raise the LOQ of every
analyte that is clean in most blanks and dirty in one.

Fewer than two usable blanks yields no LOQ and censors nothing. With one blank there is no
spread to measure, so "mean + 3SD" is that blank's own mass. Censoring against it would delete
every sample result at or below a single injection's background. A run in that state says so on
the result and in the GUI.

### 5.2 What censoring applies to

Censoring sets an at-or-below-LOQ cell to 0. The comparison is `≤`: a mass exactly at the limit
has not cleared it. Scope is deliberately narrow:

* **Native analytes only.** A labelled standard is spiked at a known amount into every file,
  including the blanks, so its blank masses are identical by construction. This gives zero
  spread and an "LOQ" equal to the spike itself, so `mass ≤ LOQ` is true everywhere. Censoring
  the full frame zeroed 177 of 180 labelled-standard cells on one run. That deletes the one
  column a reviewer uses to check the mass chain.
* **Reported files only: samples and blanks.** Calibration and reference injections are
  quantified (they are what the response factors and the recovery denominators are made of),
  but they are not results.
* An analyte with no LOQ is left alone rather than censored against an assumed limit of zero.
  A NaN stays NaN. "Not measured" and "measured, and below the limit" are different report
  entries and must not be conflated.

### 5.3 Recovery band

Surrogate recoveries outside 50-150% are flagged, per homolog and per file. The band is drawn
on the GUI's recovery heatmap. Cells outside it carry a border as well as a darker fill. Colour
alone is not an accessible encoding, so the band stays legible to a colourblind reader and in a
black-and-white printout.

### 5.4 Why `sn` must never be used as a detection test

For an analyte the presence test rejected, the integrated region is the analyte's "RT
territory" rather than a region around its apex. As a result, the apex may lie outside the
integrated bounds. The reported height then extrapolates the region's straight baseline out to
that apex, and can come out large and positive. Absent peaks therefore routinely report
signal-to-noise far above any sensible threshold: values up to ~4.9 × 10⁴ have been observed on
a calibration injection.

An absent peak consequently does not pick up a `low_sn` flag. Any rule of the form "S/N > 3
means the analyte is present" would therefore confidently assert the presence of exactly the
analytes that are not there.

The contract is therefore to gate on `detected` first, then use `sn` to grade the peak.
`PeakResult.sn` is meaningful only when `PeakResult.detected` is true. This is stated at the
field, the property and the flag. The pipeline also exposes a `detected` frame beside `areas`,
so no consumer has to infer presence from a magnitude. In the GUI's peak table, an undetected
analyte's S/N cell is left blank rather than showing `n/a` or a number, with a tooltip
explaining why.

---

## 6. Validation

### 6.1 Methodology

The laboratory's archive is the test set. Every analysis run whose analyst workbook can be read
is scored against that analyst's own manual integration and quantification. The comparison
spans 49 runs and 41,240 detectable peak comparisons, across four years, several analysts and
two instruments.

The comparison runs through the production entry point (`core.pipeline.run_analysis`) with
production defaults and no per-run tuning. The same code path and the same parameters apply to
every run.

    RPD = (auto - manual) / ((auto + manual) / 2)

Headline metrics cover native analytes whose manual value clears a detectability threshold
(area ≥ 1000 counts·s, mass ≥ 1 ng). Labelled standards and sub-threshold rows are reported
separately rather than folded in. Near the detection limit, the analyst's own integration is
imprecise, and the result is censored anyway. Including those rows would let noise dominate a
number meant to measure the method. Recoveries are scored as absolute error in recovery
fraction, which is a ratio already.

**Acceptance.** Median |RPD| ≤ 10% and ≥ 80% of rows within ±25%, at area level and at mass
level. Surrogate-recovery mean absolute error ≤ 0.05. The ±25% band is what this
laboratory's practice treats as agreement for a trace congener measurement. The two area tests
are deliberately paired. A median alone tolerates a handful of large misses, and a coverage
figure alone tolerates a systematic small bias.

Response factors are reported but not gated. The right threshold for one is the analyst's own
cross-workbook reproducibility, a property of the ground truth rather than of the software.
Corrected and censored masses are likewise not scored. The corrected column is reproducible
from two numbers that are scored, and censoring is a reporting policy, not a measurement.

### 6.2 Results

| | rows | median &#124;RPD&#124; | within ±25% | total ratio |
|---|---|---|---|---|
| Peak areas | 41,240 | **2.13%** | **87.1%** | 0.975 |
| Congener masses | 30,549 | **4.61%** | **78.7%** | 0.991 |
| Surrogate recoveries | 6,442 | 0.030 mean absolute error | | |

Areas clear both gates. Masses clear the median gate and sit 1.3 points below the coverage
gate, for the reason given in section 7.2.

No run in the archive disagrees by more than 7.4% at the median. The best sits at 0.69%.
Agreement tracks signal:

| manual area (counts·s) | rows | median &#124;RPD&#124; | within ±25% |
|---|---|---|---|
| > 10⁶ | 2,391 | 0.7% | 88.1% |
| 10⁵ – 10⁶ | 6,316 | 1.0% | 91.2% |
| 10⁴ – 10⁵ | 13,017 | 1.7% | 89.4% |
| 10³ – 10⁴ | 19,516 | 3.3% | 84.0% |
| 10² – 10³ | 25,752 | 12.4% | 66.4% |
| < 10² | 51,333 | 109% | 17.7% |

The last two bands sit below the reporting threshold and are shown for completeness. By
homolog, area agreement runs from 1.1% (mono) to 5.2% (nona). Penta is the weakest of the
well-populated homologs at 3.0%, and carries the most crowded transition in the method.

Performance: cold start runs at 0.7-0.9 s median against a 3 s budget. A 29-file run processes
in seconds, and the interface does not stall for more than 100 ms while processing.

---

## 7. Known limitations

### 7.1 Congeners buried on a larger peak's leading edge

**1,673 detectable rows (4.1% of rows, but 12.1% of detectable area) are reported at 1% or
less of the analyst's value.** They are overwhelmingly flagged `absent` (1,576 of them), so
they are reported as flagged rather than passing silently.

Five congeners supply 675 of those rows: **83, 117, 115, 88 and 143**. Each elutes on the
*leading* edge of a much larger neighbour. The trace climbs monotonically through them, so
they produce no local maximum and no valley. The presence test rejects them. The shape fit of
section 2 either is not attempted or cannot justify its result. Where the fit does succeed,
it works: 1,701 rows are resolved that way, at 24.5% median |RPD| with 51% inside the band.
However, on these five the separation is genuinely ill-posed: two overlapping peaks can trade
amplitude and reproduce the observed envelope equally well.

This is the largest single source of remaining disagreement.

### 7.2 Mass agreement trails area agreement

78.7% of masses fall within ±25% against 87.1% of areas. A mass is a response factor times
an area times an internal-standard ratio, so it inherits the area disagreement and adds the
uncertainty of two further measurements. The gap is expected in kind. Its size is worth
attention.

### 7.3 A probable sample swap in one run

In one run a sample and a lab blank appear to have been exchanged: the file acquired as the
blank matches the analyst workbook's sample column, and vice versa. Congener 8 reads 527,071 in
the workbook's sample column against 533,898 in the file assigned as the blank, and 1,131 in
its blank column against 1,273 in the file assigned as the sample. Both `.D` files' internal
metadata and the sequence log are self-consistent, so the crossing happened at the bench or
during transcription.

Quant 209 reports the files as labelled and does not correct this. It inflates that run's
blank-derived limit of quantification roughly fivefold, since a real sample is being treated
as a blank. An analyst who agrees with the diagnosis can re-type the file in the Assign step.
Resolving it is a laboratory decision, not a software one.

### 7.4 A response-factor anomaly in one calibration injection

In one run, calibration injection `Cal03` yields 19 response factors more than 25% from the
analyst's, against 4 and 4 for the other two injections of that run. The cause is not
established. Response factors are averaged across a run's usable calibration injections, so
the effect on reported masses is diluted but not zero.

### 7.5 Near the detection limit

Below roughly 1,000 counts·s agreement degrades quickly: 12.4% median in the 100-1,000 band,
109% below 100. Manual and automated integration are both imprecise there. The laboratory
censors those results against the limit of quantification, so the disagreement does not reach
a reported number.

### 7.6 Scope of the evidence

* **Fifteen of the archive's 64 folders are not scored.** Five hold no analyst
  quantification at all, being either an instrument sequence log or a workbook superseded by
  a re-run. One holds no calibration injection. The other nine process cleanly and do have
  analyst workbooks, but their sample columns do not match the sample identifiers derived
  from the acquisition filenames, so no row joins. All fifteen are reported and skipped. The
  agreement figures cover the remaining 49.
* **The analyst is the reference, not the truth.** Response factors agree at about 1.9%
  against a ground truth whose own cross-workbook reproducibility is about 1.1%. Movement much
  below that would be fitting one analyst's hand rather than the chromatography.
* **Integration parameters were calibrated against two runs** and then measured against all
  49. The parameters are not re-fitted per run. Runs never used for calibration agree as well
  as the two that were.

---

## 8. Reproducing everything in this document

To reproduce every figure and number in this document, run:

```
.venv\Scripts\python.exe tools\master_bench.py  --out docs\validation\master
.venv\Scripts\python.exe tools\master_report.py --data docs\validation\master
.venv\Scripts\python.exe -m quant209.cli validate <run folder>
.venv\Scripts\python.exe tools\measure_startup.py --calibrate
.venv\Scripts\python.exe tools\measure_processing.py
```

`master_bench.py` runs the pipeline over every run in a directory and writes the pooled
comparison. `master_report.py` turns that into the report in `docs/validation/`. Both import
their joining and metric logic from `quant209/report/validation.py` (the same module the
application and the CLI build their validation report from), rather than reimplementing it. A
number in this document and a number on screen therefore come from one implementation.

`quant209 validate` writes a per-run directory. It holds `summary.md`, four agreement CSVs
(`area_`, `mass_`, `recovery_`, `rrf_agreement.csv`), `bias_diagnosis.csv` and three vector
figures. `bias_diagnosis.csv` lists every detectable native mass row outside the ±25% band,
with its integration flags, boundary method and a best-guess cause.
