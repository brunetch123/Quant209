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

Every calibration injection is anchored separately, and each file is integrated at the
anchors of the injection nearest it in time (`nearest_calibration`). Anchors are never averaged
across injections. An anchor is a position, and averaging positions across a sequence would
blur the drift the anchors exist to track. A calibration injection whose retention-time table
does not fit anchors nothing (section 3.2), and a file whose own labelled standards show that
part of its chromatogram moved is searched where they went (section 3.3).

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

An analyte that fails carries the **`absent`** flag, unless it is taken as a rider on a
neighbour's flank (section 2.6) or a named rule integrates it (section 2.7). Its area is
still reported, as a bounded estimate of whatever occupies its territory. This is what the analyst does with a shoulder they
chose to split off a dominant neighbour. Suppressing it would silently change
what the report means. However, `absent` is the flag every downstream consumer gates on.

### 2.3 Boundaries

Four methods, recorded per peak in `PeakResult.method`. The first two bound a peak that has
an apex of its own. The last two are for a rider, an analyte on a larger neighbour's flank
(section 2.6).

* **`valley_drop`**. At least one neighbouring analyte was itself a resolved peak, so
  the shared boundary is the valley between them (a perpendicular drop). Adjacent peaks'
  boundaries tile exactly: each peak's end is the next one's start.
* **`rt_window`**. No resolved neighbour bounded this peak. An isolated peak is bounded by
  the local minima within `max_halfwidth` = 0.25 min of its apex. A peak beside an *absent*
  neighbour extends through it to its own baseline edge. An absent analyte is bounded
  valley-to-valley inside its own midpoint territory.
* **`tangent_skim`**. A rider, bounded at its own foot and by a perpendicular drop on its
  parent's side, and integrated above a flat baseline at the level of its foot.
* **`deconvolved`**. A rider whose area is its own component's in a shape fit of the
  envelope it shares with its neighbours. The persisted boundaries say where it was looked
  for. The area is the fitted component's, not the trace between them.

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

Nine flags. Each says something different, and a peak carries its flags in the order of this
table, so the same set is always spelled the same way:

| flag | meaning |
|---|---|
| `absent` | the presence test rejected this analyte. The area is a bounded territory estimate, not a peak. **This is the detection verdict.** |
| `low_sn` | a detected peak whose baseline-corrected apex height is below `min_sn` = 3.0 × trace noise: real, but near the detection limit, so the area is imprecise. |
| `narrow` | fewer than 3 samples between the boundaries, so the area is dominated by where the boundaries landed rather than by the peak shape. |
| `clustered` | this area was measured above a cluster baseline rather than its own chord. **Not a confidence verdict.** It tells anyone recomputing the area by hand which baseline to use. |
| `shoulder` | the presence test rejected this analyte, but its neighbour's flank bulges at its retention time, so it was integrated as a rider (section 2.6). A rider's area depends on where its parent's flank was judged to take over, a softer judgement than a valley, so a reviewer should check it. |
| `forced` | neither detection test found this analyte in this file, and a named rule integrated it anyway (section 2.7). **Every other peak is evidence read off the trace. This one is a rule applied to a label.** |
| `displaced` | this analyte's peak drifted outside its own retention-time territory into a minor neighbour's, and was given back to it on the strength of published Aroclor profiles (section 2.8). The peak is measured from the trace. Which congener it belongs to is the judgement to check. |
| `apportioned` | this analyte and a neighbour of comparable abundance had merged into one envelope, and its area is its share by a two-peak fit (section 2.8). Both members of the pair carry the flag. |
| `unresolved` | this run's calibration did not separate this analyte from its neighbour, and in this file the two shared one envelope; its area is its component of a two-peak fit at the spacing the calibration measured (section 2.9). Both members of the pair carry the flag. |
| `realigned` | this analyte was searched at a retention time moved by a shift measured on this file's own labelled standards (section 3.3). The peak is read off the trace. Only where to look was corrected. |

### 2.6 Riders: the shoulder test and shape fitting

Implemented in `quant209/core/integrate.py` and `quant209/core/deconvolve.py`.

A congener eluting on the flank of a much larger neighbour has no local maximum of its own.
The trace is still climbing, or falling, through it, so its prominence is about zero and the
presence test rejects it. The analyst integrates it anyway, because the flank visibly bulges.

**The shoulder test.** The bulge is downward curvature: the second derivative of the smoothed
trace dips at the rider's apex and returns through zero at the inflections either side. An
analyte the presence test rejected is taken as a rider when the curvature at its strongest
point exceeds `shoulder_curvature_sn` = 6.0 × the trace noise and a neighbour is a resolved
peak. Without the second condition a ripple in an empty window would be promoted to a peak.
A rider carries the `shoulder` flag.

**Tangent skim.** A rider is bounded at its foot on the side away from its parent, and by a
perpendicular drop half way from its apex to the midpoint it shares with the parent. It is
integrated above a flat line at the level of its foot. A chord between its two boundaries
would run up the parent's convex flank and lie above the trace almost everywhere, and the
area would come out as zero.

**Shape fitting.** A skim approximates a rider. Fitting shapes separates it. The envelope,
the rider and every resolved peak within 0.75 min of it (three times `max_halfwidth`), is
fitted as a sum of Gaussians on a flat offset, by Levenberg-Marquardt implemented in numpy.
Three constraints carry the result, and each is a statement about the chromatography:

* **One shared width**, between 0.7 and 1.4 times the width measured on the parent. Peaks a
  few hundredths of a minute apart came off the same column under the same conditions. A
  rider allowed its own width widens until it explains its parent's whole rising edge.
  Unconstrained, the fit returns up to seven times the analyst's area on this archive.
* **Centres held within 0.4 of a peak width of where each was found**: a resolved peak at
  its apex, a rider at its point of strongest curvature, a named rider (section 2.7) at its
  anchor. The fit is asked how much of the envelope belongs to each, not where the peaks are.
* **A quality gate.** The fit is used only when its relative residual is at most 0.25 and the
  rider holds at most half of the fitted envelope. Otherwise the skim stands. A deeply
  buried rider is genuinely ill-posed: two components can trade amplitude and fit the
  envelope equally well.

An accepted fit reports the rider's own component area, and the method `deconvolved`.

### 2.7 Named rules for four congeners

Implemented in `_force`, in `quant209/core/integrate.py`. The lists are
`IntegrationParams.forced_riders` and `IntegrationParams.forced_present`.

Every rule above reads the trace and decides. These two read the label. Each applies only to
an analyte that neither the presence test nor the shoulder test found in a file, and each
names its congeners explicitly. A peak either rule produced is reported as detected and
carries the **`forced`** flag, so a reviewer can find every number that a rule, rather than
the trace, put there. Like a rider, a forced peak never moves a neighbour's boundaries.

* **`forced_riders` = 83, 117.** Fitted as riders (section 2.6) without the curvature
  evidence the shoulder test asks for, and centred on their anchors. Both elute so close to
  a larger neighbour that most files show no bulge at all. The fit's quality gate still
  applies: a fit that does not reproduce the trace is refused, and the analyte stays absent.
* **`forced_present` = 90+101+113, 118.** Bounded as if the presence test had accepted them.
  Both are dominant congeners of the technical mixtures. In a concentrated sample the peak
  broadens and its top moves off its anchor, out of reach of the apex search. The point
  examined in its place lies on the peak's flank, which has no prominence, and the largest
  peak in the chromatogram is reported absent (section 7.2). The apex recorded is
  the top reached inside the rule's boundaries.

**What they change.** Measured against version 1.0, which had neither rule, on the rows the
rules produced:

| congener | files where the rule fired | detectable rows | median &#124;RPD&#124; | within ±25% |
|---|---|---|---|---|
| 83 | 273 | 179 | 200% → 27.9% | 0% → 43% |
| 117 | 241 | 150 | 200% → 23.7% | 0% → 54% |

Version 1.0 reported these rows at a few counts·s against analyst areas in the tens of
thousands. The rule brings them to the agreement shape fitting reaches wherever it can be
justified (section 7.1). That is well short of a resolved peak's.

The other two rules fire mostly in blanks and reference injections (147 and 187 of their 176
and 243 files), on small peaks that fall below the presence test's absolute floor. Across
all their files, 152 and 223 rows below the reporting threshold agree with the analyst at a
median of 2.9% and 6.2%. Few of their rows are detectable, but they include the largest
peaks in the archive. In the five most concentrated samples the rule brings 90+101+113 to
within 0.3% of the analyst's 65 to 139 million counts·s. Version 1.0 reported every one of
them as zero. For 118, three concentrated samples come to 1.00, 1.38 and 2.00 times the
analyst's area.

**Where they do not help.** Seven detectable 90+101+113 rows and one 118 row are still
reported at under 2% of the analyst's value. In the clearest of these, the displaced apex
lies inside a neighbour's territory: the neighbour is credited with the peak, and the region
the rule integrates holds only a sliver of it. All of them carry `forced`, so they are in the
QA/QC table, and some also carry `low_sn`.

The rules add peaks. They change no other value.

### 2.8 Dominant peaks displaced across a boundary

Implemented in `_reclaim_displaced`, in `quant209/core/integrate.py`. The reference table is
`quant209/config/aroclor_abundance.csv` and the threshold is `IntegrationParams.dominance_ratio`.

In a concentrated sample a dominant peak broadens, and its top can drift up to about 0.2 min
off its anchor (section 7.2). When the top crosses the boundary an analyte shares with its
neighbour, the analyte's apex search, confined to its own territory, stops on the peak's
rising flank. That point has no prominence, so the presence test rejects the largest peak in
the chromatogram. The neighbour's search finds the top and claims the whole peak. The
neighbour is then reported at many times its real amount, and the dominant congener as
absent.

The trace alone cannot say whose peak it is. An analyst decides by knowing which congener is
the big one in these mixtures. Quant 209 takes that knowledge from published Aroclor
profiles (Rushneck et al. 2004): each analyte's largest weight fraction in Aroclors 1242,
1248, 1254 and 1260. Analyte *i* takes the top its neighbour *j* claimed only when all of
these hold:

* *i* failed the presence test, its apex sits on its own territory edge, and the trace is
  still climbing across that edge into *j*'s territory;
* climbing on from there, within `max_reach` of *i*'s anchor, ends exactly on *j*'s apex;
* *i* is a major congener (at least 1% of one of those Aroclors), and at least
  `dominance_ratio` = 10 times more abundant than *j*.

*i* then becomes present at the top, and *j* absent. Every boundary in the transition is
resolved with that verdict, so *j*'s window stops where *i*'s starts and no area is counted
twice. *j* is not offered to the shoulder test: on a large peak's flank, the flank's own
curvature clears any noise-scaled threshold, so a bulge there is no evidence of a rider. A
named rider (section 2.7) is still fitted. The reclaimed peak carries the **`displaced`** flag.

Two cases are handled explicitly:

* **Comparable congeners.** Where neither is ten times the other, 87+125 beside
  86+97+109+119 for example, the profiles cannot say whose peak it is, and often it is both:
  in a concentrated sample the two merge into one envelope. Where an absent major congener
  sits beside a present one of comparable abundance, the envelope is fitted as two peaks
  with one shared width, the spacing the calibration injection measured between them, and
  one common drift (`fit_pair`, in `quant209/core/deconvolve.py`). The fit stands only if it
  reproduces the envelope (relative residual at most 0.25), its drift stays inside
  `max_reach`, and each congener holds at least 5% of the envelope. Both congeners then take
  their fitted areas and carry the **`apportioned`** flag.
* **Congeners the profiles report as zero.** Eighteen congeners are zero in every Aroclor of
  the reference. In this archive's samples sixteen of them are genuinely trace, at a median
  of 0.014% of a sample's mass or less, and are taken as such. Two are not. The reference
  reports 91 together with 88, so the packaged table gives both their combined fraction and
  neither can take the other's peak. 102 makes up a median 0.19% of a sample's mass, so the
  reference evidently did not measure it apart, and its abundance is marked unknown. A
  congener of unknown abundance never takes a peak. It loses one to a dominant neighbour
  only if the top lies more than half the apex search window (0.05 min) from its own
  anchor, as a drifted peak does and its own peak does not.

**What it changes.** On the archive the rule gives 108 peaks back, in 80 files from 23
runs, all of them samples or blanks. Every reclaimed pair improves or is neutral, and
together they remove 1,084k ng of absolute mass error, more than a third of what version 1.1
left:

| reclaimed ← neighbour | peaks | mass error change |
|---|---|---|
| 110 ← 115 | 28 | −435k ng |
| 99 ← 112 | 11 | −186k ng |
| 129+138+163 ← 160 | 3 | −144k ng |
| 147+149 ← 134 | 2 | −122k ng |
| 95 ← 93+100 | 5 | −114k ng |
| 52 ← 73 | 5 | −71k ng |
| eight more pairs | 54 | −13k ng |

The pair fit splits 613 merged envelopes, in 308 files from 46 runs. Most involve small
amounts and change little. On the pair it was built for, 86+97+109+119 beside 87+125, it
removes 55k ng of error and reproduces the analysts' split within about 10% on the files
examined. It costs 5k ng on 92 and 90+101+113.

In the most concentrated samples the dominant congeners come back to within about 1% of the
analyst's area: 129+138+163 at 143.6 million counts·s against 143.8 million, and 147+149 at
140.7 million against 141.4 million.

**Where it does not help.** The rule gives a whole envelope to one congener. In
concentrated samples the analyst's integration sometimes divides the 110 envelope between 110
and 115 in proportions that differ from the Aroclor profiles: in 19 files across 7 runs it
gives 115 more than 110, where the profiles put 110 at about 80 times 115. Quant 209 keeps the
profiles' assignment, so it reports 115 as absent there. The 110/115 pair still carries 32% of the remaining
absolute mass error, most of it from those files. Section 7.2 covers the concentrated-sample
errors the rule does not reach.

### 2.9 Pairs the run's calibration does not resolve

Implemented in `unresolved_pairs` and `integrate_transition` (`quant209/core/integrate.py`),
`fit_unresolved_pair` (`quant209/core/deconvolve.py`) and `run_analysis`
(`quant209/core/pipeline.py`).

**The problem.** A few neighbouring congeners sit only 0.02 to 0.06 min apart in the canonical
table: 110/115, 45/51, 102/98, 134/143, 83/99 and 88/91. Most runs resolve them. Some columns
do not, and there the calibration injection shows one envelope where the standard holds two
congeners at comparable amounts. The congener without an apex of its own failed the presence
test, its whole envelope went to its neighbour, and it got no response factor: no mass in
any file of the run. In the validation dataset this took PCB 45 out of 30 runs, 83 out of 8,
51 out of 7 and 91 out of 5. It also inflated the neighbour's response factor, which was
built on both congeners' area.

**The rule.** The calibration standard holds every congener, so in its injection an analyte
with no apex, whose anchor lies inside a neighbour's integrated peak and within the apex
search window of its apex, is not absent: it eluted with that neighbour. Each usable
calibration is integrated once to find such pairs. The absent congener's anchor was only a
prediction, so the envelope is fitted as two Gaussian components with free positions and one
shared width, over the envelope alone and above its own baseline chord. Where the fit
reproduces the envelope (`MAX_RESIDUAL`) and each congener holds at least `MIN_PAIR_SHARE`,
the two fitted retention times become both congeners' anchors for that calibration.

In every file of the run, the calibration included, an unresolved pair that shows one envelope
is then split by the rigid pair fit of section 2.8 at that measured spacing, with one common
drift, over the envelope alone, and both members carry the **`unresolved`** flag. A congener
alone in a sample fits as one component and nothing, which the fit refuses, so a single peak
is never split. As for a comparable pair, the absent congener is not offered to the shoulder
test: a skim models a rider on a parent, not two congeners the calibration showed eluting
together, and on the case it was built on the skim gave the rider barely a quarter of its
area. A pair the calibration resolves, and a run without such pairs, integrate exactly as
before.

**What it changes.** On the four failing calibrations examined, the fit gives PCB 45 0.43,
0.46, 0.44 and 0.45 of the 45+51 envelope, against the analysts' 0.41, 0.46, 0.48 and 0.42, at
a fitted spacing of 0.073 to 0.085 min (the runs that resolve the pair show 0.074). On the
first of them, PCB 45's response factor comes out at 3.32 against the analyst's 3.27, where it
had none, and PCB 51's at 2.17 against 2.41, where it had 1.64.

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
  peak at all. A borrowed line can miss under a different GC program (section 7.4).
* The slope is bounded to 0.5 to 2.0. A method change re-times a chromatogram. It does not
  compress or stretch it by more than a factor of two. The extremes measured on this dataset
  are 0.791 and 1.494. Bounding the slope stops the search from "explaining" a channel with a
  degenerate near-flat map that piles every congener onto one peak.
* A fitted map displacing by less than `IDENTITY_TOL` = 0.05 min is reported as the identity
  and not applied. Fitting a line to apexes always returns something slightly off identity.
  On the run the canonical table was extracted from, the fitted maps displace by at most
  0.008 min, which is apex noise, not a method difference. Declaring identity there keeps a
  same-method run's anchors exactly what they would be without alignment.

**A near tie is settled by the standard's composition.** In a crowded channel, a line shifted by
one peak across part of the channel puts each congener beside its neighbour's peak and matches
about as many peaks as the right one. On one re-run, two of three calibrations chose such a
line in the hexa channel by a single match (29 against 28), and the run's cross-calibration
reconciliation then "repaired" the third, correct one to agree with them. Every hexa congener
from 164 to 158 took its neighbour's peak, and PCB 141's response factor came out 34 times too
large. An analyst rejects that at a glance: the 129+138+163 group holds three congeners and
must be the big peak. So among lines within `SIZE_COUNT_SLACK` = 2 matches of the best, one whose
matched peaks agree with the standard's composition (correlation of log certified mass with
log peak height) by at least `SIZE_CORR_MARGIN` = 0.3 more is taken instead. There the right
line scored +0.36 and the shifted one -0.11. The tie must be one between near neighbours: the
alternative may move no position by more than `CONSENSUS_TOL` = 1.0 min from the best line's,
and the composition is read only where at least `SIZE_MIN_MATCHED` = 12 congeners match. On
the eleven-congener di and octa channels a correlation over a handful of masses is noise, and
without these bounds it once preferred a line that re-timed the octa channel by 8.5 min and
lost the internal standard on it.

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

**Every calibration injection is asked the question** (`select_calibrations`). One whose
table does not fit, while another does, is a different chromatographic frame or a failed
injection, and its anchors would send every file integrated at them to the wrong retention
times. It therefore anchors no other file and contributes no response factor. It is named in
a run note and in the Review screen's attention table, and the run's verdict and alignment come
from the first injection that fits. Measured case: one run's first calibration sat 4.3 min late
against its other three. The run failed this guard, and the four files nearest it (a sample, a
reference injection and two blanks) were integrated at its retention times and lost both
internal standards. With it excluded, the run passes, the sample is quantified, and that run's
recoveries and its two late surrogates' masses come right. If no injection fits, nothing is
excluded and the run fails as before.

### 3.3 A file whose labelled standards moved

Each file is integrated at the retention times of the calibration injection nearest it. That
is right as long as the file's chromatography matches its calibration's. In one archive file
it did not: everything from 16 to 24 min eluted 0.2 min late. The internal standard 30D and
three surrogates were searched in the wrong place and lost, and with the standard went the
masses of about a hundred congeners. The sample came out at 170 ng against the analyst's 1,486.

The labelled standards are spiked into every file, so they show where its chromatography went.
`realign_shifted_file`, in `quant209/core/pipeline.py`, acts only when all of these hold:

* at least two labelled standards on channels that carry no native congener
  (`landmark_labels`) were lost in this file, while the run's other files find them;
* each is found within 0.4 min of where it was expected, at 15% or more of its typical height
  in the run's samples and blanks;
* they agree on the shift to within 0.1 min, and the shift is larger than the apex search
  window.

Every analyte from the first to the last moved standard, plus 1 min either side, is then
searched that much later. The new integration is kept only if it recovers those standards,
loses no standard that was found before, and finds no fewer native peaks in that span.
Otherwise the file is left exactly as integrated. Moved peaks carry the **`realigned`** flag,
and the run note says which standards moved, by how much, and how many native peaks were found
before and after.

On the archive it acts on 2 of 672 reported files and changes no other. The file above now
agrees with the analyst, at 1,511 ng against 1,486. The other is an injection the analyst left
unquantified: two of its surrogates are recovered, and its earliest region, which moved beyond
the search, keeps the missing-internal-standard error of section 4.2.

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

**No internal standard, no mass.** The internal standard is the one exception. A standard the
presence test rejects is never used, so nothing it normalises gets a mass: a number scaled by
a standard that was never measured would be a guess. The file then carries an error, as
`missing_internal_standard:<label>=<reason>`, and a run note. The reason says which remedy
applies:

* `absent`: the standard is not in the injection. Its channel's tallest point within 1 min of
  the expected time stays under 5% of the standard's typical height in the run. Measured case:
  one archive injection tops out at 152 and 223 counts on the two standard channels, against
  standard peaks of 9,000 to 14,000 in the files around it, while all its surrogates are
  present.
* `not_found`: the channel carries a peak nearby, but no standard was detected where it was
  expected. That is a retention-time problem (sections 3.2 and 3.3).
* `rejected`: detected, and refused by the consistency guard of section 7.3.

The Review screen shows it as an error: a red banner, an "internal standard missing" row in the
attention table, and a "missing IS" event in the QA/QC table.

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
spans 49 runs and 41,376 detectable peak comparisons, across four years, several analysts and
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

**Every comparison is checked to be the same data file** (`verify_pairing`, in
`quant209/report/validation.py`). An analyst sheet and a data file are matched by name, and a
name can be wrong. The evidence is the labelled standards: for the same file our areas and the
analyst's agree to a median 0.4% (99% within 1.2%), and different files differ by several
percent. The native congeners being scored are never consulted, so the check cannot pick a
pairing for its agreement. A sheet is re-paired only when it and another sheet are each a
decisive match for the other's file, a mutual crossing. Anything else that disagrees is
reported as unverified and left in place. The benchmark writes every verdict to
`pairing.csv`. Masses the program refuses because an internal standard is absent from the
injection (section 4.2) are not scored, since there is nothing of the program's to compare and
the analyst's figure rests on standard areas that were not measured. A standard the program
merely failed to find stays scored.

**The analyst's final values.** A run folder can hold several issues of a workbook covering
the same samples, and the last one read wins. Reading them in filename order made that the one
whose name sorts last, which need not be the latest: a tag run into a name sorts after a later
date stamp, because a letter sorts after a digit, and words such as `OLD`, `test` and `Copy`
sort late too. Workbooks are therefore read in the order they were issued: those marked `OLD`,
`test` or `Copy` first, then by the date stamp in the name (or the file's own date when there
is none), then by name.

**Group masses on one basis.** A co-eluting group's calibration peak holds every congener in
it, so its calibration mass is the sum of theirs, and that is the basis Quant 209 reports on. A
workbook can instead enter one member's mass for a group, which puts its group masses on a
per-congener basis, a half or a third of the summed figure. Where a workbook's own `cal_info`
mass for a group sits within 10% of 1/N of the certified sum, the masses of the samples read
from that workbook are put on the summed basis by the ratio (`rebase_group_masses`), so that
every comparison is like for like. The evidence is the workbook's own calibration masses; no
scored value is consulted.

### 6.2 Results

| | rows | median &#124;RPD&#124; | within ±25% | total ratio |
|---|---|---|---|---|
| Peak areas | 41,376 | **1.99%** | **89.8%** | 1.008 |
| Congener masses | 30,712 | **4.14%** | **85.5%** | 1.015 |
| Surrogate recoveries | 6,711 | 0.021 mean absolute error | | |

Areas and masses both clear both gates. Mass agreement still trails area agreement, for the
reason given in section 7.5.

**Since version 1.0** (2.13% and 87.1% for areas, 4.61% and 78.7% for masses), every change
comes from the named rules of section 2.7, the displaced-peak rule of section 2.8, reading
diluted injections as dilutions (section 7.3), excluding a calibration injection that does not
fit (section 3.2), searching a shifted file where its standards went (section 3.3), scoring
every analyst sheet against the file it came from, and putting the analyst's group masses on
one basis (both section 6.1), and, in version 1.3, splitting the pairs a run's calibration does
not resolve (section 2.9) and settling near ties in the retention-time alignment by the
standard's composition (section 3.1). The area median barely moves, because these act on few
rows. What moves is the large error: total absolute mass error on detectable rows fell from
3,214k ng in version 1.0 to 3,004k with the named rules, to 1,565k in version 1.2 and to 1,474k
now. No sample's native total differs from the analyst's by 2× or more, except two trace-level
samples where the difference is a few nanograms, and only 11 of 480 by 25% or more, seven of
them under 11 ng. The total ratios went from 0.975 to 1.008 for areas and from 0.991 to 1.015
for masses.

No run in the archive disagrees by more than 7.4% at the median. The best sits at 0.67%.
Agreement tracks signal:

| manual area (counts·s) | rows | median &#124;RPD&#124; | within ±25% |
|---|---|---|---|
| > 10⁶ | 2,394 | 0.6% | 93.9% |
| 10⁵ – 10⁶ | 6,321 | 0.9% | 95.0% |
| 10⁴ – 10⁵ | 13,076 | 1.6% | 92.3% |
| 10³ – 10⁴ | 19,585 | 3.1% | 86.1% |
| 10² – 10³ | 25,903 | 11.8% | 68.1% |
| < 10² | 52,084 | 108% | 18.0% |

The last two bands sit below the reporting threshold and are shown for completeness. By
homolog, area agreement runs from 1.1% (mono) to 5.0% (nona). Penta is the weakest of the
well-populated homologs at 2.7%, and carries the most crowded transition in the method.

Performance: cold start runs at 0.7-0.9 s median against a 3 s budget. A 29-file run processes
in seconds, and the interface does not stall for more than 100 ms while processing.

---

## 7. Known limitations

### 7.1 Congeners buried on a larger peak's leading edge

**1,174 detectable rows (2.8% of rows) are reported at 1% or less of the analyst's value.**
They are overwhelmingly flagged `absent` (1,082 of them), so they are reported as flagged
rather than passing silently. Version 1.0 had 1,673. Most of the difference is congeners 83
and 117, now handled by a named rule (section 2.7). The rows are two different failures. The
1,099 below 10⁶ counts·s are many but small: they carry 0.4% of the analyst's detectable
area. The 75 above it carry 3.7%, and are the subject of section 7.2.

Four congeners supply 417 of the small rows: **115, 88, 143 and 98**. Each elutes on the
*leading* edge of a much larger neighbour. The trace climbs monotonically through them, so
they produce no local maximum and no valley. The presence test rejects them, and either the
shoulder test finds too little curvature or the shape fit of section 2.6 cannot justify its
result. Where the fit does succeed, it works about as well as a fit can here: 3,007
detectable rows are resolved that way, at 16.3% median |RPD| with 63% inside the band. On
these four, however, the separation is genuinely ill-posed: two overlapping peaks can trade
amplitude and reproduce the observed envelope equally well.

### 7.2 Concentrated samples: dominant peaks the rule does not reach

**75 detectable rows, in 41 files from 13 runs, are dominant peaks the analyst integrates at
more than 10⁶ counts·s and Quant 209 reports at 1% or less.** They carry 3.7% of the
analyst's detectable area, down from 8.7% before the displaced-peak rule of section 2.8, and
they remain the largest disagreement by area. The congeners are major components of the
technical mixtures, led by 86+97+109+119 (17 rows) and 129+138+163 (15).

The mechanism was traced on the largest cases. In a concentrated sample a dominant peak
broadens, to half-height widths of 0.12 to 0.17 min on the files examined, and its top moves
off its anchor by up to 0.18 min. The apex search is confined to the analyte's own territory,
so it stops on the peak's flank, where there is no prominence, and the analyte is reported
absent. The files' total areas still agree with the analyst's, at a median of 0.99, so the
area is measured and misattributed rather than lost. Section 2.8 resolves the clear cases.
Of the 75 that remain:

* 26 have no single neighbour holding the missing area. It is spread over several analytes.
* 24 are a major congener beside a minor one, but the drift is not confirmed: the apex is
  not on the territory edge, or the climb does not end on the neighbour's apex within
  `max_reach`.
* 15 are the other half of an envelope the analyst splits, as 115 beside 110 (section 2.8).
* 10 involve congeners the profiles cannot rank: a neighbour of comparable abundance, or a
  congener that is not major.

The dominant congener carries `absent` and appears in the QA/QC table. A neighbour that took
its area carries no flag. The general remedy remains a per-file retention-time correction
for concentrated samples, which Quant 209 does not yet make.

### 7.3 Samples the pipeline declines to quantify

**One of the archive's 480 samples is refused, by design.** Its injection contains no internal
standard (section 4.2): the standard channels are flat, while every surrogate is present and
agrees with the analyst's. Quant 209 reports an error and no mass, and the benchmark does not score that sample's
masses (section 6.1). Its peak areas agree with the analyst's and stay scored.

Two other samples were refused before, for reasons that were the program's own: one was
anchored to a calibration injection that did not fit (section 3.2), and one had part of its
chromatogram shifted (section 3.3). Both are now quantified and agree with the analyst's
totals.

The consistency guard (`reject_inconsistent_standards`) drops a standard outside 0.2 to 5
times its run's median (`IS_CONSISTENCY_RANGE`), to catch a mis-integrated standard, which
would scale every congener that shares it. A mis-integrated standard moves on its own. A diluted
or concentrated injection moves every standard together, and leaves every analyte-to-standard
ratio, and so every mass, unchanged. A file whose standards all leave the band together, within
a factor of 1.5 of each other (`IS_DILUTION_AGREEMENT`) and no further than a hundredfold
(`IS_DILUTION_RANGE`), is therefore read as a dilution: its standards are kept and the reading
is noted. Six archive samples are dilutions of this kind, their two standards at 0.11 to 0.14
of the run median and agreeing within 6%. Version 1.1 refused them, which reported 137k ng of
analyst mass as unquantified.

### 7.4 Late labelled standards under a different GC program

A labelled standard is alone on its channel, so its retention time is predicted from a
neighbouring channel's map (section 3.1). On three archive runs acquired under a shorter GC
program, that borrowed prediction puts 13C 194 1.4 min from its own peak, and its area comes
out at a few percent of the analyst's. The recoveries computed from it there still agree with
the analyst's, because a recovery is a ratio against the reference injection, measured at the
same position. Masses use the internal standards and are unaffected.

Two fixes were tried and not kept, because each repaired some files and broke others:
anchoring a labelled standard on its channel's tallest calibration peak, which moved the nona
and deca standards onto interfering peaks in one run, and anchoring it beside its native twin,
which on these runs moved it onto no peak at all.

Octa recovery agreement, which this section previously reported as the weakest homolog at
0.095 mean absolute error from a single run of unknown cause, turned out to be that run's
off-frame calibration injection (section 3.2): its reference injection was anchored there. With
that injection excluded, octa recoveries agree to 0.023, in line with the other homologs.

### 7.5 Mass agreement trails area agreement

85.5% of masses fall within ±25% against 89.8% of areas. A mass is a response factor times
an area times an internal-standard ratio, so it inherits the area disagreement and adds the
uncertainty of two further measurements. The gap is expected in kind. Its size is worth
attention.

### 7.6 Pairing analyst sheets with data files

Quant 209 reports every file under the name and type it was recorded with, and re-labels
nothing: which label is right is a laboratory question. The benchmark checks every analyst
sheet against the file it is scored against, and scores a sheet against another file only on a
decisive, mutual match of the labelled standards (section 6.1). A file recorded as a blank that
holds a sample's signal raises its run's blank-derived limit of quantification. Files can be
re-labelled or re-typed in the Assign step.

### 7.7 A response-factor anomaly in one calibration injection

In one run, calibration injection `Cal03` yields 19 response factors more than 25% from the
analyst's, against 4 and 4 for the other two injections of that run. The cause is not
established. Response factors are averaged across a run's usable calibration injections, so
the effect on reported masses is diluted but not zero.

### 7.8 Near the detection limit

Below roughly 1,000 counts·s agreement degrades quickly: 12.3% median in the 100-1,000 band,
109% below 100. Manual and automated integration are both imprecise there. The laboratory
censors those results against the limit of quantification, so the disagreement does not reach
a reported number.

### 7.9 Scope of the evidence

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
* **The named rules were chosen on the archive they are scored against.** The four
  congeners of section 2.7 were picked by reading this comparison, so their contribution to
  the figures above is in-sample.
* **The thresholds of the same-file check (section 6.1) and of the per-file retention
  correction (section 3.3) were set while examining this archive.** The correction acts on two
  files in it.
* **The two ground-truth corrections of section 6.1 change what the program is scored
  against, never what it reports.** Each rests on evidence the scored values do not enter: the
  workbooks' names and dates, and their own calibration masses.
* **The displaced-peak rule rests on published Aroclor profiles, not on the analysts'
  numbers, but two of its details were settled on this archive.** The half-window position
  check for congeners the profiles report as zero, and the correction that gives 91 the
  fraction the reference reports for 88 and 91 together, were both added after inspecting
  the rule's misfires here.

### 7.10 Five penta congeners in the 2025 runs

In the 2025 calibrations, five congeners of the penta transition are identified one peak too
early: 87+125, 117, 85+116, 110 and 115. The alignment places them correctly, but the
elution-order search that follows (`assign_in_elution_order`) moves them. The four-congener
group 86+97+109+119 shows there as two sub-peaks, and 115 merges into 110's tail, so the window
has one peak too many in one place and one too few in another. The search, which values
explaining every congener above placing each one close to its prediction, shifts the five by
one peak each, 117 by 0.77 min. All four calibration injections shift alike, so comparing them
does not catch it, and every sample inherits the anchors. 87+125's response factor comes out
thirteen times too high, and the five congeners' masses in those runs are not trustworthy. They
carry most of the remaining excess in the two runs that still read 1.14 to 1.18 times high.

A fix was built and measured. Limiting how far the search may move a congener the alignment
identified on its own peak repairs these runs completely, but on fourteen other runs the same
search makes correct moves of that size (190 by 0.8 min, 203 by up to 1.2), and the limit
blocked them: area agreement fell from 88.8% to 87.6% within the band. Neither the distance
moved nor the margin by which the search's choice wins separates the two cases, so the fix was
not kept.

---

## 8. Reproducing everything in this document

To reproduce every figure and number in this document, run:

```
.venv\Scripts\python.exe tools\master_bench.py  --out docs\validation\master
.venv\Scripts\python.exe tools\master_report.py --data docs\validation\master
.venv\Scripts\python.exe tools\validation_figures.py --data docs\validation\master
.venv\Scripts\python.exe -m quant209.cli validate <run folder>
.venv\Scripts\python.exe tools\measure_startup.py --calibrate
.venv\Scripts\python.exe tools\measure_processing.py
```

`master_bench.py` runs the pipeline over every run in a directory and writes the pooled
comparison. `master_report.py` turns that into a full report, `docs/validation/MASTER-REPORT.md`,
from which `docs/validation/README.md` is taken. Both import their joining and metric logic
from `quant209/report/validation.py` (the same module the application and the CLI build their
validation report from), rather than reimplementing it. A number in this document and a number
on screen therefore come from one implementation.

`validation_figures.py` draws the validation figures from the same pooled comparison, as SVG:
agreement by method and flag, the area and mass scatters, Aroclor fingerprint agreement, and
breakdowns of where and why the two integrations differ. The Aroclor figures need a
congener-level reference table, given with `--aroclor`. The pooled comparison, the full report
and the figures all name run folders or plot per-sample results. They are laboratory data and
are not committed.

`quant209 validate` writes a per-run directory. It holds `summary.md`, four agreement CSVs
(`area_`, `mass_`, `recovery_`, `rrf_agreement.csv`), `bias_diagnosis.csv` and three vector
figures. `bias_diagnosis.csv` lists every detectable native mass row outside the ±25% band,
with its integration flags, boundary method and a best-guess cause.
