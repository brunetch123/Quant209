# Changelog

## 1.4.1 (2026-11-23)

### Integration

- **A peak's window reaches 0.25 min in every part of a trace.** A side of a peak with no
  resolved neighbour may reach `max_halfwidth`, 0.25 min, from the apex. The reach was counted in
  samples of the trace's median sampling interval, and a transition's MRM time segments do not
  all sample it at that rate: where the hexa channel is sampled every 0.305 s instead of 0.549 s,
  the count reached 0.137 min. 156+157, PCB 169, 128+166 and PCB 126 were cut on their own tails
  there, and the straight baseline from the floor to the cut took a triangle out from under them.
  169, the hexa channel's last congener, was cut in every calibration injection, so its response
  factor ran high, and 156+157 in samples whenever 169 is absent. The reach is now measured on the
  trace's own retention times. See [`docs/METHODS.md`](docs/METHODS.md) section 2.3.
- **Effect.** On the benchmark, rows within ±25% of the analyst rise from 91.5% to 91.8% for peak
  areas and from 87.3% to 87.7% for masses. 169's response factor lies within ±25% of the
  analyst's in 67.2% of runs (27.6% before), 126's in 98.3% (96.6%), and no run's area agreement
  falls by a percentage point or more. On the validation dataset, 300 area rows move inside ±25%
  of the analyst and 160 outside (masses 258 in, 55 out). 156+157 in samples of runs under the
  shorter oven program goes from 54.6% to 92.7% within ±25%, and 169's response factor lies
  within ±25% in 30 instead of 12 of 55 such runs and in 69 instead of 50 of 71 others. Some
  low-level peaks whose windows now reach further take in what the analysts leave out, mostly
  167, 170, 162, 172 and 190.

### Retention-time alignment

- **A copy of another congener is never taken for a calibration peak.** In a calibration
  standard every congener leaves a copy of its peak on the channels one and two chlorines down,
  among the analytes' own peaks, and an analyte predicted in the wrong place could be anchored on
  one. On runs under a shorter oven program, PCB 169 was anchored on PCB 170's copy and PCB 189
  on PCB 208's, and every file of those runs was integrated there. A calibration peak within
  0.012 min of a peak on a parent homolog's channel, and no larger than twice the copy that
  homolog can leave, is no longer offered to the elution-order search. See
  [`docs/METHODS.md`](docs/METHODS.md) section 3.1.
- **A channel whose anchors disagree with the other channels is searched again.** On the same
  runs, the straight line fitted per channel placed the last congeners of the hexa to nona
  channels up to 1.7 min early, and the search gave each of them the peak before its own. Each
  anchor is now compared with the drift the run's other channels show at its retention time, and
  a channel with anchors more than 0.15 min off it is searched again, predicted from that drift.
  A channel whose anchors all agree is left exactly as it was. This also places 87+125, 117,
  85+116 and 110 on their own peaks in the 2025 runs of METHODS section 7.10.
- **Effect.** Against the retention times on the analysts' calibration sheets, anchors more
  than 0.05 min off fall from 4.49% to 1.60% on runs under the shorter oven program and from
  0.05% to 0.01% on the others. On the benchmark, rows within ±25% of the analyst rise from
  90.5% to 91.5% for peak areas and from 86.5% to 87.3% for masses, and the total mass
  disagreement falls by 9.5%. 189's response factor lies within ±25% of the analyst's in 95.7%
  of runs (48.9% before), 195's in 98.3% (65.5%). On the validation dataset, areas within ±25%
  rise from 87.5% to 88.3% and masses from 86.7% to 87.5%, and response factors not issued fall
  from 88 to 6. No run's area agreement falls by a percentage point or more. Agreement that
  rested on a misplaced anchor goes: 156+157's window no longer runs on to 169's, and 115, which
  shares 110's envelope, is reported absent in more files.

### Documentation

- [`docs/METHODS.md`](docs/METHODS.md) sections 2.6 and 2.8 said a rider is fitted with every
  resolved peak within 0.75 min of it and that the climb to a displaced peak's top stays within
  0.20 min of the anchor. Both are counted in samples at the trace's median sampling interval, so
  they span those distances only where a trace is sampled at its median rate. The sections now
  say so, and section 2.3 gives the measurement behind keeping them: measured in minutes, they
  moved no headline figure and gave 182's peak to 187 in 12 benchmark files.
- The README, METHODS sections 6 and 7 and
  [`docs/validation/README.md`](docs/validation/README.md) report this version's validation
  figures: peak areas 1.86% median |RPD| and 91.8% within ±25%, congener masses 3.94% and 87.7%,
  surrogate recoveries 0.021 mean absolute error, and no run worse than 3.6% at the median.
  METHODS section 7.9's response-factor agreement had kept version 1.0's figure and is brought up
  to date.

## 1.4.0 (2026-11-22)

### Integration

- **Ghost peaks are kept out of integrations.** A PCB with n chlorines that loses HCl or Cl2 in
  the ion source leaves a copy of its peak on the channel for n-1 or n-2 chlorines, at its own
  retention time and with its shape. The internal standard 204 puts 4-6% of its area on the
  hepta channel and 14-16% on the hexa channel in every injection that holds it, where the copy
  was reported as 156+157 or taken into the windows of 172 and 171+173. Native congeners do the
  same, and in the calibration standards their copies fell inside windows that response factors
  are measured in. Each copy is now predicted from a per-congener table
  (`quant209/config/ghost_fractions.csv`) and the run's own 204 fraction, read in its blanks,
  and confirmed by fitting the parent peak's shape. A confirmed copy between analytes is fenced
  off so no window takes it. One on an analyte's expected retention time is subtracted, never
  beyond the prediction, except under PCB 58 and PCB 159, whose reported peaks usually hold the
  copy of PCB 95 and PCB 174. A copy the trace cannot confirm is left alone. Changed peaks carry
  the new `ghost` flag. `IntegrationParams.ghost_correction = False` turns it off. See
  [`docs/METHODS.md`](docs/METHODS.md) section 2.10.
- **Effect.** On the benchmark, rows within ±25% of the analyst rise from 89.8% to 90.5% for
  peak areas and from 85.5% to 86.5% for masses, sample totals within ±25% from 97.7% to 98.1%,
  and rows reported where the analyst has no peak fall from 2,742 to 2,041. 105's response
  factor lies within ±25% of the analyst's in 98.3% of runs (63.8% before), 156+157's in 82.8%
  (58.6%). On the validation dataset, areas within ±25% rise from 86.9% to 87.5% and masses from
  85.7% to 86.7%. A run takes about 30% longer.

### Retention-time alignment

- **A file's drift is followed along a curve, not by one shift.** A file whose chromatography
  drifted from its calibration's was searched at one shift, which fits a step but not a slope.
  In two concentrated runs the labelled standards sat 0.05 min late at 23 min and 0.33 min late
  at 67 min, and a search 0.1 min off credits a congener with its neighbour's peak. The file's
  drift is now a curve through every labelled standard it shows, and each analyte is searched
  where the curve puts it. See [`docs/METHODS.md`](docs/METHODS.md) section 3.3.
- **A shift is kept only if the native congeners corroborate it.** Standards that sat off their
  calibration anchors in every file of a run set the `realigned` flag on files whose natives had
  not moved. The new integration must now find more native peaks where it moved the search than
  the unshifted one, so a file with nothing in the moved region to check it against is left
  alone, and a peak several times a standard's usual height is no longer taken for the
  standard.
- **A lost internal standard starts the same search when the other standards show a drift.**
  A file whose labelled standards were all still found, but 0.1 to 0.2 min late, lost its
  internal standard 204 on a crowded channel and with it every mass that standard normalises.
  Three standards found where expected, with a median offset beyond half the apex search
  window, now start the search for it.
- **Effect.** On the 49 benchmark runs it acts on 1 file instead of 2 and changes no other:
  that file's rows within ±25% of the analyst rise from 84.5% to 90.5% for peak areas and from
  77.5% to 83.8% for masses. In archive runs added since the benchmark it acts on 18 more files.
  Across the 21 files whose results change, rows within ±25% of the analyst rise from 54.8% to
  62.0% for peak areas and from 57.6% to 67.9% for masses (the 14 changed files the analyst
  quantified), and six files that had lost their internal standard 204 get it back, with the
  masses it normalises.

### Documentation

- The README, [`docs/METHODS.md`](docs/METHODS.md) sections 6 and 7 and
  [`docs/validation/README.md`](docs/validation/README.md) report this version's validation
  figures: peak areas 1.94% median |RPD| and 90.5% within ±25%, congener masses 4.03% and 86.5%,
  surrogate recoveries 0.021 mean absolute error. METHODS sections 7.1, 7.2 and 7.8 and
  `docs/validation/README.md` had kept earlier versions' figures and are brought up to date.

## 1.3.0 (2026-11-20)

### Integration

- **Pairs the run's calibration does not resolve are split.** A few neighbouring congeners
  that most columns separate (45/51, 83/99 and 88/91 in the archive) elute as one envelope on
  some. The calibration injection then gave one of them the whole envelope and the other no
  response factor, so no mass in any file of the run. Each calibration now measures both
  congeners' retention times by a two-peak fit with free positions, and every file splits such
  a pair's shared envelope by a two-peak fit at that spacing. A single peak is never split.
  Both members carry the new `unresolved` flag. See [`docs/METHODS.md`](docs/METHODS.md)
  section 2.9. On the 58 benchmark runs, PCB 45 within ±25% of the analyst rises from 15% to
  50% of rows and PCB 51 from 13% to 28%.

### Retention-time alignment

- **A near tie is settled by the standard's composition.** In a crowded channel a line
  shifted by one peak can match as many peaks as the right one. Two of one re-run's three
  calibrations chose such a line, and reconciliation then moved the third to agree, so a whole
  window of congeners took their neighbours' peaks. Among lines within two matches of the best,
  one whose matched peak heights agree clearly better with the standard's certified masses is
  now taken, provided it moves no position by more than 1 min and at least 12 congeners match
  (section 3.1). With both changes, on the 58 benchmark runs, rows within ±25% of the analyst
  rise from 84.1% to 85.5% and absolute mass error falls by 6% (1,565k to 1,474k ng); three
  runs improve by 10k to 58k ng and none gets worse; 129+138+163 within ±25% rises from 81% to
  91%.

### Documentation

- The README and [`docs/METHODS.md`](docs/METHODS.md) section 6 report this version's
  validation figures: peak areas 1.99% median |RPD| and 89.8% within ±25%, congener masses
  4.14% and 85.5%, surrogate recoveries 0.021 mean absolute error.

## 1.2.0 (2026-11-16)

### Integration

- **Dominant peaks displaced across a boundary are given back.** In a concentrated sample a
  dominant peak broadens and its top can drift into a minor neighbour's retention-time
  window, where the neighbour was credited with all of it and the dominant congener reported
  as absent. When the drift is unmistakable, and published Aroclor profiles put the dominant
  congener at least ten times above its neighbour, the peak goes back to it and the two are
  re-split, so nothing is counted twice. Such peaks carry the new `displaced` flag. See
  [`docs/METHODS.md`](docs/METHODS.md) section 2.8.
- `IntegrationParams.dominance_ratio` (default 10) sets the abundance advantage the rule
  needs, and 0 turns it off. The reference fractions ship as
  `quant209/config/aroclor_abundance.csv`, from Rushneck et al. (2004). Congeners the
  reference reports as zero are taken as trace, except 91, which it reports with 88 (both
  carry the combined fraction), and 102, whose abundance is marked unknown.
- **Comparable congeners merged into one envelope are split.** Where neither congener is ten
  times the other, the envelope is fitted as two peaks with one width, the calibration's
  spacing and one common drift. Both take their fitted areas and carry the new `apportioned`
  flag. On 86+97+109+119 beside 87+125 the fit reproduces the analysts' split within about
  10%.

- **A file whose labelled standards moved is searched where they went.** When part of one
  file's chromatogram elutes off its calibration's retention times, the analytes there were
  searched in the wrong place and lost, internal standard included. If two or more labelled
  standards were lost and are found nearby, agreeing on the shift, that span is searched again
  at the shift, kept only if it recovers the standards and loses nothing. Moved peaks carry
  the new `realigned` flag.
- **A calibration injection that does not fit anchors nothing.** Every injection is now asked
  the run-fit question. One that fails, while another passes, anchors no other file and
  contributes no response factor. It is named in a note and in the Review screen.

### Quantification

- **No internal standard is an error, and no mass.** A file without a usable internal
  standard carries `missing_internal_standard:<label>=absent|not_found|rejected` and a run
  note, and the Review screen shows it in red. Nothing the standard normalises is given a mass.
- **A diluted injection is no longer refused.** The internal-standard check dropped any
  standard more than five times off its run's median, so a diluted sample, whose standards
  fall together, was reported unquantified. A file whose standards all move together, within
  a factor of 1.5 of each other, is now read as a dilution: every mass is an
  analyte-to-standard ratio, so the dilution cancels. Six archive samples are quantified that
  were refused before.

### Validation

- **Every comparison is checked to be the same data file**, by the labelled standards spiked
  into every file. A sheet is scored against a different file only on a decisive, mutual
  match. Every verdict is written to `pairing.csv`. Masses refused because an internal standard is absent
  from the injection are not scored.
- **The ground truth follows the analyst's final workbook.** Re-issues are read in the order
  they were issued (marked `_OLD`/`test`/`Copy` first, then by date stamp), not by how their
  names sort.
- **Group masses are compared on one basis.** Where a workbook enters one member's
  calibration mass for a co-eluting group, its group masses are put on the summed basis Quant
  209 reports on, so every comparison is like for like.
- The benchmark also writes `rrfs.csv`: every run's response factors, ours and the analyst's.
- Figure 41 (total mass per sample) sums native congeners only. It had included the spiked
  labelled standards.
- On the archive, total absolute mass error on detectable rows falls by 48%, from 3,004k to
  1,565k ng. Areas agree at 2.04% median and 88.8% within ±25%, masses at 4.24% and 84.1%,
  recoveries to 0.021. No sample's native total differs from the analyst's by 2× or more,
  apart from two trace-level samples, and 11 of 477 by 25% or more.

## 1.1.0

Released privately only. Its changes are part of 1.2.0.

### Integration

- **Named rules for four congeners.** Where neither detection test finds them, congeners 83
  and 117 are fitted as riders on their larger neighbour, and 90+101+113 and 118 are bounded
  as if present. The lists are `IntegrationParams.forced_riders` and
  `IntegrationParams.forced_present`. Setting either to `()` turns its rule off. See
  [`docs/METHODS.md`](docs/METHODS.md) section 2.7 for what the rules recover and where they
  do not help.
- **A new flag, `forced`,** marks every peak a named rule produced. It appears in the QA/QC
  table and in every export, and the workbook's run-info sheet lists the named congeners. A
  rule may add a peak. It never moves another peak's boundaries.

### Validation

- Archive agreement on 41,240 detectable peaks: areas 2.12% median, 87.4% within ±25%
  (1.0.0: 2.13%, 87.1%). Masses 4.60% and 79.0% (4.61%, 78.7%). Every change comes from the
  named rules, and no other row's value moved.
- [`docs/METHODS.md`](docs/METHODS.md) now documents the shoulder test, tangent skim and
  shape fitting (section 2.6), and two further known limitations: dominant peaks credited to
  a neighbour in concentrated samples (7.2), and samples a guard declines to quantify (7.3).

### Tools

- `tools/validation_figures.py` draws the archive validation figures from a pooled
  comparison.
- `tools/bench_areas.py` and `tools/bench_quant.py` accept label lists with `--param`, for
  example `--param forced_riders=` to measure with a rule off.

## 1.0.1

Release for archiving. No changes to the software.

## 1.0.0

First public release.
