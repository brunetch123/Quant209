# Changelog

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
