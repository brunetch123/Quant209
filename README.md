<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/logo-dark.png">
    <img src="docs/images/logo.png" alt="Quant 209" height="58">
  </picture>
</h1>

Quant 209 quantifies all 209 PCB congeners from Agilent triple-quadrupole GC-MS/MS data. It
reads an analysis run, integrates every MRM transition, and reports congener masses, surrogate
recoveries and QA/QC flags.

Quant 209 is a Windows desktop application. A command-line interface is also available for
batch and scripted use.

![The Review step, showing one sample's penta transition with integrations drawn over the
acquired trace](docs/images/review-samples.png)

## What it does

- Converts Agilent `.D` acquisitions to mzML through ProteoWizard's MSConvert.
- Integrates 185 chromatographic peaks per file. These are 173 native congeners and
  co-elution groups plus 12 labelled standards. Retention times come from the run's own
  calibration injections.
- Calculates response factors, congener masses, per-homolog surrogate recoveries and
  recovery-corrected masses.
- Applies a blank-derived limit of quantification and flags results that need review.
- Writes results to an Excel workbook. Where a manually quantified workbook is present in the
  run folder, it also writes a report comparing the two.

Each integrated peak records its retention window, baseline, boundaries and integration
method. Any reported value can be traced back to the chromatogram it came from.

## Install

### The application

Download `Quant209-<version>-windows.zip` from the
[Releases page](https://github.com/brunetch123/Quant209/releases), extract it, and
double-click **Install Quant 209.vbs**. Shortcuts are placed on the Desktop and in the Start
Menu. What changed in each version is listed in [`CHANGELOG.md`](CHANGELOG.md).

The installer does not require administrator rights or a Python installation. It writes to
`%LOCALAPPDATA%\Programs\Quant209` within the current user profile. **Uninstall Quant
209.vbs** removes the program and its shortcuts. Data files and exported reports are left in
place.

Windows shows a "Windows protected your PC" warning on first launch because the executable is
not code-signed. Click **More info**, then **Run anyway**.

If `.vbs` files are blocked on your system, open the extracted `app` folder and run
`Quant209.exe` directly. No installation is needed.

### ProteoWizard

Reading raw `.D` acquisitions requires ProteoWizard's `msconvert.exe`, which is not bundled.
The Convert step asks for its location on first use and stores it. A run folder that already
contains a `Converted` folder of `.mzML` files does not need it.

### From source

Python 3.11 or later is required. Development and validation use 3.13.

```
git clone https://github.com/brunetch123/Quant209.git
cd quant209
python -m venv .venv
.venv\Scripts\activate
pip install -e ".[app,report,dev]"
```

| extra | brings | needed for |
|---|---|---|
| *(base)* | numpy, pandas, openpyxl | the pipeline and the Excel report |
| `app` | PySide6-Essentials, pyqtgraph | the desktop application |
| `report` | matplotlib | the comparison report's figures |
| `dev` | pytest, pytest-qt | the test suite |

To build the installer:

```
pip install pyinstaller
powershell -ExecutionPolicy Bypass -File packaging\build_release.ps1
```

## Demonstration data

No laboratory data ships with this repository. To generate a run folder to open:

```
.venv\Scripts\python.exe tools\make_demo_run.py --out .
```

This writes `DEMO01-06_20240115/`, a fictional run containing a calibration injection, a
reference injection, a blank, a hexane rinse and six samples. The chromatograms are
synthesised from the canonical peak table, so every value in them is invented. The run is
otherwise valid and all six workflow steps operate on it. The screenshots in this file were
taken from it.

## Using the application

Six steps, in order down the left-hand side.

1. **Select.** Choose the run folder, or select files individually.
2. **Assign.** Each file's role (calibration, reference, sample, blank, hexane rinse) is read
   from its filename. Incorrect assignments can be changed here.
3. **Convert.** Runs MSConvert on the raw acquisitions. Files that already have an mzML are
   skipped. Truncated files are converted again. A run that arrives already converted
   skips this step.
4. **Process.** Integrates and quantifies each file, reporting progress per file.
5. **Review.** Shows the run verdict, per-sample chromatograms with integrations drawn on the
   trace, response factors, a recovery heat map and a table of flags.
6. **Export.** Writes the Excel workbook, and the comparison report where a manually
   quantified workbook is present.

A run folder that already contains converted mzML needs only Select and Process.

## Command line

The command line runs the same pipeline without the interface.

```
quant209 run <run-folder> --out output\
quant209 run <run-folder> --from-raw --msconvert "C:\ProteoWizard\msconvert.exe"
quant209 convert <run-folder> --msconvert "C:\ProteoWizard\msconvert.exe"
quant209 validate <run-folder> --out output\
```

`run` writes the Excel workbook and prints a summary. `validate` also compares the result
against a manually quantified workbook in the same folder. Exit status is 0 on success, 1 on
failure and 2 on a usage error.

## Validation

Quant 209 has been compared against manual integrations across 49 analysis runs and 41,376
peak comparisons, spanning four years, several analysts and two instruments. The laboratory
archive holds 64 run folders in total. [`docs/validation/`](docs/validation/) accounts for the
fifteen that are not scored.

| | median relative difference | within ±25% |
|---|---|---|
| Peak areas | 1.9% | 91.8% |
| Congener masses | 3.9% | 87.7% |
| Surrogate recoveries | 0.021 mean absolute error | |

No run in the archive disagrees by more than 3.6% at the median. Agreement is best on large,
well-resolved peaks, with a median of 0.5% above 10⁶ counts·s. It degrades near the detection
limit, where manual integration is also imprecise and results are censored.

[`docs/validation/`](docs/validation/) holds the full comparison.
[`docs/METHODS.md`](docs/METHODS.md) describes every rule the pipeline applies and states its
known limitations. Four congeners are integrated by named rules where the trace alone does not
find them, and each such peak carries the `forced` flag (METHODS section 2.7). In a
concentrated sample, a dominant peak that drifts into a minor neighbour's retention-time window
is given back to it on the strength of published Aroclor profiles, and carries the `displaced`
flag. Two comparable congeners merged into one envelope are split by a two-peak fit and carry
the `apportioned` flag (section 2.8). Two neighbouring congeners a run's calibration does not
separate have their retention times measured there by a two-peak fit, and wherever they share
one envelope it is split at that spacing, flagged `unresolved` (section 2.9). A congener's
copy of itself on a lower homolog's channel, left by losing HCl or Cl2 in the ion source, is
kept out of its neighbours' windows or subtracted where it sits on an analyte's retention time
(other than PCB 58's and PCB 159's), and the peaks it changed are flagged `ghost` (section
2.10). A file whose labelled standards show its chromatography drifted from its calibration's
is searched where the drift puts each analyte, flagged `realigned` (section 3.3). A file
without a usable internal standard gets an error and no mass (section 4.2).

## Layout

```
quant209/
  core/        integration and quantification (no interface code)
    mzml         SRM mzML reader
    config       the 185-analyte canonical peak table
    peaks        integration primitives
    calibrate    retention-time alignment and anchors
    integrate    peak integration with co-elution handling
    deconvolve   shape fitting for peaks no boundary can separate
    quantify     response factors, masses, recoveries
    qaqc         limit of quantification, censoring, flags
    pipeline     end-to-end orchestration
  convert/     MSConvert wrapper
  report/      Excel workbook, comparison report, workbook reader
  app/         the desktop application
  cli.py       command line
tools/         batch conversion, benchmarking, measurement, validation figures
packaging/     installer and build scripts
```

`core/` and `report/` do not import from `app/`. The analysis therefore runs headless and can
be tested without a display.

## Tests

```
.venv\Scripts\python.exe -m pytest -q
```

The suite contains 462 tests. Tests that require instrument data skip when it is absent, so a
fresh clone passes without an archive. Coverage includes the numerical primitives against
analytic values, the calculation chain against hand-computed figures, the desktop application
offscreen, and the packaging configuration.

To re-run the archive comparison:

```
.venv\Scripts\python.exe tools\master_bench.py --out docs\validation\master
.venv\Scripts\python.exe tools\master_report.py --data docs\validation\master
.venv\Scripts\python.exe tools\validation_figures.py --data docs\validation\master
```

The outputs of all three name the laboratory's runs and samples, so none of them is committed.

## Licence

Quant 209 is released under the MIT licence. See [`LICENSE`](LICENSE).

ProteoWizard's `msconvert.exe` is third-party software under its own licence and is not
bundled with Quant 209.

## A note on data

No laboratory data is distributed with this repository and no real sample or run identifier
appears in it. Tests that compare against a real archive read the names they need from
`tests/dataset.json`, which is local and untracked. Without that file those tests skip and the
remainder of the suite runs normally. See [`tests/dataset.py`](tests/dataset.py).
