# Build the Quant 209 release: a standalone application plus a click-to-install wrapper.
#
#     powershell -ExecutionPolicy Bypass -File packaging\build_release.ps1
#
# Produces packaging\release\Quant209-<version>-windows.zip, which a tester extracts and
# installs by double-clicking. Nothing on the target machine needs Python, a terminal, or
# administrator rights.

$ErrorActionPreference = "Stop"
$root = Split-Path -Parent $PSScriptRoot
$python = Join-Path $root ".venv\Scripts\python.exe"
$pyinstaller = Join-Path $root ".venv\Scripts\pyinstaller.exe"

if (-not (Test-Path $pyinstaller)) {
    throw "PyInstaller is not installed in .venv. Run: .venv\Scripts\python.exe -m pip install pyinstaller"
}

$version = & $python -c "import quant209; print(quant209.__version__)"
Write-Host "Building Quant 209 $version"

# 1. Freeze the application into a folder. One folder, not one file: a one-file build
#    unpacks itself on every launch and would break the cold-start budget.
& $pyinstaller (Join-Path $PSScriptRoot "quant209.spec") --noconfirm `
    --distpath (Join-Path $PSScriptRoot "dist") --workpath (Join-Path $PSScriptRoot "build")
if ($LASTEXITCODE -ne 0) { throw "PyInstaller failed" }

# 2. Lay out what the tester receives.
$stage = Join-Path $PSScriptRoot "release\Quant209-$version-windows"
if (Test-Path $stage) { Remove-Item $stage -Recurse -Force }
New-Item -ItemType Directory -Path $stage -Force | Out-Null

Copy-Item (Join-Path $PSScriptRoot "dist\Quant209") (Join-Path $stage "app") -Recurse
Copy-Item (Join-Path $PSScriptRoot "Install Quant 209.vbs") $stage
Copy-Item (Join-Path $PSScriptRoot "Uninstall Quant 209.vbs") $stage

@"
Quant 209 $version
==================

TO INSTALL
  Double-click "Install Quant 209.vbs".
  A shortcut is placed on your Desktop and in the Start Menu.

TO RUN
  Double-click the Quant 209 shortcut.

TO REMOVE
  Double-click "Uninstall Quant 209.vbs".

You do not need Python, a terminal, or administrator rights. The program installs into
your own user profile and can be removed at any time.

If your organisation blocks .vbs files, you can skip installing entirely: open the "app"
folder and double-click Quant209.exe.

CONVERTING RAW AGILENT FILES
  Reading .D acquisitions needs ProteoWizard's msconvert.exe, which is not included here.
  The Convert step will ask for its location once and remember it. A run folder that
  already contains a "Converted" folder of .mzML files does not need it at all.

Full documentation is in README.md in the source repository.
"@ | Set-Content (Join-Path $stage "READ ME FIRST.txt") -Encoding UTF8

# 3. Zip it.
$zip = Join-Path $PSScriptRoot "release\Quant209-$version-windows.zip"
if (Test-Path $zip) { Remove-Item $zip -Force }
Compress-Archive -Path "$stage\*" -DestinationPath $zip -CompressionLevel Optimal

$mb = [math]::Round((Get-Item $zip).Length / 1MB, 1)
Write-Host ""
Write-Host "Release: $zip ($mb MB)"
Write-Host "Give a tester that one file. They extract it and double-click Install Quant 209.vbs."
