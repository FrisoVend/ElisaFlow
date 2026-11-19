# Build the application and offline manual in an isolated environment.
[CmdletBinding()]
param([string]$Python)

$ErrorActionPreference = "Stop"
$projectDirectory = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
if (-not $Python) {
    $Python = Join-Path $projectDirectory ".venv\Scripts\python.exe"
}
$environmentDirectory = Join-Path $projectDirectory "build\.venv"
$buildPython = Join-Path $environmentDirectory "Scripts\python.exe"
$documentationDirectory = Join-Path $projectDirectory "build\docs"
$distributionDirectory = Join-Path $projectDirectory ("dist\build_" + (Get-Date -Format "yyyyMMdd_HHmmss_fff"))
$bundleDirectory = Join-Path $distributionDirectory "ElisaFlow"
$previousInstallerCache = $env:PYINSTALLER_CONFIG_DIR
$previousPipConfig = $env:PIP_CONFIG_FILE

function Invoke-BuildStep {
    param([string]$Executable, [string[]]$Arguments)
    & $Executable @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "Build step failed with exit code ${LASTEXITCODE}: $Executable $Arguments"
    }
}

function Remove-BuildDirectory {
    param([string]$Directory)
    if (-not (Test-Path -LiteralPath $Directory)) { return }
    $buildDirectory = Join-Path $projectDirectory "build"
    $resolvedDirectory = (Resolve-Path -LiteralPath $Directory).Path
    if ((Split-Path -Parent $resolvedDirectory) -ne $buildDirectory -or
        $resolvedDirectory -ne [IO.Path]::GetFullPath($Directory) -or
        ((Get-Item -LiteralPath $buildDirectory).Attributes -band [IO.FileAttributes]::ReparsePoint) -or
        ((Get-Item -LiteralPath $resolvedDirectory).Attributes -band [IO.FileAttributes]::ReparsePoint)) {
        throw "Unexpected build directory: $resolvedDirectory"
    }
    Remove-Item -LiteralPath $resolvedDirectory -Recurse -Force
}

Push-Location $projectDirectory
try {
    # Ignore machine-wide package indexes and settings for this dedicated build.
    $env:PIP_CONFIG_FILE = "nul"
    Invoke-BuildStep -Executable $Python -Arguments @("-I", "-m", "venv", $environmentDirectory)
    Invoke-BuildStep -Executable $buildPython -Arguments @(
        "-I", "-m", "pip", "--isolated", "install", "--disable-pip-version-check", "--no-cache-dir",
        "-r", "requirements-build.txt"
    )
    # Discard stale modules and wheel metadata after source or version changes.
    Remove-BuildDirectory -Directory (Join-Path $projectDirectory "build\lib")
    Get-ChildItem -LiteralPath (Join-Path $projectDirectory "build") -Directory -Filter "bdist.*" | ForEach-Object {
        Remove-BuildDirectory -Directory $_.FullName
    }
    Invoke-BuildStep -Executable $buildPython -Arguments @(
        "-I", "-m", "pip", "--isolated", "install", "--disable-pip-version-check", "--no-cache-dir",
        "--no-deps", "--no-build-isolation", "."
    )
    Invoke-BuildStep -Executable $buildPython -Arguments @("-I", "-m", "pip", "check")
    Remove-BuildDirectory -Directory $documentationDirectory
    Invoke-BuildStep -Executable $buildPython -Arguments @(
        "-I", "-m", "sphinx", "-E", "-W", "-b", "html", "docs", $documentationDirectory
    )
    $env:PYINSTALLER_CONFIG_DIR = Join-Path $projectDirectory "build\pyinstaller-cache"
    Invoke-BuildStep -Executable $buildPython -Arguments @(
        "-I", "-m", "PyInstaller", "--clean", "--distpath", $distributionDirectory, "elisaflow.spec"
    )
    Copy-Item -LiteralPath "requirements-build.txt" -Destination $bundleDirectory
    Write-Host "Copy this complete folder to the USB stick: $bundleDirectory"
}
finally {
    $env:PYINSTALLER_CONFIG_DIR = $previousInstallerCache
    $env:PIP_CONFIG_FILE = $previousPipConfig
    Pop-Location
}
