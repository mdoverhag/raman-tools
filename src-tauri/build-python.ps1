#!/usr/bin/env pwsh
# Build script for Windows - creates Python bundle for production builds

$ErrorActionPreference = "Stop"

# Script configuration
$PYTHON_VERSION = "3.13"
$SCRIPT_DIR = Split-Path -Parent $MyInvocation.MyCommand.Path
$BUNDLE_DIR = Join-Path $SCRIPT_DIR "resources\python\windows"

Write-Host "Building Python bundle for Windows..." -ForegroundColor Cyan
Write-Host "Target directory: $BUNDLE_DIR"

# Check if uv is installed
if (-not (Get-Command uv -ErrorAction SilentlyContinue)) {
    Write-Host "❌ Error: uv is not installed" -ForegroundColor Red
    Write-Host "Please install uv first:" -ForegroundColor Yellow
    Write-Host "  powershell -ExecutionPolicy ByPass -c `"irm https://astral.sh/uv/install.ps1 | iex`"" -ForegroundColor White
    Write-Host "Or via winget:" -ForegroundColor Yellow
    Write-Host "  winget install --id=astral-sh.uv -e" -ForegroundColor White
    exit 1
}

Write-Host "✅ Found uv: $(Get-Command uv | Select-Object -ExpandProperty Source)" -ForegroundColor Green

# Clean up existing bundle
if (Test-Path $BUNDLE_DIR) {
    Write-Host "Cleaning existing bundle..." -ForegroundColor Yellow
    Remove-Item -Recurse -Force $BUNDLE_DIR
}

# Create bundle directory
New-Item -ItemType Directory -Force -Path $BUNDLE_DIR | Out-Null

# Create virtual environment with Python 3.13
Write-Host "Creating virtual environment with Python $PYTHON_VERSION..." -ForegroundColor Cyan
& uv venv --python $PYTHON_VERSION --relocatable "$BUNDLE_DIR"
if ($LASTEXITCODE -ne 0) {
    Write-Host "❌ Failed to create virtual environment" -ForegroundColor Red
    exit 1
}

# Install required packages from requirements.txt
Write-Host "Installing packages from requirements.txt..." -ForegroundColor Cyan
$REQUIREMENTS = Join-Path $SCRIPT_DIR "python\requirements.txt"
if (-not (Test-Path $REQUIREMENTS)) {
    Write-Host "❌ Error: requirements.txt not found at $REQUIREMENTS" -ForegroundColor Red
    exit 1
}

& uv pip install --python "$BUNDLE_DIR" -r "$REQUIREMENTS"
if ($LASTEXITCODE -ne 0) {
    Write-Host "❌ Failed to install packages" -ForegroundColor Red
    exit 1
}

# Copy Python scripts
Write-Host "Copying Python scripts..." -ForegroundColor Cyan
$BASELINE_SCRIPT = Join-Path $SCRIPT_DIR "python\baseline_correction.py"
if (Test-Path $BASELINE_SCRIPT) {
    Copy-Item $BASELINE_SCRIPT -Destination $BUNDLE_DIR
    Write-Host "✅ Copied baseline_correction.py" -ForegroundColor Green
} else {
    Write-Host "⚠️  Warning: baseline_correction.py not found" -ForegroundColor Yellow
}

# Verify installation
Write-Host "`nVerifying installation..." -ForegroundColor Cyan
$PYTHON_EXE = Join-Path $BUNDLE_DIR "Scripts\python.exe"

if (-not (Test-Path $PYTHON_EXE)) {
    Write-Host "❌ Python executable not found at $PYTHON_EXE" -ForegroundColor Red
    exit 1
}

# Test imports
Write-Host "Testing Python imports..." -ForegroundColor Cyan
$TestScript = @'
import sys
import numpy as np
import scipy
print(f"Python {sys.version}")
print(f"NumPy {np.__version__}")
print(f"SciPy {scipy.__version__}")
'@

$TestScript | & $PYTHON_EXE 2>&1
if ($LASTEXITCODE -ne 0) {
    Write-Host "❌ Failed to import required packages" -ForegroundColor Red
    exit 1
}

Write-Host "`n✅ Python bundle built successfully!" -ForegroundColor Green
Write-Host "Bundle location: $BUNDLE_DIR" -ForegroundColor Cyan
Write-Host "Python executable: $PYTHON_EXE" -ForegroundColor Cyan

# Display bundle size
$BundleSize = (Get-ChildItem -Recurse $BUNDLE_DIR | Measure-Object -Property Length -Sum).Sum / 1MB
Write-Host "Bundle size: $([math]::Round($BundleSize, 2)) MB" -ForegroundColor Cyan