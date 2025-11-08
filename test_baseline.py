#!/usr/bin/env python3
"""
Test baseline correction and averaging with real data.
"""

from raman_lib import load_spectra, apply_baseline_correction, calculate_average

# Load spectra from a real directory
test_dir = "../spectometry-graphs/2025-07-23/MBA/MBA IgG"

print(f"Loading spectra from: {test_dir}")
spectra = load_spectra(test_dir)
print(f"✓ Loaded {len(spectra)} spectra\n")

# Apply baseline correction to all spectra
print("Applying ALS baseline correction...")
corrected_spectra = []

for i, spectrum in enumerate(spectra):
    corrected = apply_baseline_correction(spectrum)
    corrected_spectra.append(corrected)

    if i == 0:
        # Show details for first spectrum
        print(f"\nFirst spectrum: {spectrum['filename']}")
        print(f"  Original intensity range: {min(spectrum['intensities']):.1f} - {max(spectrum['intensities']):.1f}")
        print(f"  Corrected intensity range: {min(corrected['corrected']):.1f} - {max(corrected['corrected']):.1f}")
        print(f"  Baseline range: {min(corrected['baseline']):.1f} - {max(corrected['baseline']):.1f}")

print(f"\n✓ Baseline correction applied to {len(corrected_spectra)} spectra\n")

# Calculate average
print("Calculating average across all spectra...")
averaged = calculate_average(corrected_spectra)

print(f"\n✓ Average calculated from {averaged['count']} spectra")
print(f"  Data points: {len(averaged['wavenumbers'])}")
print(f"  Raw average range: {min(averaged['raw_avg']):.1f} - {max(averaged['raw_avg']):.1f}")
print(f"  Corrected average range: {min(averaged['corrected_avg']):.1f} - {max(averaged['corrected_avg']):.1f}")
print(f"  Average std dev (corrected): {sum(averaged['corrected_std']) / len(averaged['corrected_std']):.2f}")

# Find peak around MBA expected peak (1078 cm⁻¹)
mba_peak = 1078
wavenumbers = averaged['wavenumbers']
corrected_avg = averaged['corrected_avg']

# Find closest wavenumber to MBA peak
closest_idx = min(range(len(wavenumbers)), key=lambda i: abs(wavenumbers[i] - mba_peak))
peak_wavenumber = wavenumbers[closest_idx]
peak_intensity = corrected_avg[closest_idx]

print(f"\nIntensity at MBA peak (~{mba_peak} cm⁻¹):")
print(f"  Wavenumber: {peak_wavenumber:.1f} cm⁻¹")
print(f"  Corrected intensity: {peak_intensity:.1f}")

print("\n✓ All tests passed!")
