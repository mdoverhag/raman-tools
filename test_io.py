#!/usr/bin/env python3
"""
Quick test of the I/O functions with real data.
"""

from raman_lib import load_spectra, get_peak, get_color, RAMAN_MOLECULES

# Print molecule configuration
print("Configured Raman molecules:")
for molecule, props in RAMAN_MOLECULES.items():
    print(f"  {molecule}: peak={get_peak(molecule)} cm⁻¹, color={get_color(molecule)}")

print("\n" + "="*60)

# Test loading spectra from a directory
test_dir = "../spectometry-graphs/2025-07-23/MBA/MBA IgG"

print(f"\nAttempting to load spectra from: {test_dir}")

try:
    spectra = load_spectra(test_dir)

    print(f"✓ Successfully loaded {len(spectra)} spectra")

    # Show details of first spectrum
    if spectra:
        first = spectra[0]
        print(f"\nFirst spectrum: {first['filename']}")
        print(f"  Data points: {len(first['wavenumbers'])}")
        print(f"  Wavenumber range: {first['wavenumbers'][0]:.1f} - {first['wavenumbers'][-1]:.1f} cm⁻¹")
        print(f"  Intensity range: {min(first['intensities']):.1f} - {max(first['intensities']):.1f}")

except FileNotFoundError as e:
    print(f"✗ Directory not found: {e}")
    print("\nTo test with your data, edit test_io.py and set test_dir to a valid path")
except Exception as e:
    print(f"✗ Error: {e}")
