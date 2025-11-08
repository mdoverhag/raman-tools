#!/usr/bin/env python3
"""
Test plotting functions by generating reference and sample plots.
"""

from raman_lib import (
    load_spectra,
    apply_baseline_correction,
    calculate_average,
    plot_reference,
    plot_sample,
    create_output_dir,
    ensure_output_subdir
)

# Create output directory for test plots
output = create_output_dir("test-plots")
print(f"Output directory: {output}\n")

# Test 1: Plot a reference spectrum (MBA)
print("="*60)
print("Test 1: Reference Spectrum (MBA)")
print("="*60)

mba_dir = "../spectometry-graphs/2025-07-23/MBA/MBA IgG"
print(f"Loading MBA spectra from: {mba_dir}")

mba_spectra = load_spectra(mba_dir)
print(f"✓ Loaded {len(mba_spectra)} spectra")

# Apply baseline correction
print("Applying baseline correction...")
mba_corrected = [apply_baseline_correction(s) for s in mba_spectra]
print(f"✓ Baseline correction applied")

# Calculate average
print("Calculating average...")
mba_averaged = calculate_average(mba_corrected)
print(f"✓ Average calculated from {mba_averaged['count']} spectra")

# Create plot
refs_dir = ensure_output_subdir(output, "references")
mba_plot_path = f"{refs_dir}/MBA.png"
print(f"Creating plot: {mba_plot_path}")
plot_reference(mba_averaged, molecule="MBA", output_path=mba_plot_path)
print(f"✓ Plot saved!\n")


# Test 2: Plot a sample spectrum (multiplex)
print("="*60)
print("Test 2: Sample Spectrum (Multiplex)")
print("="*60)

# Try to find a multiplex directory
multiplex_dir = "../spectometry-graphs/2025-07-23/Multiplex/Multiplex BSA"
print(f"Loading multiplex spectra from: {multiplex_dir}")

try:
    multiplex_spectra = load_spectra(multiplex_dir)
    print(f"✓ Loaded {len(multiplex_spectra)} spectra")

    # Apply baseline correction
    print("Applying baseline correction...")
    multiplex_corrected = [apply_baseline_correction(s) for s in multiplex_spectra]
    print(f"✓ Baseline correction applied")

    # Calculate average
    print("Calculating average...")
    multiplex_averaged = calculate_average(multiplex_corrected)
    print(f"✓ Average calculated from {multiplex_averaged['count']} spectra")

    # Create plot (multiplex has all three molecules)
    samples_dir = ensure_output_subdir(output, "samples")
    multiplex_plot_path = f"{samples_dir}/Multiplex_BSA.png"
    print(f"Creating plot: {multiplex_plot_path}")
    plot_sample(
        multiplex_averaged,
        sample_name="Multiplex BSA",
        molecules=["MBA", "DTNB", "TFMBA"],
        output_path=multiplex_plot_path
    )
    print(f"✓ Plot saved!\n")

except FileNotFoundError as e:
    print(f"✗ Multiplex directory not found, skipping sample plot test")
    print(f"  (This is OK - we tested the main functionality with references)\n")

print("="*60)
print("✓ Plotting tests complete!")
print(f"✓ Check output in: {output}")
print("="*60)
