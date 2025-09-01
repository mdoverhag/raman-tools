import { invoke } from "@tauri-apps/api/core";

export interface BaselineParameters {
  /**
   * Whether to apply denoising before baseline correction
   */
  denoise?: boolean;

  /**
   * Window size for moving average denoising (must be odd)
   * @default 5
   */
  windowSize?: number;

  /**
   * ALS smoothing parameter (larger = smoother baseline)
   * @default 1e7
   */
  lambdaParam?: number;

  /**
   * ALS asymmetry parameter (smaller = more weight to points below baseline)
   * @default 0.01
   */
  p?: number;

  /**
   * ALS order of differences for penalty matrix
   * @default 2
   */
  d?: number;
}

export interface BaselineResult {
  /**
   * The baseline-corrected spectrum
   */
  corrected: number[];

  /**
   * The extracted baseline
   */
  baseline: number[];

  /**
   * The denoised spectrum (if denoising was applied)
   */
  denoised?: number[];
}

/**
 * Apply baseline correction to a spectrum using the ALS algorithm
 * @param spectrum Array of intensity values
 * @param params Optional parameters for baseline correction
 * @returns Promise with corrected spectrum and baseline
 */
export async function applyBaselineCorrection(
  spectrum: number[],
  params: BaselineParameters = {}
): Promise<BaselineResult> {
  const { denoise = false, windowSize = 5, lambdaParam = 1e7, p = 0.01, d = 2 } = params;

  try {
    const result = await invoke<BaselineResult>("apply_baseline_correction", {
      spectrum,
      denoise,
      windowSize,
      lambdaParam,
      p,
      d,
    });

    return result;
  } catch (error) {
    console.error("Baseline correction failed:", error);
    throw new Error(
      `Failed to apply baseline correction: ${error instanceof Error ? error.message : String(error)}`
    );
  }
}

/**
 * Check if the Python runtime is available for baseline correction
 * @returns Promise<boolean> indicating if Python is available
 */
export async function checkPythonRuntime(): Promise<boolean> {
  try {
    return await invoke<boolean>("check_python_runtime");
  } catch (error) {
    console.error("Failed to check Python runtime:", error);
    return false;
  }
}

/**
 * Get information about the Python runtime for diagnostics
 * @returns Promise<string> with Python version and path information
 */
export async function getPythonInfo(): Promise<string> {
  try {
    return await invoke<string>("get_python_info");
  } catch (error) {
    console.error("Failed to get Python info:", error);
    throw error;
  }
}

/**
 * Apply baseline correction to multiple spectra in parallel
 * @param spectra Array of spectra (each spectrum is an array of intensities)
 * @param params Parameters for baseline correction
 * @returns Promise with array of baseline results
 */
export async function batchBaselineCorrection(
  spectra: number[][],
  params: BaselineParameters = {}
): Promise<BaselineResult[]> {
  // Process spectra in parallel with a concurrency limit
  const BATCH_SIZE = 10; // Process 10 spectra at a time
  const results: BaselineResult[] = [];

  for (let i = 0; i < spectra.length; i += BATCH_SIZE) {
    const batch = spectra.slice(i, i + BATCH_SIZE);
    const batchPromises = batch.map((spectrum) => applyBaselineCorrection(spectrum, params));

    const batchResults = await Promise.all(batchPromises);
    results.push(...batchResults);
  }

  return results;
}

/**
 * Helper function to subtract baseline from original spectrum
 * @param original Original spectrum intensities
 * @param baseline Baseline intensities
 * @returns Corrected spectrum
 */
export function subtractBaseline(original: number[], baseline: number[]): number[] {
  if (original.length !== baseline.length) {
    throw new Error("Original and baseline must have the same length");
  }

  return original.map((value, index) => value - baseline[index]);
}

/**
 * Helper function to add baseline back to corrected spectrum
 * @param corrected Corrected spectrum intensities
 * @param baseline Baseline intensities
 * @returns Original spectrum
 */
export function addBaseline(corrected: number[], baseline: number[]): number[] {
  if (corrected.length !== baseline.length) {
    throw new Error("Corrected and baseline must have the same length");
  }

  return corrected.map((value, index) => value + baseline[index]);
}
