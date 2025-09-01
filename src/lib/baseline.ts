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
