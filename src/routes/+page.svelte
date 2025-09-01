<script lang="ts">
  import { listen } from "@tauri-apps/api/event";
  import { invoke } from "@tauri-apps/api/core";
  import { onMount } from "svelte";
  import SpectrumChart from "$lib/SpectrumChart.svelte";
  import { applyBaselineCorrection, type BaselineResult } from "$lib/baseline";

  interface Spectrum {
    filename: string;
    filepath: string;
    wavenumbers: number[];
    intensities: number[];
    baseline?: number[];
    corrected?: number[];
  }

  let isDragging = $state(false);
  let spectra = $state<Spectrum[]>([]);
  let isLoading = $state(false);
  let error = $state<string | null>(null);
  let selectedSpectrum = $state<Spectrum | null>(null);
  let baselineParams = $state({
    denoise: true,
    windowSize: 5,
    lambdaParam: 1e7,
    p: 0.01,
    d: 2,
  });

  async function handleFileDrop(paths: string[]) {
    console.log("Processing files:", paths);
    isLoading = true;
    error = null;

    try {
      // Filter for .txt files only
      const txtFiles = paths.filter((path) => path.toLowerCase().endsWith(".txt"));

      if (txtFiles.length === 0) {
        error = "No .txt files found in dropped files";
        return;
      }

      // Call Rust command to parse the files
      const parsedSpectra = await invoke<Spectrum[]>("parse_spectrum_files", {
        filepaths: txtFiles,
      });

      console.log("Parsed spectra:", parsedSpectra);

      // Check if some files were skipped
      const skippedCount = txtFiles.length - parsedSpectra.length;
      if (skippedCount > 0) {
        error = `Warning: ${skippedCount} file(s) contained no valid spectrum data`;
      }

      if (parsedSpectra.length === 0) {
        error = "No valid spectrum data found in any of the files";
      } else {
        // Apply baseline correction to each spectrum
        console.log("Applying baseline correction to spectra...");
        const correctedSpectra = await Promise.all(
          parsedSpectra.map(async (spectrum) => {
            try {
              const result = await applyBaselineCorrection(spectrum.intensities, baselineParams);
              return {
                ...spectrum,
                baseline: result.baseline,
                corrected: result.corrected,
              };
            } catch (e) {
              console.error(`Failed to apply baseline correction to ${spectrum.filename}:`, e);
              // Return spectrum without baseline correction on error
              return spectrum;
            }
          })
        );
        spectra = correctedSpectra;
      }
    } catch (e) {
      console.error("Error parsing files:", e);
      error = e instanceof Error ? e.message : String(e);
    } finally {
      isLoading = false;
    }
  }

  onMount(() => {
    // Listen for file drop events from Tauri v2
    const unlisten = listen<{ paths: string[] }>("tauri://drag-drop", (event) => {
      isDragging = false;
      if (event.payload?.paths) {
        handleFileDrop(event.payload.paths);
      }
    });

    // Listen for drag over events
    const unlistenHover = listen("tauri://drag-over", () => {
      isDragging = true;
    });

    // Listen for drag leave
    const unlistenCancelled = listen("tauri://drag-leave", () => {
      isDragging = false;
    });

    // Cleanup listeners
    return () => {
      unlisten.then((fn) => fn());
      unlistenHover.then((fn) => fn());
      unlistenCancelled.then((fn) => fn());
    };
  });
</script>

<main class="min-h-screen bg-gray-900 text-gray-100">
  <div class="max-w-7xl mx-auto p-6">
    <h1 class="text-3xl font-bold text-gray-100 mb-8">Raman Spectrum Viewer</h1>

    <!-- Drop Zone -->
    <div
      class="border-2 border-dashed rounded-lg p-12 text-center transition-all duration-300 {isDragging
        ? 'border-blue-500 bg-blue-500/10'
        : 'border-gray-600 bg-gray-800/50 hover:border-gray-500'}"
    >
      <div class="space-y-2">
        <svg
          class="mx-auto h-12 w-12 transition-colors duration-300 {isDragging
            ? 'text-blue-400'
            : 'text-gray-500'}"
          stroke="currentColor"
          fill="none"
          viewBox="0 0 48 48"
        >
          <path
            d="M28 8H12a4 4 0 00-4 4v20m32-12v8m0 0v8a4 4 0 01-4 4H12a4 4 0 01-4-4v-4m32-4l-3.172-3.172a4 4 0 00-5.656 0L28 28M8 32l9.172-9.172a4 4 0 015.656 0L28 28m0 0l4 4m4-24h8m-4-4v8m-12 4h.02"
            stroke-width="2"
            stroke-linecap="round"
            stroke-linejoin="round"
          />
        </svg>
        {#if isLoading}
          <p class="text-gray-400">Processing files...</p>
        {:else if isDragging}
          <p class="text-blue-400 font-medium">Drop files here...</p>
        {:else}
          <p class="text-gray-400">Drag and drop spectrum files (.txt) here</p>
        {/if}
      </div>
    </div>

    <!-- Error Message -->
    {#if error}
      <div class="mt-4 bg-red-900/50 border border-red-700 rounded-lg p-4">
        <p class="text-red-400">{error}</p>
      </div>
    {/if}

    <!-- Content Grid -->
    {#if spectra.length > 0}
      <div class="grid grid-cols-1 lg:grid-cols-[350px_1fr] gap-6 mt-8">
        <!-- File List -->
        <div class="bg-gray-800 rounded-lg border border-gray-700 overflow-hidden">
          <div class="p-4 border-b border-gray-700 bg-gray-800/50">
            <h2 class="text-lg font-semibold text-gray-100">
              Loaded Spectra ({spectra.length} files)
            </h2>
          </div>
          <div class="max-h-[600px] overflow-y-auto">
            <ul class="divide-y divide-gray-700">
              {#each spectra as spectrum}
                <li>
                  <button
                    class="w-full text-left px-4 py-3 hover:bg-gray-700/50 transition-colors {selectedSpectrum ===
                    spectrum
                      ? 'bg-blue-900/30 border-l-2 border-blue-500'
                      : ''}"
                    onclick={() => (selectedSpectrum = spectrum)}
                  >
                    <div class="font-medium text-gray-200">
                      {spectrum.filename}
                    </div>
                    <div class="text-sm text-gray-500">
                      {spectrum.wavenumbers.length} data points
                    </div>
                  </button>
                </li>
              {/each}
            </ul>
          </div>
        </div>

        <!-- Chart Section -->
        {#if selectedSpectrum}
          <div class="bg-gray-800 rounded-lg border border-gray-700 p-6">
            <SpectrumChart
              wavenumbers={selectedSpectrum.wavenumbers}
              intensities={selectedSpectrum.intensities}
              baseline={selectedSpectrum.baseline}
              corrected={selectedSpectrum.corrected}
              title={selectedSpectrum.filename}
            />
          </div>
        {:else}
          <div
            class="bg-gray-800 rounded-lg border border-gray-700 p-12 flex items-center justify-center"
          >
            <div class="text-center">
              <svg
                class="mx-auto h-12 w-12 text-gray-600 mb-4"
                fill="none"
                stroke="currentColor"
                viewBox="0 0 24 24"
              >
                <path
                  stroke-linecap="round"
                  stroke-linejoin="round"
                  stroke-width="2"
                  d="M9 19v-6a2 2 0 00-2-2H5a2 2 0 00-2 2v6a2 2 0 002 2h2a2 2 0 002-2zm0 0V9a2 2 0 012-2h2a2 2 0 012 2v10m-6 0a2 2 0 002 2h2a2 2 0 002-2m0 0V5a2 2 0 012-2h2a2 2 0 012 2v14a2 2 0 01-2 2h-2a2 2 0 01-2-2z"
                />
              </svg>
              <p class="text-gray-500">Select a spectrum to view</p>
            </div>
          </div>
        {/if}
      </div>
    {/if}
  </div>
</main>
