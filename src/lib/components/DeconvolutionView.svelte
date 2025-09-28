<script lang="ts">
  import { sampleStore, type Sample, type MoleculePair } from "$lib/stores/samples.svelte";
  import { invoke } from "@tauri-apps/api/core";
  import NormalizationChart from "./NormalizationChart.svelte";

  interface Props {
    sample: Sample | null;
  }

  let { sample }: Props = $props();

  // State for normalization
  let currentRunId = $state<string | null>(null);
  let isNormalizing = $state(false);
  let normalizedData = $state<any>(null);
  let error = $state<string | null>(null);

  // State for deconvolution
  let isDeconvoluting = $state(false);
  let deconvolutionResults = $state<any>(null);

  // Check which singleplex reference samples are available
  let requiredReferences = $derived(() => {
    if (!sample?.moleculePairs) return [];

    // Get unique Raman molecules from the multiplex
    const ramanMolecules = [...new Set(sample.moleculePairs.map((p: MoleculePair) => p.raman))];

    // Check which singleplex samples exist for each Raman molecule
    return ramanMolecules.map((raman) => {
      // Find singleplex samples (exactly one molecule pair with this Raman)
      const singleplexSample = sampleStore.samples.find(
        (s) =>
          s.moleculePairs.length === 1 &&
          s.moleculePairs[0].raman === raman &&
          s.averageCorrected && // Must have average spectrum
          s.averageCorrected.length > 0
      );

      return {
        raman,
        available: !!singleplexSample,
        sampleName: singleplexSample?.name,
        sampleId: singleplexSample?.id,
      };
    });
  });

  let allReferencesAvailable = $derived(
    requiredReferences().length > 0 && requiredReferences().every((r) => r.available)
  );

  async function runNormalization() {
    if (!sample || !allReferencesAvailable) return;

    isNormalizing = true;
    error = null;

    try {
      // Get reference sample IDs
      const referenceSampleIds = requiredReferences()
        .filter((r) => r.sampleId)
        .map((r) => r.sampleId);

      // Create a deconvolution run
      const run = await invoke("create_deconvolution_run", {
        data: {
          multiplexSampleId: sample.id,
          referenceSampleIds: referenceSampleIds,
          wavenumberRange: [1000, 1500], // Default range
        },
      });

      currentRunId = (run as any).id;

      // Calculate normalization
      const result: any = await invoke("calculate_normalization", {
        deconvolutionId: currentRunId,
      });

      // Add reference sample IDs to the result for the chart
      normalizedData = {
        ...result,
        referenceSampleIds,
      };
    } catch (e) {
      error = String(e);
      console.error("Normalization failed:", e);
    } finally {
      isNormalizing = false;
    }
  }

  async function runDeconvolution() {
    if (!currentRunId || !normalizedData) return;

    isDeconvoluting = true;
    error = null;

    try {
      // Perform NNLS deconvolution
      const result = await invoke("perform_deconvolution", {
        deconvolutionId: currentRunId,
      });

      deconvolutionResults = result;
    } catch (e) {
      error = String(e);
      console.error("Deconvolution failed:", e);
    } finally {
      isDeconvoluting = false;
    }
  }
</script>

<div class="bg-gray-800 rounded-lg border border-gray-700 p-6">
  {#if normalizedData}
    <!-- Show visualization when we have normalized data -->
    <div class="space-y-4">
      <NormalizationChart {normalizedData} {sample} {sampleStore} {deconvolutionResults} />

      <!-- Deconvolution section -->
      <div class="bg-gray-900/50 rounded-lg p-4 border border-gray-700">
        <h3 class="text-lg font-semibold text-gray-200 mb-4">NNLS Deconvolution</h3>

        {#if deconvolutionResults}
          <!-- Show deconvolution results -->
          <div class="space-y-4">
            <!-- Contributions -->
            <div class="bg-gray-800/50 rounded-lg p-4">
              <h4 class="text-sm font-semibold text-gray-300 mb-3">Component Contributions</h4>
              <div class="space-y-2">
                {#each Object.entries(deconvolutionResults.contributions || {}) as [name, percentage]}
                  <div class="flex items-center justify-between p-2 bg-gray-900/50 rounded">
                    <span class="text-sm text-gray-300">{name}</span>
                    <span class="text-sm font-mono text-gray-200">
                      {(percentage as number).toFixed(2)}%
                    </span>
                  </div>
                {/each}
              </div>
            </div>

            <!-- Metrics -->
            <div class="bg-gray-800/50 rounded-lg p-4">
              <h4 class="text-sm font-semibold text-gray-300 mb-3">Quality Metrics</h4>
              <div class="grid grid-cols-2 gap-4 text-sm">
                <div>
                  <span class="text-gray-400">R²:</span>
                  <span class="text-gray-200 ml-2">
                    {deconvolutionResults.metrics?.rSquared.toFixed(4)}
                  </span>
                </div>
                <div>
                  <span class="text-gray-400">RMSE:</span>
                  <span class="text-gray-200 ml-2">
                    {deconvolutionResults.metrics?.rmse.toFixed(4)}
                  </span>
                </div>
              </div>
            </div>
          </div>
        {:else}
          <!-- Show deconvolution button -->
          <div class="flex flex-col items-center">
            <p class="text-gray-400 mb-4">
              Ready to perform NNLS deconvolution on normalized spectra
            </p>
            <button
              onclick={runDeconvolution}
              disabled={isDeconvoluting}
              class="px-6 py-2 bg-purple-600 hover:bg-purple-700 disabled:bg-gray-600
                     text-white rounded-lg transition-colors flex items-center justify-center gap-2"
            >
              {#if isDeconvoluting}
                <svg class="animate-spin h-5 w-5" fill="none" viewBox="0 0 24 24">
                  <circle
                    class="opacity-25"
                    cx="12"
                    cy="12"
                    r="10"
                    stroke="currentColor"
                    stroke-width="4"
                  ></circle>
                  <path
                    class="opacity-75"
                    fill="currentColor"
                    d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"
                  ></path>
                </svg>
                Deconvoluting...
              {:else}
                Run Deconvolution
              {/if}
            </button>
          </div>
        {/if}
      </div>
    </div>
  {:else}
    <!-- Show status and controls when no data -->
    <div class="flex flex-col items-center justify-center min-h-[24rem]">
      <svg
        class="h-16 w-16 {allReferencesAvailable ? 'text-purple-400' : 'text-gray-500'} mb-4"
        fill="none"
        stroke="currentColor"
        viewBox="0 0 24 24"
      >
        <path
          stroke-linecap="round"
          stroke-linejoin="round"
          stroke-width="2"
          d="M4 4v5h.582m15.356 2A8.001 8.001 0 004.582 9m0 0H9m11 11v-5h-.581m0 0a8.003 8.003 0 01-15.357-2m15.357 2H15"
        />
      </svg>
      <h2
        class="text-2xl font-bold {allReferencesAvailable
          ? 'text-purple-300'
          : 'text-gray-400'} mb-4"
      >
        Deconvolution Analysis
      </h2>

      {#if sample?.moleculePairs}
        <div class="w-full max-w-md">
          <!-- Reference samples status -->
          <div class="mb-6">
            <h3 class="text-sm font-medium text-gray-300 mb-3">Required Reference Samples:</h3>
            <div class="space-y-2">
              {#each requiredReferences() as ref}
                <div class="flex items-center justify-between p-3 bg-gray-700/30 rounded-lg">
                  <div class="flex items-center gap-3">
                    {#if ref.available}
                      <svg
                        class="w-5 h-5 text-green-400"
                        fill="none"
                        stroke="currentColor"
                        viewBox="0 0 24 24"
                      >
                        <path
                          stroke-linecap="round"
                          stroke-linejoin="round"
                          stroke-width="2"
                          d="M5 13l4 4L19 7"
                        />
                      </svg>
                    {:else}
                      <svg
                        class="w-5 h-5 text-red-400"
                        fill="none"
                        stroke="currentColor"
                        viewBox="0 0 24 24"
                      >
                        <path
                          stroke-linecap="round"
                          stroke-linejoin="round"
                          stroke-width="2"
                          d="M6 18L18 6M6 6l12 12"
                        />
                      </svg>
                    {/if}
                    <span class="text-gray-200 font-medium">{ref.raman}</span>
                  </div>
                  <span class="text-sm {ref.available ? 'text-gray-400' : 'text-red-400'}">
                    {ref.available ? ref.sampleName : "Missing"}
                  </span>
                </div>
              {/each}
            </div>
          </div>

          <!-- Status message and normalization -->
          {#if normalizedData}
            <!-- Show normalized data -->
            <div class="p-4 bg-green-900/20 border border-green-600/30 rounded-lg mb-4">
              <p class="text-green-300 text-center">✓ Normalization complete</p>
              <p class="text-gray-400 text-sm text-center mt-2">
                Spectra normalized in range 1000-1500 cm⁻¹
              </p>
            </div>
          {:else if allReferencesAvailable}
            <div class="p-4 bg-purple-900/20 border border-purple-600/30 rounded-lg">
              <p class="text-purple-300 text-center">✓ All reference samples available</p>
              <p class="text-gray-400 text-sm text-center mt-2">Ready for deconvolution analysis</p>

              <!-- Normalization button -->
              <button
                onclick={runNormalization}
                disabled={isNormalizing}
                class="mt-4 w-full px-4 py-2 bg-purple-600 hover:bg-purple-700 disabled:bg-gray-600
                     text-white rounded-lg transition-colors flex items-center justify-center gap-2"
              >
                {#if isNormalizing}
                  <svg class="animate-spin h-5 w-5" fill="none" viewBox="0 0 24 24">
                    <circle
                      class="opacity-25"
                      cx="12"
                      cy="12"
                      r="10"
                      stroke="currentColor"
                      stroke-width="4"
                    ></circle>
                    <path
                      class="opacity-75"
                      fill="currentColor"
                      d="M4 12a8 8 0 018-8V0C5.373 0 0 5.373 0 12h4zm2 5.291A7.962 7.962 0 014 12H0c0 3.042 1.135 5.824 3 7.938l3-2.647z"
                    ></path>
                  </svg>
                  Normalizing...
                {:else}
                  Run Normalization
                {/if}
              </button>
            </div>

            {#if error}
              <div class="mt-4 p-3 bg-red-900/20 border border-red-600/30 rounded-lg">
                <p class="text-red-300 text-sm">Error: {error}</p>
              </div>
            {/if}
          {:else}
            <div class="p-4 bg-red-900/20 border border-red-600/30 rounded-lg">
              <p class="text-red-300 text-center">⚠ Missing reference samples</p>
              <p class="text-gray-400 text-sm text-center mt-2">
                Please import singleplex samples for all Raman molecules
              </p>
            </div>
          {/if}
        </div>
      {:else}
        <p class="text-gray-400">No molecule pairs defined</p>
      {/if}
    </div>
  {/if}
</div>
