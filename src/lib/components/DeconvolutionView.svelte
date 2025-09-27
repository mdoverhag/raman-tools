<script lang="ts">
  import { sampleStore, type Sample, type MoleculePair } from "$lib/stores/samples.svelte";

  interface Props {
    sample: Sample | null;
  }

  let { sample }: Props = $props();

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
      };
    });
  });

  let allReferencesAvailable = $derived(
    requiredReferences().length > 0 && requiredReferences().every((r) => r.available)
  );
</script>

<div class="bg-gray-800 rounded-lg border border-gray-700 p-6">
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
      class="text-2xl font-bold {allReferencesAvailable ? 'text-purple-300' : 'text-gray-400'} mb-4"
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

        <!-- Status message -->
        {#if allReferencesAvailable}
          <div class="p-4 bg-purple-900/20 border border-purple-600/30 rounded-lg">
            <p class="text-purple-300 text-center">✓ All reference samples available</p>
            <p class="text-gray-400 text-sm text-center mt-2">Ready for deconvolution analysis</p>
          </div>
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
</div>
