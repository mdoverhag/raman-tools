<script lang="ts">
  import { sampleStore } from "$lib/stores/samples.svelte";
  import { BeakerIcon, PlusIcon, TrashIcon } from "@babeard/svelte-heroicons/mini";
  import { getMoleculeClasses } from "$lib/utils/molecule-colors";

  async function handleCreateSample() {
    const sample = await sampleStore.createSample("New Sample");
    if (sample) {
      sampleStore.selectSample(sample.id);
    }
  }
</script>

<div class="flex h-full w-64 flex-col border-r border-gray-700 bg-gray-900">
  <div class="flex-1 overflow-y-auto">
    <div class="sticky top-0 z-10 bg-gray-900 px-4 pt-6 pb-2">
      <h2 class="text-sm font-medium text-gray-500">Samples</h2>
    </div>

    <div class="space-y-0.5 px-2 pb-2">
      {#each sampleStore.samples as sample}
        <div
          class="group flex w-full items-start rounded px-2 py-2 text-left text-sm hover:bg-gray-800 {sampleStore.selectedSampleId ===
          sample.id
            ? 'bg-blue-900/30 text-blue-400'
            : 'text-gray-300'}"
        >
          <button
            onclick={() => sampleStore.selectSample(sample.id)}
            class="flex min-w-0 flex-1 items-start text-left"
          >
            <BeakerIcon class="mr-2 mt-0.5 h-4 w-4 flex-shrink-0" />
            <div class="min-w-0 flex-1">
              <div class="font-medium">{sample.name}</div>
              <div class="mt-1 space-y-1">
                <div class="flex flex-wrap gap-1">
                  {#if sample.ramanMolecules.length > 0}
                    {#each sample.ramanMolecules as molecule}
                      <span
                        class="inline-flex items-center rounded-md px-1.5 py-0.5 text-[10px] font-medium {getMoleculeClasses(
                          molecule
                        )}"
                      >
                        {molecule}
                      </span>
                    {/each}
                  {:else}
                    <span
                      class="inline-flex items-center rounded-md border border-dashed border-gray-700 px-1.5 py-0.5 text-[10px] text-gray-600"
                    >
                      No Raman
                    </span>
                  {/if}
                </div>
                <div class="flex flex-wrap gap-1">
                  {#if sample.targetMolecules.length > 0}
                    {#each sample.targetMolecules as molecule}
                      <span
                        class="inline-flex items-center rounded-md px-1.5 py-0.5 text-[10px] font-medium {getMoleculeClasses(
                          molecule
                        )}"
                      >
                        {molecule}
                      </span>
                    {/each}
                  {:else}
                    <span
                      class="inline-flex items-center rounded-md border border-dashed border-gray-700 px-1.5 py-0.5 text-[10px] text-gray-600"
                    >
                      No Target
                    </span>
                  {/if}
                </div>
              </div>
              {#if sample.spectrumIds.length > 0}
                <div class="mt-0.5 text-xs text-gray-600">
                  {sample.spectrumIds.length} spectra
                </div>
              {/if}
            </div>
          </button>
          <button
            onclick={async (e) => {
              e.stopPropagation();
              await sampleStore.deleteSample(sample.id);
            }}
            class="ml-1 rounded p-0.5 opacity-0 hover:bg-red-900/50 hover:text-red-400 group-hover:opacity-100"
            type="button"
          >
            <TrashIcon class="h-3 w-3" />
          </button>
        </div>
      {/each}

      {#if sampleStore.samples.length === 0}
        <div class="px-2 py-8 text-center text-sm text-gray-500">
          No samples yet. Create one to get started.
        </div>
      {/if}
    </div>
  </div>

  <div class="border-t border-gray-700 p-3">
    <button
      onclick={handleCreateSample}
      class="flex w-full items-center justify-center gap-2 rounded-md bg-gray-800 px-3 py-2 text-sm font-medium text-gray-300 hover:bg-gray-700 hover:text-gray-100"
    >
      <PlusIcon class="h-4 w-4" />
      New Sample
    </button>
  </div>
</div>
