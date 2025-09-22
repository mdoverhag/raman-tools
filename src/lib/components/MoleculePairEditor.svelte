<script lang="ts">
  import { PlusIcon, XMarkIcon } from "@babeard/svelte-heroicons/mini";
  import { sampleStore } from "$lib/stores/samples.svelte";

  let showAddDialog = $state(false);
  let selectedRaman = $state("");
  let selectedTarget = $state("");
  let isUpdating = $state(false);

  const availableRamanMolecules = ["DTNB", "MBA", "TFMBA"];
  const availableTargetMolecules = ["IgG", "BSA", "HER2", "EpCAM", "TROP2"];

  async function addPair() {
    if (!selectedRaman || !selectedTarget || !sampleStore.selectedSample) return;

    // Check for duplicate pairs
    const exists = sampleStore.selectedSample.moleculePairs.some(
      (p) => p.raman === selectedRaman && p.target === selectedTarget
    );
    if (exists) return;

    isUpdating = true;
    try {
      const newPairs = [
        ...sampleStore.selectedSample.moleculePairs,
        { raman: selectedRaman, target: selectedTarget },
      ];
      await sampleStore.updateSample(sampleStore.selectedSample.id, { moleculePairs: newPairs });

      // Reset selection
      selectedRaman = "";
      selectedTarget = "";
      showAddDialog = false;
    } finally {
      isUpdating = false;
    }
  }

  async function removePair(index: number) {
    if (!sampleStore.selectedSample) return;

    isUpdating = true;
    try {
      const newPairs = sampleStore.selectedSample.moleculePairs.filter((_, i) => i !== index);
      await sampleStore.updateSample(sampleStore.selectedSample.id, { moleculePairs: newPairs });
    } finally {
      isUpdating = false;
    }
  }

  function handleClickOutside(event: MouseEvent) {
    const target = event.target as HTMLElement;
    // Close only if clicking outside both the button and the dropdown
    if (!target.closest(".molecule-pair-editor")) {
      showAddDialog = false;
    }
  }

  $effect(() => {
    if (showAddDialog) {
      document.addEventListener("click", handleClickOutside);
      return () => document.removeEventListener("click", handleClickOutside);
    }
  });
</script>

<div class="molecule-pair-editor relative">
  <div class="flex flex-wrap gap-2">
    {#if sampleStore.selectedSample && sampleStore.selectedSample.moleculePairs.length > 0}
      {#each sampleStore.selectedSample.moleculePairs as pair, index}
        <span
          class="inline-flex items-center gap-1 rounded-md bg-pink-400/10 px-2 py-1 text-xs font-medium text-pink-400"
        >
          <span class="font-light opacity-75">{pair.raman}</span>
          <span class="mx-1 opacity-50">|</span>
          <span class="font-bold">{pair.target}</span>
          <button
            onclick={() => removePair(index)}
            disabled={isUpdating}
            class="ml-1 hover:text-pink-300 disabled:opacity-50"
            type="button"
            aria-label="Remove pair"
          >
            <XMarkIcon class="h-3 w-3" />
          </button>
        </span>
      {/each}
    {:else}
      <span
        class="inline-flex items-center rounded-md border border-dashed border-gray-700 px-2 py-1 text-xs text-gray-600"
      >
        No molecule pairs
      </span>
    {/if}

    <button
      onclick={(e) => {
        e.stopPropagation();
        showAddDialog = !showAddDialog;
      }}
      class="inline-flex items-center rounded-md border border-dashed border-gray-600 px-2 py-1 text-xs font-medium text-gray-500 hover:border-gray-500 hover:text-gray-400"
      type="button"
    >
      <PlusIcon class="h-3 w-3 mr-1" />
      Add Pair
    </button>
  </div>

  {#if showAddDialog}
    <div
      class="absolute left-0 top-10 z-20 mt-1 w-64 rounded-md bg-gray-800 shadow-lg ring-1 ring-black ring-opacity-5 p-4"
    >
      <div class="space-y-3">
        <div>
          <label for="raman-select" class="block text-xs font-medium text-gray-400 mb-1"
            >Raman Molecule</label
          >
          <select
            id="raman-select"
            bind:value={selectedRaman}
            class="w-full rounded-md bg-gray-700 border border-gray-600 px-3 py-1.5 text-sm text-gray-200 focus:border-blue-500 focus:outline-none"
          >
            <option value="">Select Raman...</option>
            {#each availableRamanMolecules as raman}
              <option value={raman}>{raman}</option>
            {/each}
          </select>
        </div>

        <div>
          <label for="target-select" class="block text-xs font-medium text-gray-400 mb-1"
            >Target Molecule</label
          >
          <select
            id="target-select"
            bind:value={selectedTarget}
            class="w-full rounded-md bg-gray-700 border border-gray-600 px-3 py-1.5 text-sm text-gray-200 focus:border-blue-500 focus:outline-none"
          >
            <option value="">Select Target...</option>
            {#each availableTargetMolecules as target}
              <option value={target}>{target}</option>
            {/each}
          </select>
        </div>

        {#if selectedRaman && selectedTarget}
          <div class="text-xs text-gray-500 text-center">
            <span class="opacity-75">{selectedRaman}</span>
            <span class="mx-1">|</span>
            <span class="font-semibold">{selectedTarget}</span>
          </div>
        {/if}

        <div class="flex gap-2 pt-2 border-t border-gray-700">
          <button
            onclick={addPair}
            disabled={!selectedRaman || !selectedTarget || isUpdating}
            class="flex-1 rounded-md bg-pink-600 px-3 py-1.5 text-sm font-medium text-white hover:bg-pink-700 disabled:opacity-50 disabled:cursor-not-allowed"
            type="button"
          >
            Add Pair
          </button>
          <button
            onclick={() => {
              showAddDialog = false;
              selectedRaman = "";
              selectedTarget = "";
            }}
            class="flex-1 rounded-md bg-gray-700 px-3 py-1.5 text-sm font-medium text-gray-300 hover:bg-gray-600"
            type="button"
          >
            Cancel
          </button>
        </div>
      </div>
    </div>
  {/if}
</div>
