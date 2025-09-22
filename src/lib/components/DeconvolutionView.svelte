<script lang="ts">
  import type { Sample } from '$lib/types';

  interface Props {
    sample: Sample | null;
  }

  let { sample }: Props = $props();
</script>

<div class="bg-gray-800 rounded-lg border border-gray-700 p-6">
  <div class="flex flex-col items-center justify-center h-96">
    <svg
      class="h-16 w-16 text-purple-400 mb-4"
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
    <h2 class="text-2xl font-bold text-purple-300 mb-2">Deconvolution View</h2>
    <p class="text-gray-400 text-center max-w-md">
      This view will show the deconvoluted spectra for the multiplex sample
      containing {sample?.moleculePairs?.map((p) => p.raman).join(', ')}.
    </p>
    <p class="text-gray-500 text-sm mt-4">Feature coming soon</p>
  </div>
</div>