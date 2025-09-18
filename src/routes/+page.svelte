<script lang="ts">
  import { listen } from "@tauri-apps/api/event";
  import { invoke } from "@tauri-apps/api/core";
  import { onMount } from "svelte";
  import SpectrumChart from "$lib/SpectrumChart.svelte";
  import SampleSidebar from "$lib/components/SampleSidebar.svelte";
  import MoleculePairEditor from "$lib/components/MoleculePairEditor.svelte";
  import { sampleStore } from "$lib/stores/samples.svelte";

  let editingHeaderName = $state(false);
  let editingName = $state("");
  let headerNameInput = $state<HTMLInputElement | null>(null);

  function startEditingName() {
    if (!sampleStore.selectedSample) return;
    editingHeaderName = true;
    editingName = sampleStore.selectedSample.name;
  }

  async function saveNameEdit() {
    if (!sampleStore.selectedSample || !editingName.trim()) return;
    await sampleStore.updateSample(sampleStore.selectedSample.id, { name: editingName });
    editingHeaderName = false;
    editingName = "";
  }

  function cancelNameEdit() {
    editingHeaderName = false;
    editingName = "";
  }

  function handleNameKeydown(event: KeyboardEvent) {
    if (event.key === "Enter") {
      saveNameEdit();
    } else if (event.key === "Escape") {
      cancelNameEdit();
    }
  }

  let isDragging = $state(false);
  let isLoading = $state(false);
  let error = $state<string | null>(null);
  let importProgress = $state<{
    stage: string;
    current: number;
    total: number;
    filename: string;
  } | null>(null);

  type ImportResult = {
    success: boolean;
    count: number;
    sampleId: string;
  };

  async function handleFileDrop(paths: string[]) {
    // Check if we have a selected sample first
    if (!sampleStore.selectedSampleId) {
      error = "Please select a sample before importing spectra";
      return;
    }

    isLoading = true;
    error = null;
    importProgress = null;

    try {
      // Filter for .txt files only
      const txtFiles = paths.filter((path) => path.toLowerCase().endsWith(".txt"));

      if (txtFiles.length === 0) {
        error = "No .txt files found in dropped files";
        return;
      }

      // Get the selected sample ID
      const sampleId = sampleStore.selectedSampleId;

      // Call Rust command to parse the files and apply baseline correction
      const result = await invoke<ImportResult>("parse_spectrum_files", {
        filepaths: txtFiles,
        sampleId: sampleId,
      });

      if (!result.success || result.count === 0) {
        error = "No valid spectrum data found in any of the files";
      }
      // No need to reload - the store updates automatically via events
    } catch (e) {
      console.error("Error parsing files:", e);
      error = e instanceof Error ? e.message : String(e);
    } finally {
      isLoading = false;
      importProgress = null;
    }
  }

  onMount(() => {
    // Listen for file drop events from Tauri v2
    const unlisten = listen<{ paths: string[] }>("tauri://drag-drop", (event) => {
      isDragging = false;
      if (event.payload?.paths && sampleStore.selectedSampleId) {
        handleFileDrop(event.payload.paths);
      }
    });

    // Listen for drag over events
    const unlistenHover = listen("tauri://drag-over", () => {
      // Only show dragging state if a sample is selected
      if (sampleStore.selectedSampleId) {
        isDragging = true;
      }
    });

    // Listen for drag leave
    const unlistenCancelled = listen("tauri://drag-leave", () => {
      isDragging = false;
    });

    // Listen for import progress events
    const unlistenProgress = listen<any>("import:progress", (event) => {
      if (event.payload) {
        importProgress = {
          stage: event.payload.stage,
          current: event.payload.current,
          total: event.payload.total,
          filename: event.payload.filename,
        };
      }
    });

    // Note: Spectrum ready/updated events are no longer needed
    // since we fetch spectra from the backend when sample changes

    // Listen for import error events
    const unlistenError = listen<any>("import:error", (event) => {
      if (event.payload) {
        console.error(`Import error for ${event.payload.filename}: ${event.payload.error}`);
        // Could accumulate errors or show them in UI
      }
    });

    // Cleanup listeners
    return () => {
      unlisten.then((fn) => fn());
      unlistenHover.then((fn) => fn());
      unlistenCancelled.then((fn) => fn());
      unlistenProgress.then((fn) => fn());
      unlistenError.then((fn) => fn());
    };
  });
</script>

<div class="flex h-screen bg-gray-900">
  <!-- Sample Sidebar -->
  <SampleSidebar />

  <!-- Main Content -->
  <main class="flex-1 overflow-auto text-gray-100">
    {#if !sampleStore.selectedSample}
      <!-- No Sample Selected -->
      <div class="flex h-full items-center justify-center">
        <div class="text-center">
          <svg
            class="mx-auto h-16 w-16 text-gray-600 mb-4"
            fill="none"
            stroke="currentColor"
            viewBox="0 0 24 24"
          >
            <path
              stroke-linecap="round"
              stroke-linejoin="round"
              stroke-width="1.5"
              d="M19.428 15.428a2 2 0 00-1.022-.547l-2.387-.477a6 6 0 00-3.86.517l-.318.158a6 6 0 01-3.86.517L6.05 15.21a2 2 0 00-1.806.547M8 4h8l-1 1v5.172a2 2 0 00.586 1.414l5 5c1.26 1.26.367 3.414-1.415 3.414H4.828c-1.782 0-2.674-2.154-1.414-3.414l5-5A2 2 0 009 10.172V5L8 4z"
            />
          </svg>
          <h2 class="text-xl font-medium text-gray-400 mb-2">No Sample Selected</h2>
          <p class="text-gray-500">Choose a sample from the sidebar to get started</p>
        </div>
      </div>
    {:else}
      <div class="p-6">
        <!-- Breadcrumb and Header -->
        <div class="mb-6">
          <nav aria-label="Breadcrumb" class="flex">
            <ol role="list" class="flex items-center space-x-4">
              <li>
                <div class="flex">
                  <span class="text-sm font-medium text-gray-500">Samples</span>
                </div>
              </li>
              <li>
                <div class="flex items-center">
                  <svg
                    viewBox="0 0 20 20"
                    fill="currentColor"
                    aria-hidden="true"
                    class="size-5 shrink-0 text-gray-500"
                  >
                    <path
                      d="M8.22 5.22a.75.75 0 0 1 1.06 0l4.25 4.25a.75.75 0 0 1 0 1.06l-4.25 4.25a.75.75 0 0 1-1.06-1.06L11.94 10 8.22 6.28a.75.75 0 0 1 0-1.06Z"
                      clip-rule="evenodd"
                      fill-rule="evenodd"
                    />
                  </svg>
                  <span class="ml-4 text-sm font-medium text-gray-400" aria-current="page"
                    >{sampleStore.selectedSample.name}</span
                  >
                </div>
              </li>
            </ol>
          </nav>
          <div class="mt-2">
            {#if editingHeaderName}
              <input
                bind:this={headerNameInput}
                bind:value={editingName}
                onkeydown={handleNameKeydown}
                onblur={saveNameEdit}
                class="text-2xl font-bold bg-gray-800 text-white border border-gray-600 rounded px-2 py-1 focus:border-blue-500 focus:outline-none"
              />
            {:else}
              <h2 class="text-2xl font-bold text-white cursor-text" ondblclick={startEditingName}>
                {sampleStore.selectedSample.name}
              </h2>
            {/if}
            <div class="mt-3">
              <MoleculePairEditor />
            </div>
          </div>
        </div>

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
            {#if isLoading && importProgress}
              <div class="space-y-2">
                <p class="text-gray-400">
                  {importProgress.stage === "parsing"
                    ? `Parsing files... ${importProgress.current}/${importProgress.total}`
                    : importProgress.stage === "preparing"
                      ? "Initializing baseline correction..."
                      : importProgress.stage === "baseline"
                        ? `Applying baseline correction... ${importProgress.current}/${importProgress.total}`
                        : importProgress.stage === "averaging"
                          ? "Calculating average spectrum..."
                          : `Processing... ${importProgress.current}/${importProgress.total}`}
                </p>
                <p class="text-sm text-gray-500 h-5">
                  {importProgress.filename || "\u00A0"}
                </p>
                <div class="w-64 mx-auto bg-gray-700 rounded-full h-2">
                  <div
                    class="bg-blue-500 h-2 rounded-full"
                    style:width="{(importProgress.current / importProgress.total) * 100}%"
                  ></div>
                </div>
              </div>
            {:else if isLoading}
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
        {#if sampleStore.spectra.length > 0}
          <div class="grid grid-cols-1 lg:grid-cols-[350px_1fr] gap-6 mt-8">
            <!-- File List -->
            <div class="bg-gray-800 rounded-lg border border-gray-700 overflow-hidden">
              <div class="p-4 border-b border-gray-700 bg-gray-800/50">
                <h2 class="text-lg font-semibold text-gray-100">
                  Loaded Spectra ({sampleStore.spectra.length} files)
                </h2>
              </div>
              <div class="max-h-[600px] overflow-y-auto">
                <ul class="divide-y divide-gray-700">
                  {#each sampleStore.spectra as spectrum}
                    <li>
                      <button
                        class="w-full text-left px-4 py-3 hover:bg-gray-700/50 transition-colors {sampleStore.selectedSpectrumId ===
                        spectrum.id
                          ? 'bg-blue-900/30 border-l-2 border-blue-500'
                          : ''}"
                        onclick={() => sampleStore.selectSpectrum(spectrum.id)}
                      >
                        <div class="font-medium text-gray-200">
                          {spectrum.filename}
                        </div>
                        <div class="text-sm text-gray-500">
                          {spectrum.intensities.length} data points
                        </div>
                      </button>
                    </li>
                  {/each}
                </ul>
              </div>
            </div>

            <!-- Chart Section -->
            {#if sampleStore.selectedSpectrum}
              <div class="bg-gray-800 rounded-lg border border-gray-700 p-6">
                <SpectrumChart
                  spectrum={sampleStore.selectedSpectrum}
                  title={sampleStore.selectedSpectrum.filename}
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
    {/if}
  </main>
</div>
