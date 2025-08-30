<script lang="ts">
  import { listen } from "@tauri-apps/api/event";
  import { invoke } from "@tauri-apps/api/core";
  import { onMount } from "svelte";
  import SpectrumChart from "$lib/SpectrumChart.svelte";

  interface Spectrum {
    filename: string;
    filepath: string;
    wavenumbers: number[];
    intensities: number[];
  }

  let isDragging = $state(false);
  let spectra = $state<Spectrum[]>([]);
  let isLoading = $state(false);
  let error = $state<string | null>(null);
  let selectedSpectrum = $state<Spectrum | null>(null);

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
        spectra = parsedSpectra;
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

<main>
  <h1>Raman Spectrum Viewer</h1>

  <div class="drop-zone" class:dragging={isDragging}>
    {#if isLoading}
      <p>Processing files...</p>
    {:else if isDragging}
      <p>Drop files here...</p>
    {:else}
      <p>Drag and drop spectrum files (.txt) here</p>
    {/if}
  </div>

  {#if error}
    <div class="error">
      <p>Error: {error}</p>
    </div>
  {/if}

  {#if spectra.length > 0}
    <div class="content-container">
      <div class="file-list">
        <h2>Loaded Spectra ({spectra.length} files)</h2>
        <ul>
          {#each spectra as spectrum}
            <li>
              <button
                class="file-item"
                class:selected={selectedSpectrum === spectrum}
                onclick={() => (selectedSpectrum = spectrum)}
              >
                {spectrum.filename}
                <span class="data-points">({spectrum.wavenumbers.length} points)</span>
              </button>
            </li>
          {/each}
        </ul>
      </div>

      {#if selectedSpectrum}
        <div class="chart-section">
          <SpectrumChart
            wavenumbers={selectedSpectrum.wavenumbers}
            intensities={selectedSpectrum.intensities}
            title={selectedSpectrum.filename}
          />
        </div>
      {/if}
    </div>
  {/if}
</main>

<style>
  main {
    max-width: 1200px;
    margin: 0 auto;
    padding: 2rem;
    font-family:
      system-ui,
      -apple-system,
      sans-serif;
  }

  h1 {
    color: #333;
    margin-bottom: 2rem;
  }

  .drop-zone {
    border: 2px dashed #ccc;
    border-radius: 8px;
    padding: 4rem 2rem;
    text-align: center;
    transition: all 0.3s ease;
    background-color: #fafafa;
  }

  .drop-zone.dragging {
    border-color: #4a90e2;
    background-color: #e8f2ff;
  }

  .drop-zone p {
    color: #666;
    font-size: 1.1rem;
    margin: 0;
  }

  .drop-zone.dragging p {
    color: #4a90e2;
    font-weight: 500;
  }

  .error {
    background-color: #fee;
    border: 1px solid #fcc;
    border-radius: 4px;
    padding: 1rem;
    margin: 1rem 0;
  }

  .error p {
    color: #c00;
    margin: 0;
  }

  .content-container {
    display: grid;
    grid-template-columns: 350px 1fr;
    gap: 2rem;
    margin-top: 2rem;
  }

  .file-list {
    max-height: 600px;
    overflow-y: auto;
  }

  .file-list h2 {
    color: #333;
    margin-bottom: 1rem;
    position: sticky;
    top: 0;
    background: white;
    padding: 0.5rem 0;
  }

  .file-list ul {
    list-style: none;
    padding: 0;
    margin: 0;
  }

  .file-list li {
    margin-bottom: 0.5rem;
  }

  .file-item {
    width: 100%;
    text-align: left;
    padding: 0.75rem 1rem;
    background: white;
    border: 1px solid #ddd;
    border-radius: 4px;
    cursor: pointer;
    transition: all 0.2s;
    font-size: 1rem;
    font-family: inherit;
  }

  .file-item:hover {
    background: #f5f5f5;
    border-color: #4a90e2;
  }

  .file-item.selected {
    background: #e8f2ff;
    border-color: #4a90e2;
    font-weight: 500;
  }

  .data-points {
    color: #666;
    font-size: 0.9rem;
    margin-left: 0.5rem;
  }

  .chart-section {
    background: white;
    border: 1px solid #ddd;
    border-radius: 4px;
    padding: 1rem;
  }

  @media (max-width: 768px) {
    .content-container {
      grid-template-columns: 1fr;
    }
  }
</style>
