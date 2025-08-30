<script lang="ts">
  import { listen } from "@tauri-apps/api/event";
  import { onMount } from "svelte";

  let isDragging = $state(false);

  onMount(() => {
    // Listen for file drop events from Tauri v2
    const unlisten = listen<{ paths: string[] }>("tauri://drag-drop", (event) => {
      isDragging = false;
      console.log("Files dropped:", event.payload);
      // event.payload.paths is an array of file paths
      if (event.payload?.paths) {
        event.payload.paths.forEach((path) => {
          console.log("File path:", path);
        });
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
    <p>
      {#if isDragging}
        Drop files here...
      {:else}
        Drag and drop spectrum files (.txt) here
      {/if}
    </p>
  </div>
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
</style>
