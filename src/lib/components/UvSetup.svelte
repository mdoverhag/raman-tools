<script lang="ts">
  import { onMount } from "svelte";
  import { invoke } from "@tauri-apps/api/core";
  import { listen } from "@tauri-apps/api/event";

  interface UvStatus {
    installed: boolean;
    version: string | null;
    path: string | null;
  }

  interface DownloadProgress {
    stage: string;
    message: string;
    percentage: number | null;
  }

  let showModal = $state(false);
  let isDownloading = $state(false);
  let progress = $state<DownloadProgress | null>(null);
  let error = $state<string | null>(null);
  let uvStatus = $state<UvStatus | null>(null);

  onMount(async () => {
    // Check if uv is installed
    uvStatus = await invoke<UvStatus>("check_uv_status");

    if (!uvStatus.installed) {
      // Show setup modal
      showModal = true;
    }

    // Listen for download progress
    const unlisten = await listen<DownloadProgress>("uv-download-progress", (event) => {
      progress = event.payload;

      if (event.payload.stage === "complete") {
        // Download complete
        setTimeout(async () => {
          showModal = false;
          isDownloading = false;
          progress = null;
          // Re-check status
          uvStatus = await invoke<UvStatus>("check_uv_status");
        }, 1500);
      }
    });

    return () => {
      unlisten();
    };
  });

  async function startDownload() {
    isDownloading = true;
    error = null;

    try {
      await invoke("download_uv");
    } catch (e) {
      error = String(e);
      isDownloading = false;
    }
  }
</script>

{#if showModal}
  <div class="fixed inset-0 bg-gray-600 bg-opacity-50 overflow-y-auto h-full w-full z-50">
    <div class="relative top-20 mx-auto p-5 border w-96 shadow-lg rounded-md bg-white">
      <div class="mt-3">
        <div class="flex items-center justify-center w-12 h-12 mx-auto bg-blue-100 rounded-full">
          <svg class="w-6 h-6 text-blue-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path
              stroke-linecap="round"
              stroke-linejoin="round"
              stroke-width="2"
              d="M7 16a4 4 0 01-.88-7.903A5 5 0 1115.9 6L16 6a5 5 0 011 9.9M9 19l3 3m0 0l3-3m-3 3V10"
            />
          </svg>
        </div>

        <h3 class="text-lg leading-6 font-medium text-gray-900 text-center mt-4">
          Python Runtime Setup Required
        </h3>

        <div class="mt-2 px-7 py-3">
          {#if !isDownloading && !error}
            <p class="text-sm text-gray-500 text-center">
              The app needs to download the Python package manager (uv) to enable baseline
              correction.
            </p>
            <p class="text-xs text-gray-400 text-center mt-2">
              This is a one-time setup (~30MB download).
            </p>
          {/if}

          {#if isDownloading && progress}
            <div class="mt-4">
              <p class="text-sm text-gray-600 text-center mb-2">
                {progress.message}
              </p>
              {#if progress.percentage !== null}
                <div class="w-full bg-gray-200 rounded-full h-2.5">
                  <div
                    class="bg-blue-600 h-2.5 rounded-full transition-all duration-300"
                    style="width: {progress.percentage}%"
                  ></div>
                </div>
                <p class="text-xs text-gray-400 text-center mt-2">
                  {progress.percentage}%
                </p>
              {/if}
            </div>
          {/if}

          {#if error}
            <div class="mt-4 p-3 bg-red-50 border border-red-200 rounded">
              <p class="text-sm text-red-600">
                Download failed: {error}
              </p>
            </div>
          {/if}
        </div>

        <div class="items-center px-4 py-3">
          {#if !isDownloading}
            <button
              onclick={startDownload}
              class="w-full px-4 py-2 bg-blue-600 text-white text-base font-medium rounded-md hover:bg-blue-700 focus:outline-none focus:ring-2 focus:ring-blue-300"
            >
              {error ? "Retry Download" : "Download uv"}
            </button>
          {/if}
        </div>
      </div>
    </div>
  </div>
{/if}
