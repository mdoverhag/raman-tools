<script lang="ts">
  import { onMount } from "svelte";
  import { invoke } from "@tauri-apps/api/core";
  import { listen } from "@tauri-apps/api/event";

  interface UvStatus {
    installed: boolean;
    version: string | null;
    path: string | null;
  }

  interface PythonStatus {
    installed: boolean;
    python_version: string | null;
    numpy_version: string | null;
    scipy_version: string | null;
  }

  interface Progress {
    stage: string;
    message: string;
    percentage: number | null;
  }

  let showModal = $state(false);
  let currentStep = $state<"checking" | "uv" | "python" | "complete">("checking");
  let isWorking = $state(false);
  let progress = $state<Progress | null>(null);
  let error = $state<string | null>(null);

  onMount(async () => {
    // Check what needs to be installed
    const uvStatus = await invoke<UvStatus>("check_uv_status");
    const pythonStatus = await invoke<PythonStatus>("check_python_status");

    if (!uvStatus.installed) {
      showModal = true;
      currentStep = "uv";
    } else if (
      !pythonStatus.installed ||
      !pythonStatus.numpy_version ||
      !pythonStatus.scipy_version
    ) {
      showModal = true;
      currentStep = "python";
    }

    // Listen for progress events
    const unlistenUv = await listen<Progress>("uv-download-progress", (event) => {
      progress = event.payload;

      if (event.payload.stage === "complete") {
        // uv installed, move to Python setup
        setTimeout(() => {
          currentStep = "python";
          progress = null;
          startPythonSetup();
        }, 1000);
      }
    });

    const unlistenPython = await listen<Progress>("python-setup-progress", (event) => {
      progress = event.payload;

      if (event.payload.stage === "complete") {
        // Everything complete
        setTimeout(() => {
          showModal = false;
          currentStep = "complete";
          progress = null;
        }, 1500);
      }
    });

    return () => {
      unlistenUv();
      unlistenPython();
    };
  });

  async function startUvDownload() {
    isWorking = true;
    error = null;

    try {
      await invoke("download_uv");
    } catch (e) {
      error = String(e);
      isWorking = false;
    }
  }

  async function startPythonSetup() {
    isWorking = true;
    error = null;

    try {
      await invoke("setup_python_env");
    } catch (e) {
      error = String(e);
      isWorking = false;
    }
  }

  async function startSetup() {
    if (currentStep === "uv") {
      await startUvDownload();
    } else if (currentStep === "python") {
      await startPythonSetup();
    }
  }

  function getStepTitle() {
    switch (currentStep) {
      case "checking":
        return "Checking Python Runtime...";
      case "uv":
        return "Python Package Manager Setup";
      case "python":
        return "Python Environment Setup";
      case "complete":
        return "Setup Complete";
      default:
        return "Setup";
    }
  }

  function getStepDescription() {
    switch (currentStep) {
      case "uv":
        return "The app needs to download the Python package manager (uv) to manage Python and packages.";
      case "python":
        return "Setting up Python environment with scientific packages for baseline correction.";
      default:
        return "";
    }
  }
</script>

{#if showModal}
  <div class="fixed inset-0 bg-gray-600 bg-opacity-50 overflow-y-auto h-full w-full z-50">
    <div class="relative top-20 mx-auto p-5 border w-[480px] shadow-lg rounded-md bg-white">
      <div class="mt-3">
        <div class="flex items-center justify-center w-12 h-12 mx-auto bg-blue-100 rounded-full">
          <svg class="w-6 h-6 text-blue-600" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path
              stroke-linecap="round"
              stroke-linejoin="round"
              stroke-width="2"
              d="M9.75 17L9 20l-1 1h8l-1-1-.75-3M3 13h18M5 17h14a2 2 0 002-2V5a2 2 0 00-2-2H5a2 2 0 00-2 2v10a2 2 0 002 2z"
            />
          </svg>
        </div>

        <h3 class="text-lg leading-6 font-medium text-gray-900 text-center mt-4">
          {getStepTitle()}
        </h3>

        <div class="mt-2 px-7 py-3">
          {#if !isWorking && !error && getStepDescription()}
            <p class="text-sm text-gray-500 text-center">
              {getStepDescription()}
            </p>
            {#if currentStep === "uv"}
              <p class="text-xs text-gray-400 text-center mt-2">
                This is a one-time download (~30MB).
              </p>
            {:else if currentStep === "python"}
              <p class="text-xs text-gray-400 text-center mt-2">
                This will download Python 3.13 and scientific packages (~150MB).
              </p>
            {/if}
          {/if}

          {#if isWorking && progress}
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
              {:else}
                <div class="flex justify-center mt-3">
                  <div class="animate-spin rounded-full h-8 w-8 border-b-2 border-blue-600"></div>
                </div>
              {/if}
            </div>
          {/if}

          {#if error}
            <div class="mt-4 p-3 bg-red-50 border border-red-200 rounded">
              <p class="text-sm text-red-600">
                Setup failed: {error}
              </p>
            </div>
          {/if}
        </div>

        <div class="items-center px-4 py-3">
          {#if !isWorking && currentStep !== "checking"}
            <button
              onclick={startSetup}
              class="w-full px-4 py-2 bg-blue-600 text-white text-base font-medium rounded-md hover:bg-blue-700 focus:outline-none focus:ring-2 focus:ring-blue-300"
            >
              {#if error}
                Retry Setup
              {:else if currentStep === "uv"}
                Download Package Manager
              {:else if currentStep === "python"}
                Setup Python Environment
              {:else}
                Continue
              {/if}
            </button>
          {/if}
        </div>
      </div>
    </div>
  </div>
{/if}
