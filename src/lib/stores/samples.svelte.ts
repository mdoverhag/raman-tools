import { invoke } from "@tauri-apps/api/core";

export interface Sample {
  id: string;
  name: string;
  ramanMolecules: string[];
  targetMolecules: string[];
  spectrumIds: string[];
}

export interface UpdateSampleData {
  name?: string;
  ramanMolecules?: string[];
  targetMolecules?: string[];
}

// This interface matches what we get from the backend
export interface Spectrum {
  id: string;
  filename: string;
  wavenumber_start: number;
  wavenumber_end: number;
  wavenumber_step: number;
  intensities: number[];
  baseline?: number[] | null;
  corrected?: number[] | null;
}

class SampleStore {
  samples = $state<Sample[]>([]);
  selectedSampleId = $state<string | null>(null);
  loading = $state(false);
  error = $state<string | null>(null);

  // New state for spectra management
  spectra = $state<Spectrum[]>([]);
  selectedSpectrumId = $state<string | null>(null);
  spectraLoading = $state(false);

  get selectedSample() {
    return this.samples.find((s) => s.id === this.selectedSampleId) || null;
  }

  get selectedSpectrum() {
    return this.spectra.find((s) => s.id === this.selectedSpectrumId) || null;
  }

  constructor() {
    // Start with empty samples, will load from backend
    this.loadSamples();
  }

  async loadSamples() {
    this.loading = true;
    this.error = null;
    try {
      this.samples = await invoke<Sample[]>("list_samples");
    } catch (err) {
      this.error = err instanceof Error ? err.message : "Failed to load samples";
    } finally {
      this.loading = false;
    }
  }

  async createSample(name: string): Promise<Sample | null> {
    this.error = null;
    try {
      const sample = await invoke<Sample>("create_sample", { name });
      // Need to reassign array for Svelte 5 reactivity
      this.samples = [...this.samples, sample];

      // Clear spectra for new sample (we know it has none)
      this.spectra = [];
      this.selectedSpectrumId = null;

      return sample;
    } catch (err) {
      this.error = err instanceof Error ? err.message : "Failed to create sample";
      return null;
    }
  }

  async updateSample(id: string, updates: UpdateSampleData): Promise<Sample | null> {
    this.error = null;
    try {
      const updated = await invoke<Sample>("update_sample", { id, updates });
      // Update local state with reassignment for reactivity
      const index = this.samples.findIndex((s) => s.id === id);
      if (index !== -1) {
        this.samples[index] = updated;
        this.samples = [...this.samples];
      }
      return updated;
    } catch (err) {
      this.error = err instanceof Error ? err.message : "Failed to update sample";
      return null;
    }
  }

  async deleteSample(id: string): Promise<boolean> {
    this.error = null;
    try {
      await invoke("delete_sample", { id });
      // Remove from local state with reassignment for reactivity
      this.samples = this.samples.filter((s) => s.id !== id);
      if (this.selectedSampleId === id) {
        this.selectedSampleId = null;
      }
      return true;
    } catch (err) {
      this.error = err instanceof Error ? err.message : "Failed to delete sample";
      return false;
    }
  }

  async selectSample(id: string | null) {
    this.selectedSampleId = id;

    // Clear spectra when no sample selected
    if (!id) {
      this.spectra = [];
      this.selectedSpectrumId = null;
      return;
    }

    // Find the sample
    const sample = this.samples.find((s) => s.id === id);
    if (!sample) {
      this.spectra = [];
      this.selectedSpectrumId = null;
      return;
    }

    // If sample has no spectra, just clear
    if (sample.spectrumIds.length === 0) {
      this.spectra = [];
      this.selectedSpectrumId = null;
      return;
    }

    // Load spectra for this sample
    await this.loadSpectraForSample(sample);
  }

  // New method to load spectra
  private async loadSpectraForSample(sample: Sample) {
    this.spectraLoading = true;
    this.spectra = []; // Clear immediately
    this.selectedSpectrumId = null;

    try {
      const loadedSpectra: Spectrum[] = [];
      for (const id of sample.spectrumIds) {
        const spectrum = await invoke<Spectrum>("get_spectrum", { id });
        loadedSpectra.push(spectrum);
      }

      this.spectra = loadedSpectra;

      // Auto-select first spectrum if we have any
      if (loadedSpectra.length > 0) {
        this.selectedSpectrumId = loadedSpectra[0].id;
      }
    } catch (err) {
      console.error("Error loading spectra:", err);
      this.error = err instanceof Error ? err.message : "Failed to load spectra";
    } finally {
      this.spectraLoading = false;
    }
  }

  selectSpectrum(id: string | null) {
    this.selectedSpectrumId = id;
  }

  // Call this after importing new spectra files
  async reloadCurrentSampleSpectra() {
    if (!this.selectedSampleId) return;

    // Reload the samples list to get updated spectrumIds
    await this.loadSamples();

    // Find the current sample with updated data
    const sample = this.samples.find((s) => s.id === this.selectedSampleId);
    if (sample && sample.spectrumIds.length > 0) {
      await this.loadSpectraForSample(sample);
    }
  }
}

export const sampleStore = new SampleStore();
