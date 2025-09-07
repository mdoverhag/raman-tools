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

class SampleStore {
  samples = $state<Sample[]>([]);
  selectedSampleId = $state<string | null>(null);
  loading = $state(false);
  error = $state<string | null>(null);

  get selectedSample() {
    return this.samples.find((s) => s.id === this.selectedSampleId) || null;
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

  selectSample(id: string | null) {
    this.selectedSampleId = id;
  }
}

export const sampleStore = new SampleStore();
