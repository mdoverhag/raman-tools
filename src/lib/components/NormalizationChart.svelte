<script lang="ts">
  import { onMount } from "svelte";
  import * as d3 from "d3";

  interface Props {
    normalizedData: any;
    sample: any;
    sampleStore: any;
  }

  let { normalizedData, sample, sampleStore }: Props = $props();

  let originalContainer = $state<HTMLDivElement>();
  let normalizedContainer = $state<HTMLDivElement>();

  // Calculate L2 norms for display
  function calculateL2Norm(spectrum: number[], startIdx: number, endIdx: number): number {
    const region = spectrum.slice(startIdx, endIdx + 1);
    const sumSquares = region.reduce((sum, val) => sum + val * val, 0);
    return Math.sqrt(sumSquares);
  }

  function getNormValues() {
    if (!normalizedData || !sample) return [];

    const wavenumbers = generateWavenumbers();
    const startIdx = wavenumbers.findIndex(w => w >= 1000);
    const endIdx = wavenumbers.findIndex(w => w > 1500) - 1;

    const results = [];

    // Add multiplex norm
    if (sample.averageCorrected) {
      const norm = calculateL2Norm(sample.averageCorrected, startIdx, endIdx);
      results.push({
        name: `${sample.name} (Multiplex)`,
        norm: norm,
        color: colors.multiplex
      });
    }

    // Add reference norms
    const referenceIds = normalizedData.referenceSampleIds || [];
    referenceIds.forEach((refId: string, index: number) => {
      const refSample = sampleStore.samples.find((s: any) => s.id === refId);
      if (refSample?.averageCorrected) {
        const norm = calculateL2Norm(refSample.averageCorrected, startIdx, endIdx);
        results.push({
          name: refSample.name,
          norm: norm,
          color: colors.references[index % colors.references.length]
        });
      }
    });

    return results;
  }

  let normValues = $derived(getNormValues());

  // Color palette for different spectra
  const colors = {
    multiplex: "rgb(168, 85, 247)", // Purple
    references: [
      "rgb(34, 197, 94)", // Green
      "rgb(59, 130, 246)", // Blue
      "rgb(251, 146, 60)", // Orange
      "rgb(244, 63, 94)", // Red
    ],
  };

  function generateWavenumbers(start = 200, end = 2000, points = 1801) {
    const step = (end - start) / (points - 1);
    return Array.from({ length: points }, (_, i) => start + i * step);
  }

  function drawOriginalChart() {
    if (!originalContainer || !sample) return;

    // Clear previous chart
    d3.select(originalContainer).selectAll("*").remove();

    const margin = { top: 40, right: 150, bottom: 50, left: 70 };
    const width = originalContainer.clientWidth - margin.left - margin.right;
    const height = 280 - margin.top - margin.bottom;

    const svg = d3
      .select(originalContainer)
      .append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom);

    const g = svg
      .append("g")
      .attr("transform", `translate(${margin.left},${margin.top})`);

    const wavenumbers = generateWavenumbers();

    // Collect all spectra data
    const spectraData = [];

    // Add multiplex spectrum
    if (sample.averageCorrected) {
      spectraData.push({
        name: `${sample.name} (Multiplex)`,
        values: sample.averageCorrected.map((y: number, i: number) => ({
          x: wavenumbers[i],
          y
        })),
        color: colors.multiplex,
        strokeWidth: 2
      });
    }

    // Add reference spectra
    const referenceIds = normalizedData.referenceSampleIds || [];
    referenceIds.forEach((refId: string, index: number) => {
      const refSample = sampleStore.samples.find((s: any) => s.id === refId);
      if (refSample?.averageCorrected) {
        spectraData.push({
          name: refSample.name,
          values: refSample.averageCorrected.map((y: number, i: number) => ({
            x: wavenumbers[i],
            y
          })),
          color: colors.references[index % colors.references.length],
          strokeWidth: 1.5
        });
      }
    });

    // Set up scales
    const xScale = d3
      .scaleLinear()
      .domain([200, 2000])
      .range([0, width]);

    const maxY = Number(d3.max(spectraData, (d: any) => d3.max(d.values, (v: any) => v.y))) || 100;
    const yScale = d3
      .scaleLinear()
      .domain([0, maxY * 1.1])
      .range([height, 0]);

    // Add axes
    g.append("g")
      .attr("transform", `translate(0,${height})`)
      .call(d3.axisBottom(xScale))
      .append("text")
      .attr("x", width / 2)
      .attr("y", 40)
      .attr("fill", "#9ca3af")
      .style("text-anchor", "middle")
      .text("Wavenumber (cm⁻¹)");

    g.append("g")
      .call(d3.axisLeft(yScale).tickFormat(d3.format(".2s")))
      .append("text")
      .attr("transform", "rotate(-90)")
      .attr("y", -50)
      .attr("x", -height / 2)
      .attr("fill", "#9ca3af")
      .style("text-anchor", "middle")
      .text("Intensity");

    // Add normalization range highlight
    g.append("rect")
      .attr("x", xScale(1000))
      .attr("y", 0)
      .attr("width", xScale(1500) - xScale(1000))
      .attr("height", height)
      .attr("fill", "rgba(168, 85, 247, 0.1)")
      .attr("stroke", "rgba(168, 85, 247, 0.3)")
      .attr("stroke-width", 1)
      .attr("stroke-dasharray", "5,5");

    // Add normalization range label
    g.append("text")
      .attr("x", xScale(1250))
      .attr("y", 15)
      .attr("fill", "rgba(168, 85, 247, 0.8)")
      .style("text-anchor", "middle")
      .style("font-size", "12px")
      .text("Normalization Range");

    // Draw lines
    const line = d3
      .line<{x: number, y: number}>()
      .x((d) => xScale(d.x))
      .y((d) => yScale(d.y))
      .curve(d3.curveMonotoneX);

    spectraData.forEach((spectrum) => {
      g.append("path")
        .datum(spectrum.values)
        .attr("fill", "none")
        .attr("stroke", spectrum.color)
        .attr("stroke-width", spectrum.strokeWidth)
        .attr("d", line);
    });

    // Add title
    svg
      .append("text")
      .attr("x", (width + margin.left + margin.right) / 2)
      .attr("y", 20)
      .attr("text-anchor", "middle")
      .style("font-size", "14px")
      .style("fill", "#e5e7eb")
      .text("Original Spectra (Average Corrected)");

    // Add legend
    const legend = svg
      .append("g")
      .attr("transform", `translate(${width + margin.left + 10}, ${margin.top})`);

    spectraData.forEach((spectrum, i) => {
      const legendItem = legend
        .append("g")
        .attr("transform", `translate(0, ${i * 20})`);

      legendItem
        .append("line")
        .attr("x1", 0)
        .attr("x2", 15)
        .attr("y1", 0)
        .attr("y2", 0)
        .attr("stroke", spectrum.color)
        .attr("stroke-width", spectrum.strokeWidth);

      legendItem
        .append("text")
        .attr("x", 20)
        .attr("y", 4)
        .style("font-size", "11px")
        .style("fill", "#e5e7eb")
        .text(spectrum.name);
    });
  }

  function drawNormalizedChart() {
    if (!normalizedContainer || !normalizedData) return;

    // Clear previous chart
    d3.select(normalizedContainer).selectAll("*").remove();

    const margin = { top: 40, right: 150, bottom: 50, left: 70 };
    const width = normalizedContainer.clientWidth - margin.left - margin.right;
    const height = 280 - margin.top - margin.bottom;

    const svg = d3
      .select(normalizedContainer)
      .append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom);

    const g = svg
      .append("g")
      .attr("transform", `translate(${margin.left},${margin.top})`);

    const wavenumbers = generateWavenumbers();
    const displayRange = [1000, 1500];

    // Filter wavenumbers and data to display range
    const startIdx = wavenumbers.findIndex(w => w >= displayRange[0]);
    const endIdx = wavenumbers.findIndex(w => w > displayRange[1]);
    const endIndex = endIdx === -1 ? wavenumbers.length - 1 : endIdx - 1;

    // Collect normalized spectra data
    const spectraData = [];

    // Add normalized multiplex
    if (normalizedData.normalizedMultiplex) {
      spectraData.push({
        name: `${sample.name} (Multiplex) - Normalized`,
        values: normalizedData.normalizedMultiplex
          .slice(startIdx, endIndex + 1)
          .map((y: number, i: number) => ({
            x: wavenumbers[startIdx + i],
            y
          })),
        color: colors.multiplex,
        strokeWidth: 2
      });
    }

    // Add normalized references
    if (normalizedData.normalizedReferences) {
      Object.entries(normalizedData.normalizedReferences).forEach(
        ([ramanName, spectrum]: [string, any], index) => {
          spectraData.push({
            name: `${ramanName} - Normalized`,
            values: spectrum
              .slice(startIdx, endIndex + 1)
              .map((y: number, i: number) => ({
                x: wavenumbers[startIdx + i],
                y
              })),
            color: colors.references[index % colors.references.length],
            strokeWidth: 1.5
          });
        }
      );
    }

    // Set up scales
    const xScale = d3
      .scaleLinear()
      .domain(displayRange)
      .range([0, width]);

    const maxY = Number(d3.max(spectraData, (d: any) => d3.max(d.values, (v: any) => v.y))) || 0.2;
    const yScale = d3
      .scaleLinear()
      .domain([0, maxY * 1.1])
      .range([height, 0]);

    // Add axes
    g.append("g")
      .attr("transform", `translate(0,${height})`)
      .call(d3.axisBottom(xScale))
      .append("text")
      .attr("x", width / 2)
      .attr("y", 40)
      .attr("fill", "#9ca3af")
      .style("text-anchor", "middle")
      .text("Wavenumber (cm⁻¹)");

    g.append("g")
      .call(d3.axisLeft(yScale).tickFormat(d3.format(".3f")))
      .append("text")
      .attr("transform", "rotate(-90)")
      .attr("y", -50)
      .attr("x", -height / 2)
      .attr("fill", "#9ca3af")
      .style("text-anchor", "middle")
      .text("Normalized Intensity");

    // Draw lines
    const line = d3
      .line<{x: number, y: number}>()
      .x((d) => xScale(d.x))
      .y((d) => yScale(d.y))
      .curve(d3.curveMonotoneX);

    spectraData.forEach((spectrum) => {
      g.append("path")
        .datum(spectrum.values)
        .attr("fill", "none")
        .attr("stroke", spectrum.color)
        .attr("stroke-width", spectrum.strokeWidth)
        .attr("d", line);
    });

    // Add title
    svg
      .append("text")
      .attr("x", (width + margin.left + margin.right) / 2)
      .attr("y", 20)
      .attr("text-anchor", "middle")
      .style("font-size", "14px")
      .style("fill", "#e5e7eb")
      .text("Normalized Spectra (L2 Normalization in 1000-1500 cm⁻¹)");

    // Add legend
    const legend = svg
      .append("g")
      .attr("transform", `translate(${width + margin.left + 10}, ${margin.top})`);

    spectraData.forEach((spectrum, i) => {
      const legendItem = legend
        .append("g")
        .attr("transform", `translate(0, ${i * 20})`);

      legendItem
        .append("line")
        .attr("x1", 0)
        .attr("x2", 15)
        .attr("y1", 0)
        .attr("y2", 0)
        .attr("stroke", spectrum.color)
        .attr("stroke-width", spectrum.strokeWidth);

      legendItem
        .append("text")
        .attr("x", 20)
        .attr("y", 4)
        .style("font-size", "11px")
        .style("fill", "#e5e7eb")
        .text(spectrum.name);
    });
  }

  // Redraw charts when data changes
  $effect(() => {
    if (normalizedData && sample && sampleStore) {
      drawOriginalChart();
      drawNormalizedChart();
    }
  });


  // Handle window resize
  onMount(() => {
    const handleResize = () => {
      if (normalizedData && sample && sampleStore) {
        drawOriginalChart();
        drawNormalizedChart();
      }
    };

    window.addEventListener("resize", handleResize);
    return () => window.removeEventListener("resize", handleResize);
  });
</script>

<div class="space-y-4">
  <!-- Controls -->
  <div class="flex justify-between items-center">
    <h3 class="text-lg font-semibold text-gray-200">Normalization Results</h3>
  </div>

  <!-- Charts -->
  <div class="grid grid-cols-1 gap-4">
    <!-- Original Spectra -->
    <div class="bg-gray-900/50 rounded-lg p-4 border border-gray-700">
      <div bind:this={originalContainer} class="w-full h-[280px]"></div>
    </div>

    <!-- Normalized Spectra -->
    <div class="bg-gray-900/50 rounded-lg p-4 border border-gray-700">
      <div bind:this={normalizedContainer} class="w-full h-[280px]"></div>
    </div>
  </div>

  <!-- L2 Norm Values Table -->
  <div class="bg-gray-800/50 rounded-lg p-4 border border-gray-700">
    <h4 class="text-sm font-semibold text-gray-300 mb-3">L2 Norm Values (1000-1500 cm⁻¹)</h4>
    <div class="space-y-2">
      {#each normValues as item}
        <div class="flex items-center justify-between p-2 bg-gray-900/50 rounded">
          <div class="flex items-center gap-2">
            <div
              class="w-3 h-3 rounded-full"
              style="background-color: {item.color}"
            ></div>
            <span class="text-sm text-gray-300">{item.name}</span>
          </div>
          <span class="text-sm font-mono text-gray-200">
            {item.norm.toFixed(2)}
          </span>
        </div>
      {/each}
    </div>

    <div class="mt-4 pt-4 border-t border-gray-700 grid grid-cols-2 gap-4 text-sm">
      <div>
        <span class="text-gray-400">Normalization Method:</span>
        <span class="text-gray-200 ml-2">L2 (Euclidean)</span>
      </div>
      <div>
        <span class="text-gray-400">Range:</span>
        <span class="text-gray-200 ml-2">1000-1500 cm⁻¹</span>
      </div>
    </div>
  </div>
</div>