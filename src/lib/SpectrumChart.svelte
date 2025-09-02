<script lang="ts">
  import { onMount } from "svelte";
  import * as d3 from "d3";

  interface Spectrum {
    id: string;
    filename: string;
    wavenumber_start: number;
    wavenumber_end: number;
    wavenumber_step: number;
    intensities: number[];
    baseline?: number[] | null;
    corrected?: number[] | null;
  }

  interface Props {
    spectrum: Spectrum;
    title?: string;
  }

  let { spectrum, title = "Spectrum" }: Props = $props();

  let container: HTMLDivElement;

  // Helper function to generate wavenumbers array from start/end/step
  function generateWavenumbers(start: number, end: number, step: number): number[] {
    const count = Math.floor((end - start) / step) + 1;
    return Array.from({ length: count }, (_, i) => start + i * step);
  }

  function drawChart() {
    if (!container || !spectrum) return;

    // Generate wavenumbers array from spectrum parameters
    const wavenumbers = generateWavenumbers(
      spectrum.wavenumber_start,
      spectrum.wavenumber_end,
      spectrum.wavenumber_step
    );
    const intensities = spectrum.intensities;
    const baseline = spectrum.baseline;
    const corrected = spectrum.corrected;

    if (!wavenumbers?.length || !intensities?.length) return;

    // Set dimensions and margins
    const margin = { top: 40, right: 30, bottom: 50, left: 70 };
    const width = container.clientWidth - margin.left - margin.right;
    const height = 400 - margin.top - margin.bottom;

    // Clear any existing chart
    d3.select(container).selectAll("*").remove();

    // Create SVG
    const svg = d3
      .select(container)
      .append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom);

    const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);

    // Prepare data
    const originalData = wavenumbers.map((wn, i) => ({
      wavenumber: wn,
      intensity: intensities[i],
    }));

    const baselineData = baseline
      ? wavenumbers.map((wn, i) => ({
          wavenumber: wn,
          intensity: baseline[i],
        }))
      : null;

    const correctedData = corrected
      ? wavenumbers.map((wn, i) => ({
          wavenumber: wn,
          intensity: corrected[i],
        }))
      : null;

    // Set scales
    const xScale = d3
      .scaleLinear()
      .domain(d3.extent(wavenumbers) as [number, number])
      .range([0, width]);

    // Find the overall min and max for y-scale
    let allIntensities = [...intensities];
    if (baseline) allIntensities = allIntensities.concat(baseline);
    if (corrected) allIntensities = allIntensities.concat(corrected);

    const yScale = d3
      .scaleLinear()
      .domain([d3.min(allIntensities) as number, d3.max(allIntensities) as number])
      .range([height, 0]);

    // Create line generator
    const line = d3
      .line<{ wavenumber: number; intensity: number }>()
      .x((d) => xScale(d.wavenumber))
      .y((d) => yScale(d.intensity));

    // Add X axis
    const xAxis = g
      .append("g")
      .attr("transform", `translate(0,${height})`)
      .call(d3.axisBottom(xScale));

    xAxis.selectAll("text").style("fill", "#9ca3af");
    xAxis.selectAll("line").style("stroke", "#4b5563");
    xAxis.select(".domain").style("stroke", "#4b5563");

    xAxis
      .append("text")
      .attr("x", width / 2)
      .attr("y", 40)
      .attr("fill", "#9ca3af")
      .style("text-anchor", "middle")
      .text("Wavenumber (cm⁻¹)");

    // Add Y axis
    const yAxis = g.append("g").call(d3.axisLeft(yScale));

    yAxis.selectAll("text").style("fill", "#9ca3af");
    yAxis.selectAll("line").style("stroke", "#4b5563");
    yAxis.select(".domain").style("stroke", "#4b5563");

    yAxis
      .append("text")
      .attr("transform", "rotate(-90)")
      .attr("y", -50)
      .attr("x", -height / 2)
      .attr("fill", "#9ca3af")
      .style("text-anchor", "middle")
      .text("Intensity");

    // Add grid lines
    g.append("g")
      .attr("class", "grid")
      .attr("transform", `translate(0,${height})`)
      .call(
        d3
          .axisBottom(xScale)
          .tickSize(-height)
          .tickFormat(() => "")
      )
      .style("stroke-dasharray", "3,3")
      .style("opacity", 0.3)
      .selectAll("line")
      .style("stroke", "#374151");

    g.append("g")
      .attr("class", "grid")
      .call(
        d3
          .axisLeft(yScale)
          .tickSize(-width)
          .tickFormat(() => "")
      )
      .style("stroke-dasharray", "3,3")
      .style("opacity", 0.3)
      .selectAll("line")
      .style("stroke", "#374151");

    // Add the baseline line (if available)
    if (baselineData) {
      g.append("path")
        .datum(baselineData)
        .attr("fill", "none")
        .attr("stroke", "#ef4444") // red for baseline
        .attr("stroke-width", 1.5)
        .attr("stroke-dasharray", "5,5") // dashed line
        .attr("d", line)
        .attr("opacity", 0.8);
    }

    // Add the original spectrum line
    g.append("path")
      .datum(originalData)
      .attr("fill", "none")
      .attr("stroke", "#60a5fa") // blue for original
      .attr("stroke-width", 1.5)
      .attr("d", line)
      .attr("opacity", correctedData ? 0.5 : 1); // Make semi-transparent if corrected exists

    // Add the corrected spectrum line (if available)
    if (correctedData) {
      g.append("path")
        .datum(correctedData)
        .attr("fill", "none")
        .attr("stroke", "#10b981") // green for corrected
        .attr("stroke-width", 2)
        .attr("d", line);
    }

    // Add legend
    const legendData: Array<{
      label: string;
      color: string;
      dasharray: string | null;
      opacity: number;
    }> = [
      { label: "Original", color: "#60a5fa", dasharray: null, opacity: correctedData ? 0.5 : 1 },
    ];
    if (baselineData) {
      legendData.push({ label: "Baseline", color: "#ef4444", dasharray: "5,5", opacity: 0.8 });
    }
    if (correctedData) {
      legendData.push({ label: "Corrected", color: "#10b981", dasharray: null, opacity: 1 });
    }

    const legend = svg.append("g").attr("transform", `translate(${width - 100}, ${margin.top})`);

    legendData.forEach((item, i) => {
      const legendRow = legend.append("g").attr("transform", `translate(0, ${i * 20})`);

      legendRow
        .append("line")
        .attr("x1", 0)
        .attr("x2", 20)
        .attr("y1", 0)
        .attr("y2", 0)
        .attr("stroke", item.color)
        .attr("stroke-width", 2)
        .attr("stroke-dasharray", item.dasharray)
        .attr("opacity", item.opacity);

      legendRow
        .append("text")
        .attr("x", 25)
        .attr("y", 4)
        .style("font-size", "12px")
        .style("fill", "#9ca3af")
        .text(item.label);
    });

    // Add title
    svg
      .append("text")
      .attr("x", (width + margin.left + margin.right) / 2)
      .attr("y", margin.top / 2 + 5)
      .attr("text-anchor", "middle")
      .style("font-size", "16px")
      .style("font-weight", "bold")
      .style("fill", "#e5e7eb")
      .text(title);
  }

  // Draw chart when component mounts
  onMount(() => {
    drawChart();
  });

  // Redraw chart when spectrum changes
  $effect(() => {
    spectrum; // Track spectrum changes
    drawChart();
  });
</script>

<div bind:this={container} class="chart-container"></div>

<style>
  .chart-container {
    width: 100%;
    min-height: 400px;
  }
</style>
