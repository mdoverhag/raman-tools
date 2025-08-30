<script lang="ts">
  import { onMount } from "svelte";
  import * as d3 from "d3";

  interface Props {
    wavenumbers: number[];
    intensities: number[];
    title?: string;
  }

  let { wavenumbers, intensities, title = "Spectrum" }: Props = $props();

  let container: HTMLDivElement;

  function drawChart() {
    if (!container || !wavenumbers?.length || !intensities?.length) return;

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
    const data = wavenumbers.map((wn, i) => ({
      wavenumber: wn,
      intensity: intensities[i],
    }));

    // Set scales
    const xScale = d3
      .scaleLinear()
      .domain(d3.extent(wavenumbers) as [number, number])
      .range([0, width]);

    const yScale = d3
      .scaleLinear()
      .domain([0, d3.max(intensities) as number])
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

    // Add the line
    g.append("path")
      .datum(data)
      .attr("fill", "none")
      .attr("stroke", "#60a5fa")
      .attr("stroke-width", 1.5)
      .attr("d", line);

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

  // Redraw chart when props change
  $effect(() => {
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
