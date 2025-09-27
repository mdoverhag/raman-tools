<script lang="ts">
  import { onMount } from "svelte";
  import * as d3 from "d3";

  interface Props {
    normalizedData: any;
    sample: any;
    sampleStore: any;
    deconvolutionResults?: any;
  }

  let { normalizedData, sample, sampleStore, deconvolutionResults }: Props = $props();

  let originalContainer = $state<HTMLDivElement>();
  let normalizedContainer = $state<HTMLDivElement>();
  let deconvolutionContainer = $state<HTMLDivElement>();
  let residualContainer = $state<HTMLDivElement>();

  // Calculate L2 norms for display
  function calculateL2Norm(spectrum: number[], startIdx: number, endIdx: number): number {
    const region = spectrum.slice(startIdx, endIdx + 1);
    const sumSquares = region.reduce((sum, val) => sum + val * val, 0);
    return Math.sqrt(sumSquares);
  }

  function getNormValues() {
    if (!normalizedData || !sample) return [];

    const wavenumbers = generateWavenumbers();
    const startIdx = wavenumbers.findIndex((w) => w >= 1000);
    const endIdx = wavenumbers.findIndex((w) => w > 1500) - 1;

    const results = [];

    // Add multiplex norm
    if (sample.averageCorrected) {
      const norm = calculateL2Norm(sample.averageCorrected, startIdx, endIdx);
      results.push({
        name: `${sample.name} (Multiplex)`,
        norm: norm,
        color: colors.multiplex,
      });
    }

    // Add reference norms
    const referenceIds = normalizedData.referenceSampleIds || [];
    referenceIds.forEach((refId: string) => {
      const refSample = sampleStore.samples.find((s: any) => s.id === refId);
      if (refSample?.averageCorrected) {
        const norm = calculateL2Norm(refSample.averageCorrected, startIdx, endIdx);

        // Get fixed color based on molecule
        let moleculeColor = colors.references[0]; // default
        if (refSample.moleculePairs && refSample.moleculePairs.length > 0) {
          const ramanMolecule = refSample.moleculePairs[0].raman;
          moleculeColor = moleculeColors[ramanMolecule] || colors.references[0];
        }

        results.push({
          name: refSample.name,
          norm: norm,
          color: moleculeColor,
        });
      }
    });

    return results;
  }

  let normValues = $derived(getNormValues());

  // Fixed color palette for Raman molecules
  const moleculeColors: Record<string, string> = {
    "MBA": "rgb(34, 197, 94)", // Green
    "DTNB": "rgb(59, 130, 246)", // Blue
    "TFMBA": "rgb(251, 146, 60)", // Orange
  };

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

  // Helper function to get color for a molecule
  function getMoleculeColor(name: string): string {
    // Extract molecule name from the string (e.g., "MBA" from "MBA (45.2%)")
    for (const molecule in moleculeColors) {
      if (name.includes(molecule)) {
        return moleculeColors[molecule];
      }
    }
    // Fallback to default colors if not found
    return colors.references[0];
  }

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

    const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);

    const wavenumbers = generateWavenumbers();

    // Collect all spectra data
    const spectraData = [];

    // Add multiplex spectrum
    if (sample.averageCorrected) {
      spectraData.push({
        name: `${sample.name} (Multiplex)`,
        values: sample.averageCorrected.map((y: number, i: number) => ({
          x: wavenumbers[i],
          y,
        })),
        color: colors.multiplex,
        strokeWidth: 2,
      });
    }

    // Add reference spectra
    const referenceIds = normalizedData.referenceSampleIds || [];
    referenceIds.forEach((refId: string) => {
      const refSample = sampleStore.samples.find((s: any) => s.id === refId);
      if (refSample?.averageCorrected) {
        // Get the Raman molecule name for color
        let moleculeColor = colors.references[0]; // default
        if (refSample.moleculePairs && refSample.moleculePairs.length > 0) {
          const ramanMolecule = refSample.moleculePairs[0].raman;
          moleculeColor = moleculeColors[ramanMolecule] || colors.references[0];
        }

        spectraData.push({
          name: refSample.name,
          values: refSample.averageCorrected.map((y: number, i: number) => ({
            x: wavenumbers[i],
            y,
          })),
          color: moleculeColor,
          strokeWidth: 1.5,
        });
      }
    });

    // Set up scales
    const xScale = d3.scaleLinear().domain([200, 2000]).range([0, width]);

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
      .line<{ x: number; y: number }>()
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
      const legendItem = legend.append("g").attr("transform", `translate(0, ${i * 20})`);

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

    const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);

    const wavenumbers = generateWavenumbers();
    const displayRange = [1000, 1500];

    // Filter wavenumbers and data to display range
    const startIdx = wavenumbers.findIndex((w) => w >= displayRange[0]);
    const endIdx = wavenumbers.findIndex((w) => w > displayRange[1]);
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
            y,
          })),
        color: colors.multiplex,
        strokeWidth: 2,
      });
    }

    // Add normalized references
    if (normalizedData.normalizedReferences) {
      // Define the desired order
      const moleculeOrder = ["DTNB", "MBA", "TFMBA"];

      // Sort entries by the defined order
      const sortedEntries = Object.entries(normalizedData.normalizedReferences).sort(
        ([nameA], [nameB]) => {
          const indexA = moleculeOrder.indexOf(nameA);
          const indexB = moleculeOrder.indexOf(nameB);
          // If not in list, put at end
          const orderA = indexA === -1 ? 999 : indexA;
          const orderB = indexB === -1 ? 999 : indexB;
          return orderA - orderB;
        }
      );

      sortedEntries.forEach(
        ([ramanName, spectrum]: [string, any]) => {
          // Use fixed color based on molecule name
          const moleculeColor = moleculeColors[ramanName] || colors.references[0];

          spectraData.push({
            name: `${ramanName} - Normalized`,
            values: spectrum.slice(startIdx, endIndex + 1).map((y: number, i: number) => ({
              x: wavenumbers[startIdx + i],
              y,
            })),
            color: moleculeColor,
            strokeWidth: 1.5,
          });
        }
      );
    }

    // Set up scales
    const xScale = d3.scaleLinear().domain(displayRange).range([0, width]);

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
      .line<{ x: number; y: number }>()
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
      const legendItem = legend.append("g").attr("transform", `translate(0, ${i * 20})`);

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

  function drawDeconvolutionChart() {
    if (!deconvolutionContainer || !deconvolutionResults || !normalizedData) return;

    // Clear previous chart
    d3.select(deconvolutionContainer).selectAll("*").remove();

    const margin = { top: 40, right: 150, bottom: 50, left: 70 };
    const width = deconvolutionContainer.clientWidth - margin.left - margin.right;
    const height = 280 - margin.top - margin.bottom;

    const svg = d3
      .select(deconvolutionContainer)
      .append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom);

    const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);

    const wavenumbers = generateWavenumbers();
    const displayRange = [1000, 1500];

    // Filter to display range
    const startIdx = wavenumbers.findIndex((w) => w >= displayRange[0]);
    const endIdx = wavenumbers.findIndex((w) => w > displayRange[1]);
    const endIndex = endIdx === -1 ? wavenumbers.length - 1 : endIdx - 1;

    // Prepare data
    const spectraData = [];

    // Original multiplex (normalized)
    if (normalizedData.normalizedMultiplex) {
      spectraData.push({
        name: "Multiplex (measured)",
        values: normalizedData.normalizedMultiplex
          .slice(startIdx, endIndex + 1)
          .map((y: number, i: number) => ({
            x: wavenumbers[startIdx + i],
            y,
          })),
        color: colors.multiplex,
        strokeWidth: 2.5,
        strokeDasharray: "5,5",
      });
    }

    // Individual component contributions (scaled by coefficients)
    if (normalizedData.normalizedReferences && deconvolutionResults.coefficients) {
      // Define the desired order
      const moleculeOrder = ["DTNB", "MBA", "TFMBA"];

      // Sort entries by the defined order
      const sortedEntries = Object.entries(deconvolutionResults.coefficients).sort(
        ([nameA], [nameB]) => {
          const indexA = moleculeOrder.indexOf(nameA);
          const indexB = moleculeOrder.indexOf(nameB);
          // If not in list, put at end
          const orderA = indexA === -1 ? 999 : indexA;
          const orderB = indexB === -1 ? 999 : indexB;
          return orderA - orderB;
        }
      );

      sortedEntries.forEach(
        ([name, coefficient]: [string, any]) => {
          const refSpectrum = (normalizedData.normalizedReferences as any)[name];
          if (refSpectrum && coefficient > 0) {
            // Use fixed color based on molecule name
            const moleculeColor = moleculeColors[name] || colors.references[0];

            spectraData.push({
              name: `${name} (${deconvolutionResults.contributions[name].toFixed(1)}%)`,
              values: refSpectrum
                .slice(startIdx, endIndex + 1)
                .map((y: number, i: number) => ({
                  x: wavenumbers[startIdx + i],
                  y: y * coefficient, // Scale by NNLS coefficient
                })),
              color: moleculeColor,
              strokeWidth: 2,
              strokeDasharray: null,
            });
          }
        }
      );
    }

    // Reconstructed spectrum (sum of components)
    if (deconvolutionResults.reconstructedSpectrum) {
      spectraData.push({
        name: "Reconstruction (fit)",
        values: deconvolutionResults.reconstructedSpectrum
          .slice(startIdx, endIndex + 1)
          .map((y: number, i: number) => ({
            x: wavenumbers[startIdx + i],
            y,
          })),
        color: "rgb(239, 68, 68)", // Red
        strokeWidth: 2,
        strokeDasharray: "5,5",
      });
    }

    // Set up scales
    const xScale = d3.scaleLinear().domain(displayRange).range([0, width]);

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
      .line<{ x: number; y: number }>()
      .x((d) => xScale(d.x))
      .y((d) => yScale(d.y))
      .curve(d3.curveMonotoneX);

    spectraData.forEach((spectrum) => {
      const path = g
        .append("path")
        .datum(spectrum.values)
        .attr("fill", "none")
        .attr("stroke", spectrum.color)
        .attr("stroke-width", spectrum.strokeWidth)
        .attr("d", line);

      if (spectrum.strokeDasharray) {
        path.attr("stroke-dasharray", spectrum.strokeDasharray);
      }
    });

    // Add title
    svg
      .append("text")
      .attr("x", (width + margin.left + margin.right) / 2)
      .attr("y", 20)
      .attr("text-anchor", "middle")
      .style("font-size", "14px")
      .style("fill", "#e5e7eb")
      .text("NNLS Deconvolution Components");

    // Add legend
    const legend = svg
      .append("g")
      .attr("transform", `translate(${width + margin.left + 10}, ${margin.top})`);

    spectraData.forEach((spectrum, i) => {
      const legendItem = legend.append("g").attr("transform", `translate(0, ${i * 20})`);

      const line = legendItem
        .append("line")
        .attr("x1", 0)
        .attr("x2", 15)
        .attr("y1", 0)
        .attr("y2", 0)
        .attr("stroke", spectrum.color)
        .attr("stroke-width", spectrum.strokeWidth);

      if (spectrum.strokeDasharray) {
        line.attr("stroke-dasharray", spectrum.strokeDasharray);
      }

      legendItem
        .append("text")
        .attr("x", 20)
        .attr("y", 4)
        .style("font-size", "11px")
        .style("fill", "#e5e7eb")
        .text(spectrum.name);
    });

    // Add R² annotation
    if (deconvolutionResults.metrics?.rSquared) {
      g.append("text")
        .attr("x", width - 10)
        .attr("y", 20)
        .attr("text-anchor", "end")
        .style("font-size", "12px")
        .style("fill", "#9ca3af")
        .text(`R² = ${deconvolutionResults.metrics.rSquared.toFixed(4)}`);
    }
  }

  function drawResidualChart() {
    if (!residualContainer || !deconvolutionResults || !deconvolutionResults.residual) return;

    // Clear previous chart
    d3.select(residualContainer).selectAll("*").remove();

    const margin = { top: 40, right: 150, bottom: 50, left: 70 };
    const width = residualContainer.clientWidth - margin.left - margin.right;
    const height = 200 - margin.top - margin.bottom;

    const svg = d3
      .select(residualContainer)
      .append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom);

    const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);

    const wavenumbers = generateWavenumbers();
    const displayRange = [1000, 1500];

    // Filter to display range
    const startIdx = wavenumbers.findIndex((w) => w >= displayRange[0]);
    const endIdx = wavenumbers.findIndex((w) => w > displayRange[1]);
    const endIndex = endIdx === -1 ? wavenumbers.length - 1 : endIdx - 1;

    // Check if residual array exists and has data
    if (!Array.isArray(deconvolutionResults.residual) || deconvolutionResults.residual.length === 0) {
      console.error("Residual data is missing or empty");
      return;
    }

    // Prepare residual data - the residual is already for the analysis range only
    // So we need to map it correctly
    const residualLength = deconvolutionResults.residual.length;
    const expectedLength = endIndex - startIdx + 1;

    let residualData;
    if (residualLength === expectedLength) {
      // Residual is already trimmed to analysis range
      residualData = deconvolutionResults.residual.map((y: number, i: number) => ({
        x: wavenumbers[startIdx + i],
        y,
      }));
    } else if (residualLength === wavenumbers.length) {
      // Residual is full spectrum
      residualData = deconvolutionResults.residual
        .slice(startIdx, endIndex + 1)
        .map((y: number, i: number) => ({
          x: wavenumbers[startIdx + i],
          y,
        }));
    } else {
      console.error(`Unexpected residual length: ${residualLength}, expected ${expectedLength} or ${wavenumbers.length}`);
      return;
    }

    // Set up scales
    const xScale = d3.scaleLinear().domain(displayRange).range([0, width]);

    // Calculate max absolute value for symmetric scale
    const maxY = Number(d3.max(residualData, (d: any) => Math.abs(d.y))) || 0.01;
    console.log("Residual max value:", maxY, "Residual data points:", residualData.length);

    // Ensure minimum scale for visibility
    const scaleMax = Math.max(maxY * 1.2, 0.01);
    const yScale = d3
      .scaleLinear()
      .domain([-scaleMax, scaleMax])
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
      .call(d3.axisLeft(yScale).tickFormat(d3.format(".4f")))
      .append("text")
      .attr("transform", "rotate(-90)")
      .attr("y", -50)
      .attr("x", -height / 2)
      .attr("fill", "#9ca3af")
      .style("text-anchor", "middle")
      .text("Residual");

    // Add zero line
    g.append("line")
      .attr("x1", 0)
      .attr("x2", width)
      .attr("y1", yScale(0))
      .attr("y2", yScale(0))
      .attr("stroke", "#6b7280")
      .attr("stroke-width", 1)
      .attr("stroke-dasharray", "2,2");

    // Draw residual as area chart for better visibility
    const area = d3
      .area<{ x: number; y: number }>()
      .x((d) => xScale(d.x))
      .y0(yScale(0))
      .y1((d) => yScale(d.y))
      .curve(d3.curveMonotoneX);

    // Add filled area
    g.append("path")
      .datum(residualData)
      .attr("fill", "rgba(251, 146, 60, 0.3)") // Orange with transparency
      .attr("d", area);

    // Draw residual line
    const line = d3
      .line<{ x: number; y: number }>()
      .x((d) => xScale(d.x))
      .y((d) => yScale(d.y))
      .curve(d3.curveMonotoneX);

    g.append("path")
      .datum(residualData)
      .attr("fill", "none")
      .attr("stroke", "rgb(251, 146, 60)") // Orange
      .attr("stroke-width", 2)
      .attr("d", line);

    // Add title
    svg
      .append("text")
      .attr("x", (width + margin.left + margin.right) / 2)
      .attr("y", 20)
      .attr("text-anchor", "middle")
      .style("font-size", "14px")
      .style("fill", "#e5e7eb")
      .text("Fit Residuals (Measured - Reconstructed)");

    // Add RMSE annotation
    if (deconvolutionResults.metrics?.rmse) {
      g.append("text")
        .attr("x", width - 10)
        .attr("y", 20)
        .attr("text-anchor", "end")
        .style("font-size", "12px")
        .style("fill", "#9ca3af")
        .text(`RMSE = ${deconvolutionResults.metrics.rmse.toFixed(5)}`);
    }
  }

  // Redraw charts when data changes
  $effect(() => {
    if (normalizedData && sample && sampleStore) {
      drawOriginalChart();
      drawNormalizedChart();
    }
    if (deconvolutionResults && normalizedData) {
      drawDeconvolutionChart();
      drawResidualChart();
    }
  });

  // Handle window resize
  onMount(() => {
    const handleResize = () => {
      if (normalizedData && sample && sampleStore) {
        drawOriginalChart();
        drawNormalizedChart();
      }
      if (deconvolutionResults && normalizedData) {
        drawDeconvolutionChart();
        drawResidualChart();
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
            <div class="w-3 h-3 rounded-full" style="background-color: {item.color}"></div>
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

  <!-- Deconvolution Results Charts -->
  {#if deconvolutionResults}
    <!-- Component Contributions Chart -->
    <div class="bg-gray-900/50 rounded-lg p-4 border border-gray-700">
      <div bind:this={deconvolutionContainer} class="w-full h-[280px]"></div>
    </div>

    <!-- Residual Chart -->
    <div class="bg-gray-900/50 rounded-lg p-4 border border-gray-700">
      <div bind:this={residualContainer} class="w-full h-[200px]"></div>
    </div>
  {/if}
</div>
