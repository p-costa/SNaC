// Check results, grid-quality summary, and MPI decomposition panels.

import { balanceDecomposition } from "./actions.js";
import { AXES, blockCenter } from "./block-model.js";
import { frameBlocks } from "./camera.js";
import { markProjectDirty, renderAll, setStatus } from "./core.js";
import { escapeHtml, formatBytes, formatInteger, formatNumber } from "./format.js";
import { els, selectedBlock, state, ui } from "./state.js";

export function renderDecompositionControls() {
  if (!state.project?.decomposition) return;
  const spec = state.project.decomposition;
  els.targetRanks.value = spec.targetRanks;
  els.minLocalCells.value = spec.minLocalCells;
  els.maxLocalAspect.value = spec.maxLocalAspect;
  for (const button of els.decompositionMode.querySelectorAll("button")) {
    button.classList.toggle("active", button.dataset.mode === spec.mode);
  }
  const fftIndex = AXES.indexOf(state.project.fftAxis);
  els.fftAxis.value = state.project.fftAxis ?? "";
  for (const input of els.decompositionAxes.querySelectorAll("input")) {
    const axisIndex = Number(input.dataset.axis);
    input.checked = spec.axes[axisIndex] && axisIndex !== fftIndex;
    input.disabled = axisIndex === fftIndex;
    input.closest("label").title = axisIndex === fftIndex ? "SNaC never splits the FFT axis" : "";
  }
  els.balanceDecomposition.disabled = !state.project.blocks.length || spec.targetRanks < state.project.blocks.length;
  renderDecompositionSummary();
}

export function renderDecompositionSummary() {
  if (!state.decompositionResult) {
    const target = state.project.decomposition.targetRanks;
    els.decompositionSummary.textContent = target > 0 ? "Balance to apply an exact decomposition." : "Set total ranks to enable balancing.";
    return;
  }
  const result = state.decompositionResult;
  if (!result.ok) {
    const nearby = rankButtons(
      "Nearby feasible totals",
      (result.nearbyRanks ?? []).map((ranks) => ({ ranks }))
    );
    els.decompositionSummary.innerHTML = `${(result.errors ?? []).map((message) => `<div class="error-text">${escapeHtml(message)}</div>`).join("")}${nearby}`;
  } else {
    const axes = result.activeAxes?.length ? result.activeAxes.map((axis) => axis.toUpperCase()).join("") : "serial";
    const quality = result.optimal ? "best found" : "bounded search";
    // Interface constraints can make a nearby total balance much better.
    const better = rankButtons("Better balanced totals", result.betterRanks ?? []);
    els.decompositionSummary.innerHTML = `
      <div><strong>${formatInteger(result.totalRanks)}</strong> ranks · ${axes} · ${quality}</div>
      <div>${formatInteger(result.minCells)}-${formatInteger(result.maxCells)} cells/rank · imbalance <strong>${formatNumber(result.imbalance)}</strong></div>
      <div>minimum local edge <strong>${formatInteger(result.minLocalExtent)}</strong> cells · maximum aspect <strong>${formatNumber(result.maxLocalAspect)}</strong></div>
      ${(result.blocks ?? [])
        .map((block) => `<div>B${block.blockId}: ${block.dims.join(" × ")} · ${block.ranks} ranks</div>`)
        .join("")}
      ${better}
    `;
  }
  for (const button of els.decompositionSummary.querySelectorAll(".nearby-rank")) {
    button.addEventListener("click", () => {
      state.project.decomposition.targetRanks = Number(button.dataset.ranks);
      markProjectDirty();
      renderDecompositionControls();
      void balanceDecomposition();
    });
  }
}

// Buttons that set the total rank count and balance again; items with an
// imbalance show it next to the count.
function rankButtons(title, items) {
  if (!items.length) return "";
  const buttons = items
    .map((item) => {
      const label = item.imbalance == null
        ? formatInteger(item.ranks)
        : `${formatInteger(item.ranks)} · imbalance ${formatNumber(item.imbalance)}`;
      return `<button type="button" class="nearby-rank" data-ranks="${item.ranks}">${label}</button>`;
    })
    .join("");
  return `<div class="nearby-ranks"><span>${title}</span><div>${buttons}</div></div>`;
}

export function renderQualitySummary() {
  const quality = state.lastCheck?.quality;
  const storage = state.lastCheck?.storage;
  const disabled = !state.project.blocks.length;
  els.downloadReportJson.disabled = disabled;
  els.downloadReportMarkdown.disabled = disabled;
  if (!quality) {
    els.qualitySummary.innerHTML = `<div class="quality-empty">${escapeHtml(state.lastCheck?.qualityError ?? "Not evaluated")}</div>`;
    return;
  }
  const summary = quality.summary;
  const rows = [
    ["Cells", formatInteger(summary.totalCells)],
    ["MPI ranks", formatInteger(summary.totalRanks)],
    ["Cells/rank", `${formatInteger(summary.cellsPerRankMin)} to ${formatInteger(summary.cellsPerRankMax)}`],
    ["Load ratio", formatNumber(summary.loadImbalance)],
    ["Spacing", `${formatNumber(summary.spacingMin)} to ${formatNumber(summary.spacingMax)}`],
    ["Adjacent ratio", formatNumber(summary.maxAdjacentRatio)],
    ["Cell aspect", formatNumber(summary.worstCellAspect)],
    ["Interface ratio", formatNumber(summary.maxInterfaceRatio)],
    ["Coordinates", formatBytes(summary.coordinateBytes)],
    ["SNaC grids", formatBytes(storage?.snacBinaryBytesTotal ?? 0)],
  ];
  els.qualitySummary.innerHTML = rows
    .map(([label, value]) => `<div><span>${label}</span><strong>${value}</strong></div>`)
    .join("");
}

export function renderCheckResults() {
  if (!els.checkResults) return;
  const errors = state.lastCheck?.errors ?? [];
  const warnings = state.lastCheck?.warnings ?? [];
  const diagnostics = state.lastCheck?.diagnostics ?? [
    ...errors.map((message) => ({ severity: "error", message, locations: [] })),
    ...warnings.map((message) => ({ severity: "warning", message, locations: [] })),
  ];
  const noBlocks = !state.project.blocks.length;
  els.exportCase.disabled = noBlocks || errors.length > 0;
  els.exportCase.title = noBlocks
    ? "Add a block before writing"
    : errors.length
      ? "Fix the check errors shown in the 3D view before writing"
      : "Write case files";
  if (!state.lastCheck || (!errors.length && !warnings.length)) {
    state.activeDiagnostic = null;
    els.checkResults.innerHTML = "";
    els.issuesPanel.hidden = true;
    return;
  }
  if (
    state.activeDiagnostic &&
    !diagnostics.some((item) => item.code === state.activeDiagnostic.code && item.message === state.activeDiagnostic.message)
  ) {
    state.activeDiagnostic = null;
  }
  // A new set of issues opens the panel again.
  const signature = diagnostics.map((item) => item.message).join("\n");
  if (signature !== els.issuesPanel.dataset.signature) ui.issuesCollapsed = false;
  els.issuesPanel.dataset.signature = signature;
  const count = (value, noun) => `${value} ${noun}${value === 1 ? "" : "s"}`;
  els.issuesTitle.textContent = [errors.length && count(errors.length, "error"), warnings.length && count(warnings.length, "warning")]
    .filter(Boolean)
    .join(" · ");
  els.issuesPanel.hidden = false;
  els.issuesPanel.classList.toggle("has-errors", errors.length > 0);
  els.issuesPanel.classList.toggle("collapsed", ui.issuesCollapsed);
  els.toggleIssues.setAttribute("aria-expanded", String(!ui.issuesCollapsed));
  els.checkResults.innerHTML = diagnostics
    .map((diagnostic, index) => {
      const actionable = (diagnostic.locations ?? []).some((location) =>
        state.project.blocks.some((block) => block.id === location.blockId)
      );
      const tag = actionable ? "button" : "div";
      const attributes = actionable
        ? `type="button" data-diagnostic-index="${index}" title="Show affected block or face"`
        : "";
      const icon = actionable ? `<i data-lucide="locate-fixed"></i>` : "";
      return `<${tag} class="check-line ${diagnostic.severity}${actionable ? " actionable" : ""}" ${attributes}>${icon}<span>${escapeHtml(diagnostic.message)}</span></${tag}>`;
    })
    .join("");
  for (const button of els.checkResults.querySelectorAll("[data-diagnostic-index]")) {
    button.addEventListener("click", () => focusDiagnostic(diagnostics[Number(button.dataset.diagnosticIndex)]));
  }
  if (window.lucide) window.lucide.createIcons();
}

export function focusDiagnostic(diagnostic) {
  const blockIds = (diagnostic.locations ?? [])
    .map((location) => location.blockId)
    .filter((blockId) => state.project.blocks.some((block) => block.id === blockId));
  if (!blockIds.length) return;
  state.activeDiagnostic = diagnostic;
  state.selectedId = blockIds[0];
  state.selectedIds = new Set(blockIds);
  if (AXES.includes(diagnostic.axis)) state.currentAxis = diagnostic.axis;
  // Zoom to the affected blocks so the outlined faces are visible.
  const visible = state.project.blocks.filter((block) => blockIds.includes(block.id) && !ui.hiddenBlockIds.has(block.id));
  if (!frameBlocks(visible)) {
    const block = selectedBlock();
    if (block) ui.controls.target.set(...blockCenter(block));
  }
  renderAll();
  setStatus(`Showing ${diagnostic.code?.replaceAll("_", " ") ?? "grid issue"}`);
}
