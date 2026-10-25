// Project-wide rendering, status, and dirty-state coordination.

import { renderBlockList } from "./block-list.js";
import {
  AXES,
  ensureBoundaryArrays,
  isCellCountDerived,
  normalizedAxisLocks,
  normalizedDecomposition,
  normalizedFftAxis,
  normalizedPeriodicAxes,
  normalizedScalarCount,
} from "./block-model.js";
import { fitView } from "./camera.js";
import { scheduleHistoryCommit, updateHistoryButtons } from "./history.js";
import { renderInspector, renderSpacingPreview } from "./inspector.js";
import { scheduleGridPreviews } from "./previews.js";
import {
  renderCheckResults,
  renderDecompositionControls,
  renderDecompositionSummary,
  renderQualitySummary,
} from "./results.js";
import { els, firstSelectedBlockId, previewCache, selectedBlock, state, ui } from "./state.js";
import { renderScene, renderViewControls, resetInspectionState } from "./viewport.js";

export function setStatus(message) {
  els.status.textContent = message;
  els.status.title = message; // the full text, if three lines are not enough
}

export function renderAll() {
  if (!state.project) return;
  state.project.blocks = state.project.blocks || [];
  state.project.nscal = normalizedScalarCount(state.project.nscal ?? state.project.nScal);
  state.project.periodicAxes = normalizedPeriodicAxes(state.project);
  state.project.fftAxis = normalizedFftAxis(state.project);
  state.project.decomposition = normalizedDecomposition(state.project);
  for (const block of state.project.blocks) {
    ensureBoundaryArrays(block, state.project.nscal);
    block.axisLocks = normalizedAxisLocks(block.axisLocks);
  }
  const blockIds = new Set(state.project.blocks.map((block) => block.id));
  ui.hiddenBlockIds = new Set(
    [...ui.hiddenBlockIds].filter((id) => blockIds.has(id))
  );
  state.selectedIds = new Set([...state.selectedIds].filter((id) => blockIds.has(id)));
  if (state.selectedId != null && !blockIds.has(state.selectedId)) {
    state.selectedId = firstSelectedBlockId() ?? state.project.blocks[0]?.id ?? null;
  }
  if (state.selectedId == null && state.selectedIds.size) state.selectedId = firstSelectedBlockId();
  if (state.selectedId != null) state.selectedIds.add(state.selectedId);
  els.projectName.value = state.project.name ?? "snac-grid";
  els.scalarCount.value = state.project.nscal;
  els.inferConnectivity.checked = state.project.inferConnectivity !== false;
  state.project.writeExternalGrid = state.project.writeExternalGrid !== false;
  els.writeExternalGrid.checked = state.project.writeExternalGrid;
  if (!["grid", "both"].includes(state.project.externalGridSource)) state.project.externalGridSource = "grid";
  els.gridSource.value = state.project.externalGridSource;
  els.gridSource.disabled = !state.project.writeExternalGrid;
  if (!["single", "double"].includes(state.project.externalGridPrecision)) state.project.externalGridPrecision = "double";
  els.gridPrecision.value = state.project.externalGridPrecision;
  els.gridPrecision.disabled = !state.project.writeExternalGrid;
  els.renumberBlocks.disabled = !state.project.blocks.length;
  els.addAdjacent.disabled = !selectedBlock();
  for (const button of [els.duplicateArray, els.mirrorBlocks, els.alignBlocks, els.snapBlock]) {
    button.disabled = !selectedBlock();
  }
  renderBlockList();
  renderDecompositionControls();
  renderQualitySummary();
  renderInspector();
  renderScene();
  renderSpacingPreview();
  renderCheckResults();
  updateHistoryButtons();
  renderViewControls();
  scheduleGridPreviews();
  if (window.lucide) window.lucide.createIcons();
}

// Make *project* the edited project and clear everything derived from the
// previous one: selection, check, repair and decomposition results, previews,
// and the viewport's hidden blocks and clipping. Callers handle undo history.
export function loadProject(project, { fit = true } = {}) {
  state.project = project;
  state.selectedId = project.blocks[0]?.id ?? null;
  state.selectedIds = new Set(state.selectedId == null ? [] : [state.selectedId]);
  state.activeDiagnostic = null;
  state.lastCheck = null;
  state.pendingRepair = null;
  state.decompositionResult = null;
  resetInspectionState();
  previewCache.clear();
  state.previewSequence += 1;
  renderAll();
  if (fit) fitView();
}

export function markProjectDirty(blockId = null) {
  if (blockId != null) previewCache.delete(blockId);
  state.pendingRepair = null;
  state.decompositionResult = null;
  clearCheckState();
  renderDecompositionSummary();
  scheduleHistoryCommit();
}

export function clearCheckState() {
  state.lastCheck = null;
  renderCheckResults();
  renderQualitySummary();
  setStatus("");
}

export function decompositionSettingEdited() {
  state.pendingRepair = null;
  state.decompositionResult = null;
  state.lastCheck = null;
}

export function projectRequestSnapshot() {
  // Previews fill in cell counts that the server derives from the grading
  // controls. The server recomputes them for every request, so they are left
  // out here and a preview landing mid-request does not discard its result.
  return JSON.stringify({
    ...state.project,
    blocks: state.project.blocks.map((block) => ({
      ...block,
      ng: block.ng.map((value, index) => (isCellCountDerived(block.axes?.[AXES[index]] ?? {}) ? null : value)),
    })),
  });
}

export function projectResponseIsCurrent(snapshot, operation) {
  if (projectRequestSnapshot() === snapshot) return true;
  setStatus(`${operation} result discarded because the project changed`);
  return false;
}
