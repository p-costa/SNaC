// Entry point: element binding, static event wiring, and startup.

import {
  addAdjacentBlock,
  addFreeBlock,
  applyBoundaryPreset,
  applyGeometryOperation,
  applyGridRepair,
  balanceDecomposition,
  checkProject,
  closeRepairDialog,
  downloadCaseReport,
  downloadProjectJson,
  exportCase,
  importCase,
  openGridRepair,
  openProjectJson,
  removeSelectedBlock,
  resetProject,
  updateStructure,
} from "./actions.js";
import { getJson } from "./api.js";
import { renumberBlocksFromOrder } from "./block-list.js";
import { AXES, ensureBoundaryArrays, normalizedScalarCount } from "./block-model.js";
import { fitView, focusSelection, resizeRenderer, setAxisView } from "./camera.js";
import {
  clearCheckState,
  decompositionSettingEdited,
  loadProject,
  markProjectDirty,
  setStatus,
} from "./core.js";
import { bindNumberField, toCamel } from "./format.js";
import {
  claimAutosave,
  discardRecovery,
  flushAutosave,
  loadAutosave,
  redoProject,
  resetHistory,
  restoreRecovery,
  scheduleHistoryCommit,
  showRecoveryDialog,
  undoProject,
} from "./history.js";
import { renderInspector } from "./inspector.js";
import {
  applySelectedLibraryItem,
  deleteSelectedLibraryItem,
  exportLibraryJson,
  importLibraryJson,
  renameSelectedLibraryItem,
  renderLibraryDialog,
  renderLibrarySelection,
  saveCurrentLibraryItem,
  showLibraryDialog,
} from "./library-ui.js";
import { scheduleGridPreviews, schedulePreview } from "./previews.js";
import { renderCheckResults, renderDecompositionControls, renderDecompositionSummary, tolerancePercent } from "./results.js";
import { builtInProjectTemplates } from "./snac-templates.js";
import { els, state, ui } from "./state.js";
import {
  animate,
  beginSceneSelection,
  cancelSceneSelection,
  completeSceneSelection,
  hideSelectedBlocks,
  initScene,
  isolateSelectedBlocks,
  renderScene,
  renderViewControls,
  setBlockTool,
  setClippingAxis,
  showAllBlocks,
  toggleClipping,
  updateClippingPlane,
} from "./viewport.js";

async function bootstrap() {
  bindElements();
  initScene();
  bindStaticEvents();
  let startupError = null;
  let project;
  try {
    project = await getJson("/api/sample");
  } catch (error) {
    startupError = error;
    project = JSON.parse(JSON.stringify(builtInProjectTemplates()[0].project));
  }
  loadProject(project);
  const recovery = await loadAutosave();
  claimAutosave();
  resetHistory(false);
  if (recovery) showRecoveryDialog(recovery);
  schedulePreview(state.selectedId, 0);
  animate();
  if (startupError) setStatus(`Could not load the sample project (${startupError.message}); started from a single block`);
}

function bindElements() {
  for (const id of [
    "project-name",
    "scalar-count",
    "output-dir",
    "export-case",
    "check-project",
    "update-structure",
    "repair-grids",
    "download-json",
    "import-case",
    "open-json",
    "open-library",
    "undo-project",
    "redo-project",
    "reset-project",
    "status",
    "block-list",
    "renumber-blocks",
    "add-block",
    "add-adjacent",
    "adjacent-face",
    "adjacent-size",
    "geometry-axis",
    "array-count",
    "array-gap",
    "duplicate-array",
    "mirror-plane",
    "mirror-blocks",
    "align-mode",
    "align-blocks",
    "snap-face",
    "snap-block",
    "target-ranks",
    "min-local-cells",
    "max-local-aspect",
    "ng-tolerance",
    "decomposition-mode",
    "decomposition-axes",
    "fft-axis",
    "balance-decomposition",
    "decomposition-summary",
    "download-report-json",
    "download-report-markdown",
    "quality-summary",
    "scene",
    "fit-view",
    "focus-selection",
    "hide-selected",
    "isolate-selected",
    "show-all-blocks",
    "toggle-clipping",
    "clip-tools",
    "clip-axis",
    "clip-position",
    "flip-clipping",
    "toggle-grid",
    "tool-select",
    "tool-move",
    "tool-scale",
    "toggle-boundaries",
    "toggle-partitions",
    "grid-lines-mode",
    "toggle-spacing-jumps",
    "block-title",
    "remove-block",
    "block-name",
    "inivel",
    "periodic-grid",
    "geometry-grid",
    "resolution-grid",
    "axis-tabs",
    "axis-editor",
    "spacing-preview",
    "spacing-stats",
    "infer-connectivity",
    "write-external-grid",
    "grid-source",
    "grid-precision",
    "boundary-table",
    "bc-preset",
    "bc-preset-face",
    "bc-preset-u",
    "bc-preset-v",
    "bc-preset-w",
    "bc-preset-p",
    "apply-bc-preset",
    "check-results",
    "issues-panel",
    "issues-title",
    "toggle-issues",
    "repair-dialog",
    "repair-summary",
    "cancel-repair",
    "confirm-repair",
    "recovery-dialog",
    "recovery-summary",
    "discard-recovery",
    "restore-recovery",
    "library-dialog",
    "close-library",
    "library-tabs",
    "library-name",
    "library-items",
    "library-summary",
    "apply-library-item",
    "save-library-item",
    "rename-library-item",
    "delete-library-item",
    "import-library",
    "export-library",
  ]) {
    els[toCamel(id)] = document.getElementById(id);
  }
}

function bindStaticEvents() {
  els.projectName.addEventListener("input", () => {
    state.project.name = els.projectName.value;
    markProjectDirty();
  });
  els.scalarCount.addEventListener("input", () => {
    state.project.nscal = normalizedScalarCount(els.scalarCount.value);
    for (const block of state.project.blocks) ensureBoundaryArrays(block, state.project.nscal);
    renderInspector();
    markProjectDirty();
  });
  els.inferConnectivity.addEventListener("change", () => {
    state.project.inferConnectivity = els.inferConnectivity.checked;
    markProjectDirty();
  });
  els.writeExternalGrid.addEventListener("change", () => {
    state.project.writeExternalGrid = els.writeExternalGrid.checked;
    els.gridSource.disabled = !state.project.writeExternalGrid;
    els.gridPrecision.disabled = !state.project.writeExternalGrid;
    markProjectDirty();
  });
  els.gridSource.addEventListener("change", () => {
    state.project.externalGridSource = els.gridSource.value;
    markProjectDirty();
  });
  els.gridPrecision.addEventListener("change", () => {
    state.project.externalGridPrecision = els.gridPrecision.value;
    markProjectDirty();
  });
  els.addBlock.addEventListener("click", addFreeBlock);
  els.renumberBlocks.addEventListener("click", renumberBlocksFromOrder);
  els.addAdjacent.addEventListener("click", addAdjacentBlock);
  els.duplicateArray.addEventListener("click", () => applyGeometryOperation("duplicate"));
  els.mirrorBlocks.addEventListener("click", () => applyGeometryOperation("mirror"));
  els.alignBlocks.addEventListener("click", () => applyGeometryOperation("align"));
  els.snapBlock.addEventListener("click", () => applyGeometryOperation("snap"));
  els.removeBlock.addEventListener("click", removeSelectedBlock);
  els.fitView.addEventListener("click", fitView);
  els.focusSelection.addEventListener("click", focusSelection);
  els.hideSelected.addEventListener("click", hideSelectedBlocks);
  els.isolateSelected.addEventListener("click", isolateSelectedBlocks);
  els.showAllBlocks.addEventListener("click", showAllBlocks);
  els.toggleClipping.addEventListener("click", toggleClipping);
  for (const button of document.querySelectorAll("[data-axis-view]")) {
    button.addEventListener("click", () => setAxisView(button.dataset.axisView));
  }
  for (const button of els.clipAxis.querySelectorAll("button")) {
    button.addEventListener("click", () => setClippingAxis(button.dataset.clipAxis));
  }
  els.clipPosition.addEventListener("input", () => {
    ui.clippingFraction = Number(els.clipPosition.value);
    updateClippingPlane();
  });
  els.flipClipping.addEventListener("click", () => {
    ui.clippingDirection *= -1;
    updateClippingPlane();
    renderViewControls();
  });
  els.toggleGrid.addEventListener("click", () => {
    ui.groundGridVisible = !ui.groundGridVisible;
    ui.gridHelper.visible = ui.groundGridVisible;
    els.toggleGrid.classList.toggle("active", ui.groundGridVisible);
    els.toggleGrid.setAttribute("aria-pressed", String(ui.groundGridVisible));
  });
  els.toolSelect.addEventListener("click", () => setBlockTool("select"));
  els.toolMove.addEventListener("click", () => setBlockTool("translate"));
  els.toolScale.addEventListener("click", () => setBlockTool("scale"));
  els.toggleBoundaries.addEventListener("click", () => {
    ui.boundaryColorsVisible = !ui.boundaryColorsVisible;
    renderScene();
  });
  els.togglePartitions.addEventListener("click", () => {
    ui.partitionPlanesVisible = !ui.partitionPlanesVisible;
    renderScene();
  });
  els.gridLinesMode.addEventListener("change", () => {
    ui.gridLinesMode = els.gridLinesMode.value;
    renderScene();
    scheduleGridPreviews(0);
  });
  els.toggleSpacingJumps.addEventListener("click", () => {
    ui.spacingJumpsVisible = !ui.spacingJumpsVisible;
    renderScene();
  });
  els.exportCase.addEventListener("click", exportCase);
  els.checkProject.addEventListener("click", () => checkProject());
  els.updateStructure.addEventListener("click", updateStructure);
  els.repairGrids.addEventListener("click", openGridRepair);
  bindNumberField(els.targetRanks, {
    integer: true,
    valid: (value) => value >= 0,
    read: () => state.project.decomposition.targetRanks,
    write: (value) => {
      state.project.decomposition.targetRanks = value;
      decompositionSettingEdited();
      els.balanceDecomposition.disabled = !state.project.blocks.length || value < state.project.blocks.length;
    },
  });
  bindNumberField(els.minLocalCells, {
    integer: true,
    valid: (value) => value >= 1,
    read: () => state.project.decomposition.minLocalCells,
    write: (value) => {
      state.project.decomposition.minLocalCells = value;
      decompositionSettingEdited();
    },
  });
  bindNumberField(els.maxLocalAspect, {
    valid: (value) => value === 0 || value >= 1, // zero disables the limit
    read: () => state.project.decomposition.maxLocalAspect,
    write: (value) => {
      state.project.decomposition.maxLocalAspect = value;
      decompositionSettingEdited();
    },
  });
  bindNumberField(els.ngTolerance, {
    valid: (value) => value >= 0 && value <= 50,
    read: () => tolerancePercent(state.project.decomposition.ngTolerance),
    write: (value) => {
      state.project.decomposition.ngTolerance = value / 100;
      decompositionSettingEdited();
    },
  });
  for (const input of [els.targetRanks, els.minLocalCells, els.maxLocalAspect, els.ngTolerance]) {
    input.addEventListener("change", () => {
      clearCheckState();
      renderDecompositionSummary();
      scheduleHistoryCommit();
    });
  }
  for (const button of els.decompositionMode.querySelectorAll("button")) {
    button.addEventListener("click", () => {
      state.project.decomposition.mode = button.dataset.mode;
      markProjectDirty();
      renderDecompositionControls();
    });
  }
  for (const input of els.decompositionAxes.querySelectorAll("input")) {
    input.addEventListener("change", () => {
      const axisIndex = Number(input.dataset.axis);
      const fftIndex = AXES.indexOf(state.project.fftAxis);
      // The FFT axis is never partitioned, so it does not count as enabled.
      const enabledCount = state.project.decomposition.axes.filter((enabled, index) => enabled && index !== fftIndex).length;
      if (!input.checked && enabledCount === 1) {
        input.checked = true;
        setStatus("Keep at least one partition axis enabled");
        return;
      }
      state.project.decomposition.axes[axisIndex] = input.checked;
      markProjectDirty();
      renderDecompositionControls();
    });
  }
  els.fftAxis.addEventListener("change", () => {
    state.project.fftAxis = els.fftAxis.value || null;
    const fftIndex = AXES.indexOf(state.project.fftAxis);
    if (fftIndex >= 0) {
      // SNaC never splits the FFT axis; keep at least one other axis available.
      const axes = state.project.decomposition.axes;
      axes[fftIndex] = false;
      if (!axes.some(Boolean)) state.project.decomposition.axes = AXES.map((_, index) => index !== fftIndex);
    }
    markProjectDirty();
    renderDecompositionControls();
  });
  els.balanceDecomposition.addEventListener("click", balanceDecomposition);
  els.cancelRepair.addEventListener("click", closeRepairDialog);
  els.confirmRepair.addEventListener("click", applyGridRepair);
  els.applyBcPreset.addEventListener("click", applyBoundaryPreset);
  els.toggleIssues.addEventListener("click", () => {
    ui.issuesCollapsed = !ui.issuesCollapsed;
    renderCheckResults();
  });
  els.downloadJson.addEventListener("click", downloadProjectJson);
  els.downloadReportJson.addEventListener("click", () => downloadCaseReport("json"));
  els.downloadReportMarkdown.addEventListener("click", () => downloadCaseReport("markdown"));
  els.importCase.addEventListener("click", importCase);
  els.openJson.addEventListener("change", openProjectJson);
  els.openLibrary.addEventListener("click", showLibraryDialog);
  els.closeLibrary.addEventListener("click", () => els.libraryDialog.close());
  for (const button of els.libraryTabs.querySelectorAll("button")) {
    button.addEventListener("click", () => {
      state.currentLibraryTab = button.dataset.libraryTab;
      state.selectedLibraryItem = null;
      renderLibraryDialog();
    });
  }
  els.libraryItems.addEventListener("change", () => {
    state.selectedLibraryItem = els.libraryItems.value;
    renderLibrarySelection();
  });
  els.applyLibraryItem.addEventListener("click", applySelectedLibraryItem);
  els.saveLibraryItem.addEventListener("click", saveCurrentLibraryItem);
  els.renameLibraryItem.addEventListener("click", renameSelectedLibraryItem);
  els.deleteLibraryItem.addEventListener("click", deleteSelectedLibraryItem);
  els.importLibrary.addEventListener("change", importLibraryJson);
  els.exportLibrary.addEventListener("click", exportLibraryJson);
  els.undoProject.addEventListener("click", undoProject);
  els.redoProject.addEventListener("click", redoProject);
  els.resetProject.addEventListener("click", resetProject);
  els.discardRecovery.addEventListener("click", discardRecovery);
  els.restoreRecovery.addEventListener("click", restoreRecovery);
  els.scene.addEventListener("pointerdown", beginSceneSelection);
  els.scene.addEventListener("pointerup", completeSceneSelection);
  els.scene.addEventListener("pointercancel", cancelSceneSelection);
  window.addEventListener("resize", resizeRenderer);
  window.addEventListener("pagehide", flushAutosave);
  document.addEventListener("visibilitychange", () => {
    if (document.visibilityState === "hidden") flushAutosave();
  });
}

bootstrap();
