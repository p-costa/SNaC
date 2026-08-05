import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { TransformControls } from "three/addons/controls/TransformControls.js";
import { postJson } from "./api.js";
import {
  AXES,
  FACE_ORDER,
  blockCenter,
  blockSize,
  defaultBlock,
  ensureAxis,
  ensureBoundaryArrays,
  isCellCountDerived,
  nextBlockId,
  normalizedAxisLocks,
  normalizedDecomposition,
  normalizedPeriodicAxes,
  normalizedScalarCount,
  propagateMpiPartition,
  rescaleExplicitFaces,
  resetFriendBoundaries,
} from "./block-model.js";
import {
  SnapshotHistory,
  clearStoredProject,
  loadStoredProject,
  saveStoredProject,
} from "./project-state.js";
import {
  LIBRARY_SCHEMA_VERSION,
  applyAxisPreset,
  captureAxisPreset,
  captureProjectTemplate,
  loadLibrary,
  mergeLibraries,
  normalizeLibrary,
  removeLibraryItem,
  renameLibraryItem,
  replaceLibraryItem,
  saveLibrary,
} from "./library.js";
import {
  axisLabel,
  clearDisposableGroup,
  partitionFaceIndices,
  planeOutlinePoints,
  projectBox,
} from "./scene-helpers.js";
import {
  SNAC_GRID_FUNCTIONS as GRID_FUNCTIONS,
  SNAC_INITIAL_VELOCITY_FIELDS as INIT_FIELDS,
} from "./snac.js";
import { builtInAxisPresets, builtInProjectTemplates } from "./snac-templates.js";

const PROFILE_OPTIONS = [
  ["geometric", "Geometric"],
  ["tanh", "tanh"],
  ["erf", "erf"],
];
const MONOTONE_CONTROL_LABELS = {
  n: "Cells",
  ratio: "Upper/lower ratio",
  cell_ratio: "Cell-to-cell ratio",
  width_start: "Lower spacing",
  width_end: "Upper spacing",
};
const MONOTONE_CONTROL_PAIRS = [
  ["n", "ratio"],
  ["n", "cell_ratio"],
  ["n", "width_start"],
  ["n", "width_end"],
  ["ratio", "cell_ratio"],
  ["ratio", "width_start"],
  ["ratio", "width_end"],
  ["cell_ratio", "width_start"],
  ["cell_ratio", "width_end"],
  ["width_start", "width_end"],
];
const COLORS = [0x43c989, 0x76a9ff, 0xd28a3e, 0xc789e8, 0xe46868, 0x82c6c2, 0xd6bf5b, 0x9ba7ff];
const BOUNDARY_COLORS = { D: 0xe46868, F: 0x43c989, N: 0x76a9ff, mixed: 0xd6bf5b };

let project = null;
let selectedId = 1;
let selectedIds = new Set([1]);
let currentAxis = "x";
let groundGridVisible = true;
let boundaryColorsVisible = true;
let partitionPlanesVisible = true;
let gridLinesVisible = true;
let spacingJumpsVisible = true;
let blockTool = "select";
let lastCheck = null;
let activeDiagnostic = null;
let pendingRepair = null;
let decompositionResult = null;
let previewTimer = null;
let previewSequence = 0;
const previewCache = new Map();
let historyTimer = null;
const projectHistory = new SnapshotHistory(100);
let pendingRecovery = null;
let userLibrary = loadLibrary();
let currentLibraryTab = "axisPresets";
let selectedLibraryItem = "builtin:builtin-uniform-32";
const AUTOSAVE_KEY = "snac-grid-generator-autosave-v1";

const els = {};
let scene;
let camera;
let renderer;
let controls;
let raycaster;
let pointer;
let blockGroup;
let boundaryFaceGroup;
let gridHelper;
let axesGroup;
let partitionGroup;
let selectedGridGroup;
let spacingJumpGroup;
let diagnosticGroup;
let transformControls;
let selectedBlockMesh = null;
let transformDrag = null;
let isTransformDragging = false;

bootstrap();

async function bootstrap() {
  bindElements();
  initScene();
  bindStaticEvents();
  const response = await fetch("/api/sample");
  project = await response.json();
  selectedId = project.blocks[0]?.id ?? null;
  selectedIds = new Set(selectedId == null ? [] : [selectedId]);
  renderAll();
  const recovery = loadAutosave();
  resetHistory(false);
  if (recovery) showRecoveryDialog(recovery);
  schedulePreview(selectedId, 0);
  animate();
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
    "decomposition-mode",
    "decomposition-axes",
    "balance-decomposition",
    "decomposition-summary",
    "download-report-json",
    "download-report-markdown",
    "quality-summary",
    "scene",
    "fit-view",
    "toggle-grid",
    "tool-select",
    "tool-move",
    "tool-scale",
    "toggle-boundaries",
    "toggle-partitions",
    "toggle-grid-lines",
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
    "boundary-table",
    "bc-preset",
    "bc-preset-face",
    "bc-preset-u",
    "bc-preset-v",
    "bc-preset-w",
    "bc-preset-p",
    "apply-bc-preset",
    "check-results",
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
    project.name = els.projectName.value;
    markProjectDirty();
  });
  els.scalarCount.addEventListener("input", () => {
    project.nscal = normalizedScalarCount(els.scalarCount.value);
    for (const block of project.blocks) ensureBoundaryArrays(block, project.nscal);
    renderInspector();
    markProjectDirty();
  });
  els.inferConnectivity.addEventListener("change", () => {
    project.inferConnectivity = els.inferConnectivity.checked;
    markProjectDirty();
  });
  els.writeExternalGrid.addEventListener("change", () => {
    project.writeExternalGrid = els.writeExternalGrid.checked;
    els.gridSource.disabled = !project.writeExternalGrid;
    markProjectDirty();
  });
  els.gridSource.addEventListener("change", () => {
    project.externalGridSource = els.gridSource.value;
    markProjectDirty();
  });
  els.addBlock.addEventListener("click", addFreeBlock);
  els.addAdjacent.addEventListener("click", addAdjacentBlock);
  els.duplicateArray.addEventListener("click", () => applyGeometryOperation("duplicate"));
  els.mirrorBlocks.addEventListener("click", () => applyGeometryOperation("mirror"));
  els.alignBlocks.addEventListener("click", () => applyGeometryOperation("align"));
  els.snapBlock.addEventListener("click", () => applyGeometryOperation("snap"));
  els.removeBlock.addEventListener("click", removeSelectedBlock);
  els.fitView.addEventListener("click", fitView);
  els.toggleGrid.addEventListener("click", () => {
    groundGridVisible = !groundGridVisible;
    gridHelper.visible = groundGridVisible;
    els.toggleGrid.classList.toggle("active", groundGridVisible);
  });
  els.toolSelect.addEventListener("click", () => setBlockTool("select"));
  els.toolMove.addEventListener("click", () => setBlockTool("translate"));
  els.toolScale.addEventListener("click", () => setBlockTool("scale"));
  els.toggleBoundaries.addEventListener("click", () => {
    boundaryColorsVisible = !boundaryColorsVisible;
    renderScene();
  });
  els.togglePartitions.addEventListener("click", () => {
    partitionPlanesVisible = !partitionPlanesVisible;
    renderScene();
  });
  els.toggleGridLines.addEventListener("click", () => {
    gridLinesVisible = !gridLinesVisible;
    renderScene();
  });
  els.toggleSpacingJumps.addEventListener("click", () => {
    spacingJumpsVisible = !spacingJumpsVisible;
    renderScene();
  });
  els.exportCase.addEventListener("click", exportCase);
  els.checkProject.addEventListener("click", () => {
    if (project.inferConnectivity) updateStructure({ successMessage: "Grid checks passed" });
    else checkProject();
  });
  els.updateStructure.addEventListener("click", updateStructure);
  els.repairGrids.addEventListener("click", previewGridRepair);
  els.targetRanks.addEventListener("input", () => {
    project.decomposition.targetRanks = Math.max(0, Math.round(Number(els.targetRanks.value) || 0));
    decompositionSettingEdited();
    els.balanceDecomposition.disabled = !project.blocks.length || project.decomposition.targetRanks < project.blocks.length;
  });
  els.minLocalCells.addEventListener("input", () => {
    project.decomposition.minLocalCells = Math.max(1, Math.round(Number(els.minLocalCells.value) || 4));
    decompositionSettingEdited();
  });
  els.maxLocalAspect.addEventListener("input", () => {
    project.decomposition.maxLocalAspect = Math.max(0, Number(els.maxLocalAspect.value) || 0);
    decompositionSettingEdited();
  });
  for (const input of [els.targetRanks, els.minLocalCells, els.maxLocalAspect]) {
    input.addEventListener("change", () => {
      clearCheckState();
      renderDecompositionSummary();
      scheduleHistoryCommit();
    });
  }
  for (const button of els.decompositionMode.querySelectorAll("button")) {
    button.addEventListener("click", () => {
      project.decomposition.mode = button.dataset.mode;
      markProjectDirty();
      renderDecompositionControls();
    });
  }
  for (const input of els.decompositionAxes.querySelectorAll("input")) {
    input.addEventListener("change", () => {
      const axisIndex = Number(input.dataset.axis);
      const enabledCount = project.decomposition.axes.filter(Boolean).length;
      if (!input.checked && enabledCount === 1) {
        input.checked = true;
        setStatus("Keep at least one partition axis enabled");
        return;
      }
      project.decomposition.axes[axisIndex] = input.checked;
      markProjectDirty();
      renderDecompositionControls();
    });
  }
  els.balanceDecomposition.addEventListener("click", balanceDecomposition);
  els.cancelRepair.addEventListener("click", closeRepairDialog);
  els.confirmRepair.addEventListener("click", applyGridRepair);
  els.applyBcPreset.addEventListener("click", applyBoundaryPreset);
  els.downloadJson.addEventListener("click", downloadProjectJson);
  els.downloadReportJson.addEventListener("click", () => downloadCaseReport("json"));
  els.downloadReportMarkdown.addEventListener("click", () => downloadCaseReport("markdown"));
  els.importCase.addEventListener("click", importCase);
  els.openJson.addEventListener("change", openProjectJson);
  els.openLibrary.addEventListener("click", showLibraryDialog);
  els.closeLibrary.addEventListener("click", () => els.libraryDialog.close());
  for (const button of els.libraryTabs.querySelectorAll("button")) {
    button.addEventListener("click", () => {
      currentLibraryTab = button.dataset.libraryTab;
      selectedLibraryItem = null;
      renderLibraryDialog();
    });
  }
  els.libraryItems.addEventListener("change", () => {
    selectedLibraryItem = els.libraryItems.value;
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
  els.scene.addEventListener("pointerdown", selectFromScene);
  window.addEventListener("resize", resizeRenderer);
}

function initScene() {
  scene = new THREE.Scene();
  scene.background = new THREE.Color(0x101216);
  camera = new THREE.PerspectiveCamera(45, 1, 0.01, 10000);
  camera.position.set(3.0, -4.5, 3.0);
  renderer = new THREE.WebGLRenderer({ canvas: els.scene, antialias: true, alpha: false });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  controls.dampingFactor = 0.08;
  controls.screenSpacePanning = true;
  raycaster = new THREE.Raycaster();
  pointer = new THREE.Vector2();
  blockGroup = new THREE.Group();
  boundaryFaceGroup = new THREE.Group();
  selectedGridGroup = new THREE.Group();
  spacingJumpGroup = new THREE.Group();
  diagnosticGroup = new THREE.Group();
  partitionGroup = new THREE.Group();
  axesGroup = new THREE.Group();
  gridHelper = new THREE.GridHelper(10, 20, 0x46515c, 0x252b32);
  gridHelper.rotation.x = Math.PI / 2;
  scene.add(gridHelper);
  scene.add(axesGroup);
  scene.add(blockGroup);
  scene.add(boundaryFaceGroup);
  scene.add(selectedGridGroup);
  scene.add(spacingJumpGroup);
  scene.add(diagnosticGroup);
  scene.add(partitionGroup);

  transformControls = new TransformControls(camera, renderer.domElement);
  transformControls.setSpace("world");
  transformControls.setSize(0.82);
  transformControls.addEventListener("dragging-changed", handleTransformDragging);
  transformControls.addEventListener("objectChange", handleTransformObjectChange);
  scene.add(transformControls.getHelper ? transformControls.getHelper() : transformControls);

  const ambient = new THREE.HemisphereLight(0xffffff, 0x1a2028, 1.9);
  const key = new THREE.DirectionalLight(0xffffff, 1.4);
  key.position.set(4, -5, 6);
  scene.add(ambient, key);
  resizeRenderer();
}

function renderAll() {
  if (!project) return;
  project.schemaVersion = project.schemaVersion ?? 2;
  project.blocks = project.blocks || [];
  project.nscal = normalizedScalarCount(project.nscal ?? project.nScal);
  project.periodicAxes = normalizedPeriodicAxes(project);
  project.decomposition = normalizedDecomposition(project);
  for (const block of project.blocks) {
    ensureBoundaryArrays(block, project.nscal);
    block.axisLocks = normalizedAxisLocks(block.axisLocks);
  }
  project.blocks.sort((a, b) => a.id - b.id);
  if (!selectedBlock()) selectedId = project.blocks[0]?.id ?? null;
  selectedIds = new Set([...selectedIds].filter((id) => project.blocks.some((block) => block.id === id)));
  els.projectName.value = project.name ?? "snac-grid";
  els.scalarCount.value = project.nscal;
  els.inferConnectivity.checked = project.inferConnectivity !== false;
  project.writeExternalGrid = project.writeExternalGrid !== false;
  els.writeExternalGrid.checked = project.writeExternalGrid;
  if (!["grid", "both"].includes(project.externalGridSource)) project.externalGridSource = "grid";
  els.gridSource.value = project.externalGridSource;
  els.gridSource.disabled = !project.writeExternalGrid;
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
  if (selectedId != null && !previewCache.has(selectedId)) schedulePreview(selectedId);
  if (window.lucide) window.lucide.createIcons();
}

function renderBlockList() {
  els.blockList.innerHTML = "";
  if (!project.blocks.length) {
    els.blockList.innerHTML = `<div class="empty-state">No blocks</div>`;
    return;
  }
  for (const block of project.blocks) {
    const item = document.createElement("div");
    item.className = `block-item${block.id === selectedId ? " selected" : ""}`;
    item.innerHTML = `
      <input type="checkbox" aria-label="Select ${escapeHtml(blockLabel(block))}" ${selectedIds.has(block.id) ? "checked" : ""} />
      <span class="swatch" style="background:#${COLORS[(block.id - 1) % COLORS.length].toString(16).padStart(6, "0")}"></span>
      <button type="button" class="block-main"><strong>${escapeHtml(blockLabel(block))}</strong><span>${block.lmin.join(", ")} -> ${block.lmax.join(", ")}</span></button>
    `;
    item.querySelector("input").addEventListener("change", (event) => {
      if (event.target.checked) selectedIds.add(block.id);
      else selectedIds.delete(block.id);
      renderBlockList();
    });
    item.querySelector("button").addEventListener("click", () => {
      selectedId = block.id;
      renderAll();
    });
    els.blockList.appendChild(item);
  }
}

function renderInspector() {
  const block = selectedBlock();
  document.querySelector(".inspector").classList.toggle("empty", !block);
  if (!block) {
    els.blockTitle.textContent = "No block";
    els.blockName.value = "";
    els.inivel.innerHTML = "";
    els.periodicGrid.innerHTML = "";
    els.geometryGrid.innerHTML = "";
    els.resolutionGrid.innerHTML = "";
    els.axisEditor.innerHTML = "";
    els.boundaryTable.innerHTML = "";
    els.removeBlock.disabled = true;
    return;
  }
  els.removeBlock.disabled = false;
  els.blockTitle.textContent = blockLabel(block);
  els.blockName.value = block.name ?? "";
  els.blockName.oninput = () => {
    block.name = els.blockName.value;
    markProjectDirty();
    renderBlockList();
    renderScene();
  };

  els.inivel.innerHTML = INIT_FIELDS.map((value) => `<option value="${value}">${value}</option>`).join("");
  els.inivel.value = block.inivel ?? "zer";
  els.inivel.onchange = () => {
    block.inivel = els.inivel.value;
    markProjectDirty();
  };

  renderPeriodicInputs();
  renderGeometryInputs(block);
  renderAxisTabs();
  renderAxisEditor(block);
  renderBoundaryTable(block);
}

function renderPeriodicInputs() {
  els.periodicGrid.innerHTML = `
    <span>Periodic</span>
    ${AXES.map(
      (axis, index) => `
        <label class="check periodic-check">
          <input type="checkbox" data-axis-index="${index}" ${project.periodicAxes[index] ? "checked" : ""} />
          <span>${axis.toUpperCase()}</span>
        </label>
      `
    ).join("")}
  `;
  for (const input of els.periodicGrid.querySelectorAll("input")) {
    input.onchange = () => {
      project.periodicAxes[Number(input.dataset.axisIndex)] = input.checked;
      updateStructure();
    };
  }
}

function renderGeometryInputs(block) {
  els.geometryGrid.innerHTML = `
    <span></span>${AXES.map((axis) => `<span class="axis-label">${axis.toUpperCase()}</span>`).join("")}
    <span class="axis-label">min</span>${AXES.map((axis, index) => numberInput("lmin", index, block.lmin[index])).join("")}
    <span class="axis-label">max</span>${AXES.map((axis, index) => numberInput("lmax", index, block.lmax[index])).join("")}
  `;
  els.resolutionGrid.innerHTML = `
    <span></span>${AXES.map((axis) => `<span class="axis-label">${axis.toUpperCase()}</span>`).join("")}
    <span class="axis-label">ng</span>${AXES.map((axis, index) =>
      numberInput("ng", index, block.ng[index], 1, isCellCountDerived(ensureAxis(block, axis)))
    ).join("")}
    <span class="axis-label">mpi</span>${AXES.map((axis, index) => numberInput("dims", index, block.dims[index], 1)).join("")}
  `;
  for (const input of [...els.geometryGrid.querySelectorAll("input"), ...els.resolutionGrid.querySelectorAll("input")]) {
    input.addEventListener("input", () => {
      const field = input.dataset.field;
      const index = Number(input.dataset.index);
      const value = field === "ng" || field === "dims" ? Math.max(1, Math.round(Number(input.value) || 1)) : Number(input.value);
      const oldMin = block.lmin[index];
      const oldMax = block.lmax[index];
      block[field][index] = value;
      if (field === "lmin" || field === "lmax") {
        rescaleExplicitFaces(block, index, oldMin, oldMax, block.lmin[index], block.lmax[index]);
      }
      if (field === "dims") propagateMpiPartition(project.blocks, block.id, index, value);
      markProjectDirty(field === "ng" || field === "lmin" || field === "lmax" ? block.id : null);
      if (field === "ng" || field === "lmin" || field === "lmax") schedulePreview(block.id);
      renderScene();
      renderBlockList();
    });
    if (input.dataset.field === "dims") {
      input.addEventListener("change", () => updateStructure({ silent: true }));
    }
  }
}

function renderAxisTabs() {
  for (const button of els.axisTabs.querySelectorAll("button")) {
    button.classList.toggle("active", button.dataset.axis === currentAxis);
    button.onclick = () => {
      currentAxis = button.dataset.axis;
      renderInspector();
      renderSpacingPreview();
      renderScene();
    };
  }
}

function renderAxisEditor(block) {
  const axis = ensureAxis(block, currentAxis);
  const axisIndex = AXES.indexOf(currentAxis);
  const locked = block.axisLocks[axisIndex];
  const selectedKind = axis.kind ?? "snac";
  els.axisEditor.innerHTML = `
    <label class="field">
      <span>Mode</span>
      <select id="grid-kind">
        <option value="snac">Native SNaC</option>
        <option value="simple_ratio">Monotone grading</option>
        <option value="max_min">Symmetric grading</option>
        <option value="multi">Multi-region grading</option>
        ${selectedKind === "explicit" ? '<option value="explicit">Explicit coordinates</option>' : ""}
      </select>
    </label>
    <div id="grid-kind-body"></div>
    <div class="axis-actions">
      <button id="axis-lock" class="button icon${locked ? " active" : ""}" type="button" title="${locked ? "Unlock" : "Lock"} ${currentAxis.toUpperCase()} grid"><i data-lucide="${locked ? "lock" : "unlock"}"></i></button>
      <button id="apply-axis" class="button" type="button"><i data-lucide="copy-check"></i><span>Apply to aligned blocks</span></button>
    </div>
  `;
  const kindSelect = document.getElementById("grid-kind");
  kindSelect.value = selectedKind;
  kindSelect.onchange = () => {
    axis.kind = kindSelect.value;
    if (axis.kind === "simple_ratio") {
      axis.controls = axis.controls?.length === 2 ? axis.controls : ["n", "ratio"];
      setMonotoneControlsFromPreview(block, axis, axis.controls);
    }
    if (axis.kind === "max_min" && !["both", "middle"].includes(axis.side)) axis.side = "both";
    renderAxisEditor(block);
    renderGeometryInputs(block);
    gridDefinitionChanged(block);
  };
  renderAxisKindBody(block, axis);
  document.getElementById("axis-lock").onclick = () => {
    block.axisLocks[axisIndex] = !block.axisLocks[axisIndex];
    markProjectDirty();
    renderAxisEditor(block);
    if (window.lucide) window.lucide.createIcons();
  };
  document.getElementById("apply-axis").onclick = applyAxisToAlignedBlocks;
}

function profileSelect(axis) {
  const value = axis.profile ?? "geometric";
  return `
    <label class="field">
      <span>Profile</span>
      <select id="grid-profile">
        ${PROFILE_OPTIONS.map(([key, label]) => `<option value="${key}" ${key === value ? "selected" : ""}>${label}</option>`).join("")}
      </select>
    </label>
  `;
}

function bindProfileSelect(block, axis) {
  const select = document.getElementById("grid-profile");
  if (!select) return;
  select.onchange = () => {
    axis.profile = select.value;
    if (axis.kind === "simple_ratio") {
      if (axis.profile !== "geometric" && axis.controls?.includes("cell_ratio")) {
        setMonotoneControlsFromPreview(block, axis, ["n", "ratio"]);
      }
      renderAxisKindBody(block, axis);
    }
    gridDefinitionChanged(block);
  };
}

function renderAxisKindBody(block, axis) {
  const body = document.getElementById("grid-kind-body");
  const preview = previewCache.get(block.id)?.[currentAxis];
  if (axis.kind === "explicit") {
    const faces = axis.faces ?? [];
    const spacings = faces.slice(1).map((face, index) => face - faces[index]);
    const minimum = spacings.reduce((value, spacing) => Math.min(value, spacing), Number.POSITIVE_INFINITY);
    const maximum = spacings.reduce((value, spacing) => Math.max(value, spacing), Number.NEGATIVE_INFINITY);
    body.innerHTML = `
      <div class="explicit-grid-summary">
        <span>cells <strong>${Math.max(0, faces.length - 1)}</strong></span>
        <span>minimum <strong>${formatNumber(minimum)}</strong></span>
        <span>maximum <strong>${formatNumber(maximum)}</strong></span>
      </div>
    `;
  } else if (axis.kind === "multi") {
    axis.segments = axis.segments?.length ? axis.segments : [{ length: 1, cells: 1, ratio: 1, continuous: false }];
    body.innerHTML = `
      ${profileSelect(axis)}
      <div class="segment-table">
        ${axis.segments
          .map(
            (segment, index) => `
          <div class="segment-row" data-segment="${index}">
            <label class="field"><span>Length weight</span><input data-key="length" type="number" step="any" min="1e-12" value="${segment.length}"></label>
            <label class="field"><span>Cell weight</span><input data-key="cells" type="number" step="any" min="1e-12" value="${segment.cells}"></label>
            <label class="field"><span>End/start</span><input data-key="ratio" type="number" step="any" min="1e-12" value="${segment.ratio}" ${index > 0 && segment.continuous ? "disabled" : ""}></label>
            <label class="segment-link" title="Match the previous segment spacing"><input data-key="continuous" type="checkbox" ${index > 0 && segment.continuous ? "checked" : ""} ${index === 0 ? "disabled" : ""}><i data-lucide="link"></i></label>
            <button class="button icon segment-remove" title="Remove segment"><i data-lucide="minus"></i></button>
          </div>
          ${segmentMetricMarkup(preview?.segments?.[index], index)}
        `
          )
          .join("")}
      </div>
      <button id="add-segment" class="button"><i data-lucide="list-plus"></i><span>Add Segment</span></button>
    `;
    bindProfileSelect(block, axis);
    for (const row of body.querySelectorAll(".segment-row")) {
      const index = Number(row.dataset.segment);
      for (const input of row.querySelectorAll('input[type="number"]')) {
        input.oninput = () => {
          axis.segments[index][input.dataset.key] = positiveNumber(input.value, 1);
          gridDefinitionChanged(block);
        };
      }
      const continuity = row.querySelector('input[data-key="continuous"]');
      continuity.onchange = () => {
        axis.segments[index].continuous = continuity.checked;
        renderAxisKindBody(block, axis);
        gridDefinitionChanged(block);
      };
      row.querySelector(".segment-remove").onclick = () => {
        if (axis.segments.length > 1) axis.segments.splice(index, 1);
        renderAxisKindBody(block, axis);
        gridDefinitionChanged(block);
      };
    }
    document.getElementById("add-segment").onclick = () => {
      axis.segments.push({ length: 1, cells: 1, ratio: 1, continuous: false });
      renderAxisKindBody(block, axis);
      gridDefinitionChanged(block);
    };
  } else if (axis.kind === "simple_ratio") {
    normalizeMonotoneAxis(axis);
    const pairs = MONOTONE_CONTROL_PAIRS.filter((pair) => axis.profile === "geometric" || !pair.includes("cell_ratio"));
    const pairValue = axis.controls.join(",");
    body.innerHTML = `
      ${profileSelect(axis)}
      <label class="field"><span>Controls</span><select id="monotone-controls">
        ${pairs.map((pair) => `<option value="${pair.join(",")}" ${pair.join(",") === pairValue ? "selected" : ""}>${pair.map((key) => MONOTONE_CONTROL_LABELS[key]).join(" + ")}</option>`).join("")}
      </select></label>
      <div class="form-grid two">${axis.controls.map((key) => monotoneControlInput(block, axis, key)).join("")}</div>
    `;
    bindProfileSelect(block, axis);
    document.getElementById("monotone-controls").onchange = (event) => {
      setMonotoneControlsFromPreview(block, axis, event.target.value.split(","));
      renderAxisKindBody(block, axis);
      renderGeometryInputs(block);
      gridDefinitionChanged(block);
    };
    for (const input of body.querySelectorAll("input[data-control]")) {
      input.oninput = () => updateMonotoneControl(block, axis, input);
    }
  } else if (axis.kind === "max_min") {
    body.innerHTML = `
      ${profileSelect(axis)}
      <div class="form-grid two">
        <label class="field"><span>Minimum spacing</span><input id="min-dx" type="number" step="any" min="1e-12" value="${axis.min ?? 0.01}"></label>
        <label class="field"><span>Maximum spacing</span><input id="max-dx" type="number" step="any" min="1e-12" value="${axis.max ?? 0.02}"></label>
      </div>
      <label class="field"><span>Cluster</span>
        <select id="max-min-side">
          <option value="both">Both ends</option>
          <option value="middle">Middle</option>
        </select>
      </label>
    `;
    bindProfileSelect(block, axis);
    if (!["both", "middle"].includes(axis.side)) axis.side = "both";
    document.getElementById("max-min-side").value = axis.side;
    for (const id of ["min-dx", "max-dx"]) {
      document.getElementById(id).oninput = () => {
        axis[id === "min-dx" ? "min" : "max"] = positiveNumber(document.getElementById(id).value, 0.01);
        gridDefinitionChanged(block);
      };
    }
    document.getElementById("max-min-side").onchange = (event) => {
      axis.side = event.target.value;
      gridDefinitionChanged(block);
    };
  } else {
    body.innerHTML = `
      <label class="field">
        <span>Function</span>
        <select id="snac-gt">
          ${GRID_FUNCTIONS.map((label, index) => `<option value="${index}">${index} ${label}</option>`).join("")}
        </select>
      </label>
      <label class="field"><span>gr</span><input id="snac-gr" type="number" step="any" min="0" value="${axis.gr ?? 0}"></label>
    `;
    document.getElementById("snac-gt").value = axis.gt ?? 0;
    document.getElementById("snac-gt").onchange = (event) => {
      axis.gt = Number(event.target.value);
      gridDefinitionChanged(block);
    };
    document.getElementById("snac-gr").oninput = (event) => {
      axis.gr = Math.max(0, Number(event.target.value) || 0);
      gridDefinitionChanged(block);
    };
  }
  if (window.lucide) window.lucide.createIcons();
}

function normalizeMonotoneAxis(axis) {
  axis.controls = Array.isArray(axis.controls) && axis.controls.length === 2 ? axis.controls : ["n", "ratio"];
  if (axis.profile !== "geometric" && axis.controls.includes("cell_ratio")) axis.controls = ["n", "ratio"];
  axis.ratio = positiveNumber(axis.ratio, 1);
  axis.cell_ratio = positiveNumber(axis.cell_ratio, 1);
  axis.width_start = positiveNumber(axis.width_start, 0.01);
  axis.width_end = positiveNumber(axis.width_end, 0.02);
  axis.side = "end";
}

function monotoneControlInput(block, axis, key) {
  const axisIndex = AXES.indexOf(currentAxis);
  const value = key === "n" ? block.ng[axisIndex] : axis[key];
  const step = key === "n" ? "1" : "any";
  const min = key === "n" ? "1" : "1e-12";
  return `<label class="field"><span>${MONOTONE_CONTROL_LABELS[key]}</span><input data-control="${key}" type="number" step="${step}" min="${min}" value="${value}"></label>`;
}

function updateMonotoneControl(block, axis, input) {
  const key = input.dataset.control;
  if (key === "n") {
    const axisIndex = AXES.indexOf(currentAxis);
    block.ng[axisIndex] = Math.max(1, Math.round(Number(input.value) || 1));
    const overview = els.resolutionGrid.querySelector(`input[data-field="ng"][data-index="${axisIndex}"]`);
    if (overview) overview.value = block.ng[axisIndex];
  } else {
    axis[key] = positiveNumber(input.value, axis[key] ?? 1);
  }
  gridDefinitionChanged(block);
}

function setMonotoneControlsFromPreview(block, axis, controls) {
  const preview = previewCache.get(block.id)?.[currentAxis];
  axis.controls = controls;
  const axisIndex = AXES.indexOf(currentAxis);
  const values = {
    n: preview?.n ?? block.ng[axisIndex],
    ratio: preview?.expansion ?? axis.ratio ?? 1,
    cell_ratio: preview?.cellRatio ?? axis.cell_ratio ?? 1,
    width_start: preview?.lower ?? axis.width_start ?? 0.01,
    width_end: preview?.upper ?? axis.width_end ?? 0.02,
  };
  for (const key of controls) {
    if (key === "n") block.ng[axisIndex] = Math.max(1, Math.round(values.n));
    else axis[key] = positiveNumber(values[key], axis[key] ?? 1);
  }
  axis.side = "end";
}

function segmentMetricMarkup(metric, fallbackIndex = "") {
  const index = metric?.index ?? fallbackIndex;
  if (!metric) return `<div class="segment-metrics" data-segment-metric="${index}"></div>`;
  return `<div class="segment-metrics" data-segment-metric="${index}">
    <span>N <strong>${metric.cells}</strong></span>
    <span>lower <strong>${formatNumber(metric.widthStart)}</strong></span>
    <span>upper <strong>${formatNumber(metric.widthEnd)}</strong></span>
    <span>jump <strong>${formatNumber(metric.jump)}</strong></span>
  </div>`;
}

function updateSegmentMetrics(preview) {
  for (const metric of preview?.segments ?? []) {
    const container = els.axisEditor.querySelector(`[data-segment-metric="${metric.index}"]`);
    if (container) container.outerHTML = segmentMetricMarkup(metric);
    if (metric.continuous) {
      const ratio = els.axisEditor.querySelector(`.segment-row[data-segment="${metric.index}"] input[data-key="ratio"]`);
      if (ratio) ratio.value = formatNumber(metric.ratio);
    }
  }
}

async function applyAxisToAlignedBlocks() {
  const block = selectedBlock();
  if (!block) return;
  setStatus(`Applying ${currentAxis.toUpperCase()} grid...`);
  try {
    const payload = await postJson("/api/apply-axis", {
      project,
      sourceBlockId: block.id,
      axis: currentAxis,
    });
    project = payload.project;
    previewCache.clear();
    previewSequence += 1;
    lastCheck = null;
    decompositionResult = null;
    renderAll();
    scheduleHistoryCommit();
    const count = payload.changedBlockIds?.length ?? 0;
    setStatus(count ? `Applied ${currentAxis.toUpperCase()} grid to ${count} aligned block${count === 1 ? "" : "s"}` : "No other aligned blocks");
  } catch (error) {
    setStatus(error.message);
  }
}

function renderDecompositionControls() {
  if (!project?.decomposition) return;
  const spec = project.decomposition;
  els.targetRanks.value = spec.targetRanks;
  els.minLocalCells.value = spec.minLocalCells;
  els.maxLocalAspect.value = spec.maxLocalAspect;
  for (const button of els.decompositionMode.querySelectorAll("button")) {
    button.classList.toggle("active", button.dataset.mode === spec.mode);
  }
  for (const input of els.decompositionAxes.querySelectorAll("input")) {
    input.checked = spec.axes[Number(input.dataset.axis)];
  }
  els.balanceDecomposition.disabled = !project.blocks.length || spec.targetRanks < project.blocks.length;
  renderDecompositionSummary();
}

function renderDecompositionSummary() {
  if (!decompositionResult) {
    const target = project.decomposition.targetRanks;
    els.decompositionSummary.textContent = target > 0 ? "Balance to apply an exact decomposition." : "Set total ranks to enable balancing.";
    return;
  }
  const result = decompositionResult;
  if (!result.ok) {
    const nearby = result.nearbyRanks?.length
      ? `<div class="nearby-ranks"><span>Nearby feasible totals</span><div>${result.nearbyRanks
          .map((ranks) => `<button type="button" class="nearby-rank" data-ranks="${ranks}">${formatInteger(ranks)}</button>`)
          .join("")}</div></div>`
      : "";
    els.decompositionSummary.innerHTML = `${(result.errors ?? []).map((message) => `<div class="error-text">${escapeHtml(message)}</div>`).join("")}${nearby}`;
    for (const button of els.decompositionSummary.querySelectorAll(".nearby-rank")) {
      button.addEventListener("click", () => {
        project.decomposition.targetRanks = Number(button.dataset.ranks);
        markProjectDirty();
        renderDecompositionControls();
        void balanceDecomposition();
      });
    }
    return;
  }
  const axes = result.activeAxes?.length ? result.activeAxes.map((axis) => axis.toUpperCase()).join("") : "serial";
  const quality = result.optimal ? "best found" : "bounded search";
  els.decompositionSummary.innerHTML = `
    <div><strong>${formatInteger(result.totalRanks)}</strong> ranks · ${axes} · ${quality}</div>
    <div>${formatInteger(result.minCells)}-${formatInteger(result.maxCells)} cells/rank · imbalance <strong>${formatNumber(result.imbalance)}</strong></div>
    <div>minimum local edge <strong>${formatInteger(result.minLocalExtent)}</strong> cells · maximum aspect <strong>${formatNumber(result.maxLocalAspect)}</strong></div>
    ${(result.blocks ?? [])
      .map((block) => `<div>B${block.blockId}: ${block.dims.join(" × ")} · ${block.ranks} ranks</div>`)
      .join("")}
  `;
}

function renderQualitySummary() {
  const quality = lastCheck?.quality;
  const storage = lastCheck?.storage;
  const disabled = !project.blocks.length;
  els.downloadReportJson.disabled = disabled;
  els.downloadReportMarkdown.disabled = disabled;
  if (!quality) {
    els.qualitySummary.innerHTML = `<div class="quality-empty">${escapeHtml(lastCheck?.qualityError ?? "Not evaluated")}</div>`;
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

async function previewGridRepair() {
  if (!project.blocks.length) {
    setStatus("Add a block before repairing grids");
    return;
  }
  setStatus("Checking grid congruence...");
  try {
    const payload = await postJson(
      "/api/repair",
      { project, sourceBlockId: selectedId },
      { allowInvalid: true }
    );
    if (!payload.ok) {
      pendingRepair = null;
      lastCheck = payload;
      renderCheckResults();
      setStatus(`${payload.errors.length} repair issue${payload.errors.length === 1 ? "" : "s"}`);
      return;
    }
    if (!payload.changes?.length) {
      pendingRepair = null;
      setStatus("Grids and interface spacing are already congruent");
      return;
    }
    pendingRepair = payload;
    els.repairSummary.innerHTML = `${payload.changes
      .map((change) => {
        if (change.kind === "spacing") {
          return `<div class="repair-change"><strong>Block ${change.blockId} · ${change.face.toUpperCase()}</strong><br>${formatNumber(change.oldSpacing)} → ${formatNumber(change.newSpacing)}, matched to block ${change.sourceBlockId}</div>`;
        }
        return `<div class="repair-change"><strong>Block ${change.blockId} · ${change.axis.toUpperCase()}</strong><br>${change.oldN} → ${change.newN} cells, copied from block ${change.sourceBlockId}</div>`;
      })
      .join("")}${(payload.warnings ?? []).map((warning) => `<div class="repair-warning">${escapeHtml(warning)}</div>`).join("")}`;
    els.repairDialog.showModal();
    setStatus(`${payload.changes.length} proposed grid repair${payload.changes.length === 1 ? "" : "s"}`);
  } catch (error) {
    setStatus(error.message);
  }
}

function closeRepairDialog() {
  pendingRepair = null;
  if (els.repairDialog.open) els.repairDialog.close();
}

async function applyGridRepair() {
  if (!pendingRepair?.project) return;
  const count = pendingRepair.changes?.length ?? 0;
  project = pendingRepair.project;
  pendingRepair = null;
  decompositionResult = null;
  if (els.repairDialog.open) els.repairDialog.close();
  previewCache.clear();
  previewSequence += 1;
  lastCheck = null;
  renderAll();
  scheduleHistoryCommit();
  const options = { successMessage: `Applied ${count} grid repair${count === 1 ? "" : "s"}` };
  if (project.inferConnectivity) await updateStructure(options);
  else {
    const valid = await checkProject({ silent: true });
    setStatus(valid ? options.successMessage : "Grid repair applied; checks still need attention");
  }
}

async function balanceDecomposition() {
  setStatus("Balancing MPI decomposition...");
  try {
    const payload = await postJson("/api/decompose", { project }, { allowInvalid: true });
    decompositionResult = payload;
    if (!payload.ok) {
      renderDecompositionSummary();
      scheduleHistoryCommit();
      setStatus(payload.errors?.[0] ?? "No exact decomposition found");
      return;
    }
    project = payload.project;
    previewCache.clear();
    previewSequence += 1;
    lastCheck = null;
    renderAll();
    scheduleHistoryCommit();
    const warningText = payload.warnings?.length ? `, ${payload.warnings.length} warning${payload.warnings.length === 1 ? "" : "s"}` : "";
    setStatus(`Balanced ${payload.totalRanks} MPI ranks${warningText}`);
    if (project.inferConnectivity) await updateStructure({ silent: true });
    else await checkProject({ silent: true });
  } catch (error) {
    setStatus(error.message);
  }
}

async function applyBoundaryPreset() {
  const blockIds = selectedBlockIds();
  if (!blockIds.length) return;
  setStatus("Applying boundary preset...");
  try {
    const payload = await postJson("/api/bc-preset", {
      project,
      blockIds,
      preset: els.bcPreset.value,
      face: els.bcPresetFace.value,
      velocity: [els.bcPresetU.value, els.bcPresetV.value, els.bcPresetW.value].map(Number),
      pressure: Number(els.bcPresetP.value) || 0,
    });
    project = payload.project;
    markProjectDirty();
    renderAll();
    setStatus(`Applied preset to ${blockIds.length} block${blockIds.length === 1 ? "" : "s"}`);
  } catch (error) {
    setStatus(error.message);
  }
}

function renderBoundaryTable(block) {
  ensureBoundaryArrays(block, project.nscal);
  els.boundaryTable.innerHTML = FACE_ORDER.map((face, index) => {
    return `
      <div class="boundary-row">
        <span class="axis-label">${face.toUpperCase()}</span>
        <div class="boundary-components">
          ${["u", "v", "w"]
            .map((name, component) => boundaryComponent(name, index, block.cbcvel[component][index], block.bcvel[component][index], component))
            .join("")}
          ${boundaryComponent("p", index, block.cbcpre[index], block.bcpre[index])}
          ${Array.from({ length: project.nscal }, (_, iscal) =>
            boundaryComponent(`s${iscal + 1}`, index, block.cbcscal[iscal][index], block.bcscal[iscal][index], iscal, "scal")
          ).join("")}
        </div>
      </div>
    `;
  }).join("");
  for (const select of els.boundaryTable.querySelectorAll("select")) {
    select.onchange = () => {
      const index = Number(select.dataset.face);
      if (select.dataset.kind === "vel") {
        block.cbcvel[Number(select.dataset.component)][index] = select.value;
      } else if (select.dataset.kind === "scal") {
        block.cbcscal[Number(select.dataset.component)][index] = select.value;
      } else {
        block.cbcpre[index] = select.value;
      }
      markProjectDirty();
      renderScene();
    };
  }
  for (const input of els.boundaryTable.querySelectorAll("input")) {
    input.oninput = () => {
      const face = Number(input.dataset.face);
      const value = Number(input.value);
      if (!Number.isFinite(value)) return;
      if (input.dataset.kind === "pre") {
        block.bcpre[face] = value;
      } else if (input.dataset.kind === "scal") {
        block.bcscal[Number(input.dataset.component)][face] = value;
      } else {
        block.bcvel[Number(input.dataset.component)][face] = value;
      }
      markProjectDirty();
    };
  }
}

function renderCheckResults() {
  if (!els.checkResults) return;
  const errors = lastCheck?.errors ?? [];
  const warnings = lastCheck?.warnings ?? [];
  const diagnostics = lastCheck?.diagnostics ?? [
    ...errors.map((message) => ({ severity: "error", message, locations: [] })),
    ...warnings.map((message) => ({ severity: "warning", message, locations: [] })),
  ];
  els.exportCase.disabled = !project.blocks.length || errors.length > 0;
  if (!lastCheck) {
    activeDiagnostic = null;
    els.checkResults.innerHTML = "";
    return;
  }
  if (
    activeDiagnostic &&
    !diagnostics.some((item) => item.code === activeDiagnostic.code && item.message === activeDiagnostic.message)
  ) {
    activeDiagnostic = null;
  }
  if (!errors.length && !warnings.length) {
    activeDiagnostic = null;
    els.checkResults.innerHTML = `<div class="check-line ok">Structured grid checks passed</div>`;
    return;
  }
  els.checkResults.innerHTML = diagnostics
    .map((diagnostic, index) => {
      const actionable = (diagnostic.locations ?? []).some((location) =>
        project.blocks.some((block) => block.id === location.blockId)
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
}

function focusDiagnostic(diagnostic) {
  const blockIds = (diagnostic.locations ?? [])
    .map((location) => location.blockId)
    .filter((blockId) => project.blocks.some((block) => block.id === blockId));
  if (!blockIds.length) return;
  activeDiagnostic = diagnostic;
  selectedId = blockIds[0];
  selectedIds = new Set(blockIds);
  if (AXES.includes(diagnostic.axis)) currentAxis = diagnostic.axis;
  const block = selectedBlock();
  if (block) controls.target.set(...blockCenter(block));
  renderAll();
  setStatus(`Showing ${diagnostic.code?.replaceAll("_", " ") ?? "grid issue"}`);
}

function renderScene() {
  transformControls.detach();
  selectedBlockMesh = null;
  boundaryFaceGroup.visible = true;
  selectedGridGroup.visible = gridLinesVisible;
  partitionGroup.visible = true;
  spacingJumpGroup.visible = spacingJumpsVisible;
  clearDisposableGroup(blockGroup);
  clearDisposableGroup(boundaryFaceGroup);
  clearDisposableGroup(selectedGridGroup);
  clearDisposableGroup(partitionGroup);
  clearDisposableGroup(spacingJumpGroup);
  clearDisposableGroup(diagnosticGroup);
  const diagnosticIds = new Set(
    (activeDiagnostic?.locations ?? []).map((location) => location.blockId)
  );
  for (const block of project.blocks) {
    const size = blockSize(block);
    const center = blockCenter(block);
    const color = COLORS[(block.id - 1) % COLORS.length];
    const material = new THREE.MeshStandardMaterial({
      color,
      transparent: true,
      opacity: diagnosticIds.has(block.id) ? 0.52 : block.id === selectedId ? 0.43 : 0.24,
      roughness: 0.62,
      metalness: 0.0,
      side: THREE.DoubleSide,
    });
    const mesh = new THREE.Mesh(new THREE.BoxGeometry(size[0], size[1], size[2]), material);
    mesh.position.set(center[0], center[1], center[2]);
    mesh.userData.blockId = block.id;
    mesh.userData.baseSize = size;
    if (block.id === selectedId) selectedBlockMesh = mesh;
    blockGroup.add(mesh);

    const edges = new THREE.LineSegments(
      new THREE.EdgesGeometry(mesh.geometry),
      new THREE.LineBasicMaterial({
        color: diagnosticIds.has(block.id)
          ? activeDiagnostic.severity === "error" ? 0xff6b6b : 0xf0c36a
          : block.id === selectedId ? 0xffffff : 0x9aa6b2,
        transparent: true,
        opacity: diagnosticIds.has(block.id) || block.id === selectedId ? 0.95 : 0.48,
      })
    );
    mesh.add(edges);
  }
  drawBoundaryFaces();
  drawSelectedGrid();
  drawPartitionPlanes();
  drawInterfaceSpacing();
  drawDiagnosticHighlight();
  updateAxesHelper();
  attachSelectedTransform();
  renderViewControls();
}

function drawDiagnosticHighlight() {
  if (!activeDiagnostic) return;
  const color = activeDiagnostic.severity === "error" ? 0xff6b6b : 0xf0c36a;
  for (const location of activeDiagnostic.locations ?? []) {
    if (!location.face) continue;
    const block = project.blocks.find((item) => item.id === location.blockId);
    const face = FACE_ORDER.indexOf(location.face);
    if (!block || face < 0) continue;
    const axisIndex = Math.floor(face / 2);
    const value = face % 2 === 0 ? block.lmin[axisIndex] : block.lmax[axisIndex];
    const line = new THREE.Line(
      new THREE.BufferGeometry().setFromPoints(planeOutlinePoints(block, axisIndex, value)),
      new THREE.LineBasicMaterial({ color, depthTest: false })
    );
    line.renderOrder = 8;
    diagnosticGroup.add(line);
  }
}

function drawInterfaceSpacing() {
  if (!spacingJumpsVisible) return;
  for (const metric of lastCheck?.interfaces ?? []) {
    const block = project.blocks.find((item) => item.id === metric.blockId);
    if (!block) continue;
    const axisIndex = AXES.indexOf(metric.axis);
    const value = metric.face.endsWith("+") ? block.lmax[axisIndex] : block.lmin[axisIndex];
    const color = metric.ratio > 3 ? 0xe46868 : metric.ratio > 1.5 ? 0xd6bf5b : 0x43c989;
    const line = new THREE.Line(
      new THREE.BufferGeometry().setFromPoints(planeOutlinePoints(block, axisIndex, value)),
      new THREE.LineBasicMaterial({ color, transparent: true, opacity: 0.95, depthTest: false })
    );
    line.renderOrder = 5;
    spacingJumpGroup.add(line);
  }
}

function drawBoundaryFaces() {
  const block = selectedBlock();
  if (!block || !boundaryColorsVisible) return;
  const size = blockSize(block);
  const center = blockCenter(block);
  for (let face = 0; face < 6; face += 1) {
    const axisIndex = Math.floor(face / 2);
    const side = face % 2;
    const dimensions = axisIndex === 0 ? [size[1], size[2]] : axisIndex === 1 ? [size[0], size[2]] : [size[0], size[1]];
    const geometry = new THREE.PlaneGeometry(dimensions[0], dimensions[1]);
    if (axisIndex === 0) geometry.rotateY(Math.PI / 2);
    if (axisIndex === 1) geometry.rotateX(Math.PI / 2);
    const material = new THREE.MeshBasicMaterial({
      color: boundaryFaceColor(block, face),
      transparent: true,
      opacity: 0.2,
      side: THREE.DoubleSide,
      depthWrite: false,
      polygonOffset: true,
      polygonOffsetFactor: -2,
    });
    const mesh = new THREE.Mesh(geometry, material);
    mesh.position.set(...center);
    mesh.position.setComponent(axisIndex, side === 0 ? block.lmin[axisIndex] : block.lmax[axisIndex]);
    mesh.renderOrder = 2;
    boundaryFaceGroup.add(mesh);
  }
}

function boundaryFaceColor(block, face) {
  const velocityCodes = block.cbcvel.map((component) => component[face]);
  const allCodes = [block.cbcpre[face], ...velocityCodes, ...block.cbcscal.map((scalar) => scalar[face])];
  if (allCodes.every((code) => code === "F")) return BOUNDARY_COLORS.F;
  const velocityCode = velocityCodes[0];
  if (["D", "N"].includes(velocityCode) && velocityCodes.every((code) => code === velocityCode)) {
    return BOUNDARY_COLORS[velocityCode];
  }
  return BOUNDARY_COLORS.mixed;
}

function drawPartitionPlanes() {
  const block = selectedBlock();
  const preview = block ? previewCache.get(block.id) : null;
  if (!block || !preview || !partitionPlanesVisible) return;
  const scale = Math.max(...blockSize(block), 1.0);
  const material = new THREE.LineDashedMaterial({
    color: 0xd6bf5b,
    transparent: true,
    opacity: 0.82,
    dashSize: scale * 0.025,
    gapSize: scale * 0.014,
  });
  for (let axisIndex = 0; axisIndex < 3; axisIndex += 1) {
    const faces = preview[AXES[axisIndex]]?.faces ?? [];
    for (const faceIndex of partitionFaceIndices(block.ng[axisIndex], block.dims[axisIndex])) {
      if (faceIndex <= 0 || faceIndex >= faces.length - 1) continue;
      const geometry = new THREE.BufferGeometry().setFromPoints(planeOutlinePoints(block, axisIndex, faces[faceIndex]));
      const line = new THREE.Line(geometry, material);
      line.computeLineDistances();
      line.renderOrder = 3;
      partitionGroup.add(line);
    }
  }
}

function updateAxesHelper() {
  clearDisposableGroup(axesGroup, false);
  const box = projectBox(project.blocks);
  const size = box.getSize(new THREE.Vector3());
  const origin = box.min.clone();
  if (!Number.isFinite(origin.x)) origin.set(0, 0, 0);
  const length = Math.max(size.x, size.y, size.z, 1) * 0.22;
  const headLength = length * 0.18;
  const headWidth = length * 0.08;
  const axes = [
    ["X", new THREE.Vector3(1, 0, 0), 0xe46868],
    ["Y", new THREE.Vector3(0, 1, 0), 0x43c989],
    ["Z", new THREE.Vector3(0, 0, 1), 0x76a9ff],
  ];
  for (const [label, direction, color] of axes) {
    const arrow = new THREE.ArrowHelper(direction, origin, length, color, headLength, headWidth);
    const sprite = axisLabel(label, color);
    sprite.position.copy(origin.clone().add(direction.clone().multiplyScalar(length * 1.16)));
    sprite.scale.setScalar(length * 0.16);
    axesGroup.add(arrow, sprite);
  }
}

function drawSelectedGrid() {
  const block = selectedBlock();
  if (!block || !gridLinesVisible) return;
  const preview = previewCache.get(block.id);
  if (!preview) return;
  const material = new THREE.LineBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.32 });
  for (const axis of AXES) {
    const idx = AXES.indexOf(axis);
    const faces = preview[axis]?.faces ?? [];
    const step = Math.max(1, Math.floor(faces.length / 34));
    for (let i = 0; i < faces.length; i += step) {
      const value = faces[i];
      const points = planeOutlinePoints(block, idx, value);
      selectedGridGroup.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(points), material));
    }
  }
}

function renderSpacingPreview() {
  const block = selectedBlock();
  const canvas = els.spacingPreview;
  const ctx = canvas.getContext("2d");
  const width = canvas.width;
  const height = canvas.height;
  ctx.clearRect(0, 0, width, height);
  ctx.fillStyle = "#11161b";
  ctx.fillRect(0, 0, width, height);
  if (!block) {
    els.spacingStats.innerHTML = "";
    return;
  }
  const preview = previewCache.get(block.id)?.[currentAxis];
  if (!preview?.faces?.length) {
    els.spacingStats.innerHTML = "";
    return;
  }
  const faces = preview.faces;
  const widths = faces.slice(1).map((face, index) => face - faces[index]);
  const min = preview.min;
  const max = preview.max;
  const span = faces[faces.length - 1] - faces[0];
  ctx.strokeStyle = "#2fa872";
  ctx.lineWidth = 1;
  for (const face of faces) {
    const x = 8 + ((face - faces[0]) / span) * (width - 16);
    ctx.beginPath();
    ctx.moveTo(x, 12);
    ctx.lineTo(x, height - 12);
    ctx.stroke();
  }
  ctx.strokeStyle = "#d28a3e";
  ctx.beginPath();
  widths.forEach((value, index) => {
    const x = 8 + (index / Math.max(1, widths.length - 1)) * (width - 16);
    const y = height - 12 - ((value - min) / Math.max(max - min, 1e-16)) * (height - 28);
    if (index === 0) ctx.moveTo(x, y);
    else ctx.lineTo(x, y);
  });
  ctx.stroke();
  const residuals = Object.values(preview.residuals ?? {}).map((value) => Math.abs(Number(value)));
  const fitError = residuals.length ? Math.max(...residuals) : 0;
  const idealCells = Number.isFinite(preview.idealN) ? formatNumber(preview.idealN) : "-";
  els.spacingStats.innerHTML = `
    <span>cells<strong>${preview.n ?? widths.length}</strong></span>
    <span>lower<strong>${formatNumber(preview.lower ?? widths[0])}</strong></span>
    <span>upper<strong>${formatNumber(preview.upper ?? widths[widths.length - 1])}</strong></span>
    <span>min<strong>${formatNumber(min)}</strong></span>
    <span>max<strong>${formatNumber(max)}</strong></span>
    <span>upper/lower<strong>${formatNumber(preview.expansion ?? 1)}</strong></span>
    <span>max growth<strong>${formatNumber(preview.maxGrowth ?? 1)}</strong></span>
    <span>ideal cells<strong>${idealCells}</strong></span>
    <span>fit error<strong>${formatPercent(fitError)}</strong></span>
  `;
}

function schedulePreview(blockId = selectedId, delay = 120) {
  if (blockId == null) return;
  window.clearTimeout(previewTimer);
  const sequence = ++previewSequence;
  previewTimer = window.setTimeout(() => requestPreview(blockId, sequence), delay);
}

async function requestPreview(blockId, sequence) {
  try {
    const payload = await postJson("/api/preview", { project, blockId });
    if (sequence !== previewSequence) return;
    const block = project.blocks.find((item) => item.id === payload.blockId);
    if (!block) return;
    const cellCountChanged = block.ng.some((value, index) => value !== payload.ng[index]);
    block.ng = payload.ng;
    if (cellCountChanged) scheduleHistoryCommit();
    previewCache.set(block.id, payload.axes);
    for (const input of els.resolutionGrid.querySelectorAll('input[data-field="ng"]')) {
      const index = Number(input.dataset.index);
      if (block.id === selectedId) input.value = block.ng[index];
    }
    const monotoneCells = els.axisEditor.querySelector('input[data-control="n"]');
    if (monotoneCells && block.id === selectedId) monotoneCells.value = block.ng[AXES.indexOf(currentAxis)];
    updateSegmentMetrics(payload.axes[currentAxis]);
    renderSpacingPreview();
    renderScene();
  } catch (error) {
    if (sequence !== previewSequence) return;
    previewCache.delete(blockId);
    renderSpacingPreview();
    renderScene();
    setStatus(error.message);
  }
}

function gridDefinitionChanged(block) {
  if (!block) return;
  markProjectDirty(block.id);
  schedulePreview(block.id);
  renderSpacingPreview();
  renderScene();
}

function setBlockTool(tool) {
  blockTool = tool;
  renderScene();
}

function attachSelectedTransform() {
  if (!selectedBlockMesh || blockTool === "select") return;
  transformControls.setMode(blockTool);
  transformControls.attach(selectedBlockMesh);
}

function handleTransformDragging(event) {
  if (event.value) {
    if (!selectedBlockMesh) return;
    commitHistory();
    isTransformDragging = true;
    transformDrag = {
      blockId: selectedBlockMesh.userData.blockId,
      baseSize: [...selectedBlockMesh.userData.baseSize],
      baseLmin: [...selectedBlock().lmin],
      baseLmax: [...selectedBlock().lmax],
      baseFaces: AXES.map((axis) => {
        const spec = ensureAxis(selectedBlock(), axis);
        return spec.kind === "explicit" ? [...spec.faces] : null;
      }),
    };
    controls.enabled = false;
    boundaryFaceGroup.visible = false;
    selectedGridGroup.visible = false;
    partitionGroup.visible = false;
    return;
  }

  controls.enabled = true;
  if (!isTransformDragging || !transformDrag) return;
  const blockId = transformDrag.blockId;
  isTransformDragging = false;
  transformDrag = null;
  commitHistory();
  window.setTimeout(() => {
    renderAll();
    schedulePreview(blockId, 0);
  }, 0);
}

function handleTransformObjectChange() {
  const mesh = transformControls.object;
  if (!isTransformDragging || !transformDrag || !mesh) return;
  const block = project.blocks.find((item) => item.id === transformDrag.blockId);
  if (!block) return;
  const center = mesh.position.toArray();
  const size = transformDrag.baseSize.map((value, index) => Math.max(1.0e-9, value * Math.abs(mesh.scale.getComponent(index))));
  block.lmin = center.map((value, index) => value - 0.5 * size[index]);
  block.lmax = center.map((value, index) => value + 0.5 * size[index]);
  for (let index = 0; index < 3; index += 1) {
    rescaleExplicitFaces(
      block,
      index,
      transformDrag.baseLmin[index],
      transformDrag.baseLmax[index],
      block.lmin[index],
      block.lmax[index],
      transformDrag.baseFaces[index]
    );
  }
  updateGeometryInputValues(block);
  markProjectDirty(block.id);
  renderBlockList();
}

function updateGeometryInputValues(block) {
  for (const field of ["lmin", "lmax"]) {
    for (const input of els.geometryGrid.querySelectorAll(`input[data-field="${field}"]`)) {
      input.value = formatInputNumber(block[field][Number(input.dataset.index)]);
    }
  }
}

function renderViewControls() {
  if (!selectedBlock() && blockTool !== "select") blockTool = "select";
  const toolButtons = {
    select: els.toolSelect,
    translate: els.toolMove,
    scale: els.toolScale,
  };
  for (const [tool, button] of Object.entries(toolButtons)) button.classList.toggle("active", blockTool === tool);
  els.toolMove.disabled = !selectedBlock();
  els.toolScale.disabled = !selectedBlock();
  els.toggleBoundaries.classList.toggle("active", boundaryColorsVisible);
  els.togglePartitions.classList.toggle("active", partitionPlanesVisible);
  els.toggleGridLines.classList.toggle("active", gridLinesVisible);
  els.toggleSpacingJumps.classList.toggle("active", spacingJumpsVisible);
}

function addFreeBlock() {
  const source = selectedBlock() ?? project.blocks[0];
  const copy = source ? structuredClone(source) : defaultBlock(nextBlockId(project.blocks));
  copy.id = nextBlockId(project.blocks);
  copy.name = `block-${copy.id}`;
  copy.axisLocks = [false, false, false];
  if (source) {
    const oldMin = copy.lmin[0];
    const oldMax = copy.lmax[0];
    const width = source.lmax[0] - source.lmin[0];
    copy.lmin[0] += width;
    copy.lmax[0] += width;
    rescaleExplicitFaces(copy, 0, oldMin, oldMax, copy.lmin[0], copy.lmax[0]);
  }
  resetFriendBoundaries(copy, project.nscal);
  project.blocks.push(copy);
  selectedId = copy.id;
  selectedIds = new Set([copy.id]);
  markProjectDirty(copy.id);
  renderAll();
}

function addAdjacentBlock() {
  const source = selectedBlock();
  if (!source) return;
  const face = els.adjacentFace.value;
  const axis = AXES.indexOf(face[0]);
  const sign = face[1] === "+" ? 1 : -1;
  const thickness = positiveNumber(els.adjacentSize.value, source.lmax[axis] - source.lmin[axis]);
  const copy = structuredClone(source);
  copy.id = nextBlockId(project.blocks);
  copy.name = `block-${copy.id}`;
  copy.axisLocks = [false, false, false];
  const oldMin = copy.lmin[axis];
  const oldMax = copy.lmax[axis];
  if (sign > 0) {
    copy.lmin[axis] = source.lmax[axis];
    copy.lmax[axis] = source.lmax[axis] + thickness;
  } else {
    copy.lmax[axis] = source.lmin[axis];
    copy.lmin[axis] = source.lmin[axis] - thickness;
  }
  rescaleExplicitFaces(copy, axis, oldMin, oldMax, copy.lmin[axis], copy.lmax[axis]);
  resetFriendBoundaries(copy, project.nscal);
  project.blocks.push(copy);
  selectedId = copy.id;
  selectedIds = new Set([copy.id]);
  markProjectDirty(copy.id);
  renderAll();
}

async function applyGeometryOperation(operation) {
  if (!selectedBlock()) return;
  let blockIds = selectedBlockIds();
  const payload = {
    project,
    operation,
    blockIds,
    sourceBlockId: selectedId,
    axis: els.geometryAxis.value,
  };
  if (operation === "duplicate") {
    payload.count = Math.max(1, Math.round(Number(els.arrayCount.value) || 1));
    payload.gap = Number(els.arrayGap.value) || 0;
  } else if (operation === "mirror") {
    payload.plane = Number(els.mirrorPlane.value) || 0;
  } else if (operation === "align") {
    payload.mode = els.alignMode.value;
  } else if (operation === "snap") {
    blockIds = blockIds.filter((id) => id !== selectedId);
    if (blockIds.length !== 1) {
      setStatus("Select one target block and make the source block primary");
      return;
    }
    payload.blockIds = blockIds;
    payload.face = els.snapFace.value;
  }
  setStatus(`${operation[0].toUpperCase()}${operation.slice(1)} blocks...`);
  try {
    const result = await postJson("/api/geometry", payload);
    project = result.project;
    const changed = result.changedBlockIds ?? [];
    if (["duplicate", "mirror"].includes(operation) && changed.length) {
      selectedIds = new Set(changed);
      selectedId = changed[0];
    }
    previewCache.clear();
    previewSequence += 1;
    markProjectDirty();
    renderAll();
    setStatus(`${operation[0].toUpperCase()}${operation.slice(1)} updated ${changed.length} block${changed.length === 1 ? "" : "s"}`);
  } catch (error) {
    setStatus(error.message);
  }
}

function removeSelectedBlock() {
  if (!selectedBlock()) return;
  project.blocks = project.blocks.filter((block) => block.id !== selectedId);
  previewCache.delete(selectedId);
  selectedIds.delete(selectedId);
  selectedId = project.blocks[0]?.id ?? null;
  markProjectDirty();
  renderAll();
}

function resetProject() {
  project.blocks = [];
  selectedId = null;
  selectedIds.clear();
  lastCheck = null;
  pendingRepair = null;
  decompositionResult = null;
  previewCache.clear();
  previewSequence += 1;
  renderAll();
  commitHistory();
  setStatus("All blocks erased");
}

async function checkProject(options = {}) {
  if (!options.silent) setStatus("Checking grid...");
  try {
    const payload = await postJson("/api/check", { project }, { allowInvalid: true });
    lastCheck = payload;
    renderCheckResults();
    renderQualitySummary();
    renderScene();
    if (!payload.ok) {
      setStatus(`${payload.errors.length} check issue${payload.errors.length === 1 ? "" : "s"}`);
      return false;
    }
    const warningText = payload.warnings?.length ? `, ${payload.warnings.length} warnings` : "";
    if (!options.silent) setStatus(`Grid checks passed${warningText}`);
    return true;
  } catch (error) {
    lastCheck = { ok: false, errors: [error.message], warnings: [] };
    renderCheckResults();
    renderQualitySummary();
    renderScene();
    setStatus(error.message);
    return false;
  }
}

async function updateStructure(options = {}) {
  if (!options.silent) setStatus("Updating structure...");
  try {
    const payload = await postJson("/api/update", { project, sourceBlockId: selectedId }, { allowInvalid: true });
    if (payload.project) {
      const oldSelectedId = selectedId;
      project = payload.project;
      selectedId = project.blocks.find((block) => block.id === oldSelectedId)?.id ?? project.blocks[0]?.id ?? null;
      previewCache.clear();
      previewSequence += 1;
    }
    lastCheck = payload;
    renderAll();
    scheduleHistoryCommit();
    const warningText = payload.warnings?.length ? `, ${payload.warnings.length} warnings` : "";
    if (!options.silent || !payload.ok) {
      const successMessage = options.successMessage ?? "Structure updated";
      setStatus(payload.ok ? `${successMessage}${warningText}` : `${payload.errors.length} check issue${payload.errors.length === 1 ? "" : "s"}`);
    }
    return payload.ok;
  } catch (error) {
    setStatus(error.message);
    return false;
  }
}

async function exportCase() {
  const valid = project.inferConnectivity ? await updateStructure({ silent: true }) : await checkProject({ silent: true });
  if (!valid) {
    setStatus("Fix grid checks before writing");
    return;
  }
  setStatus("Writing files...");
  try {
    const payload = await postJson("/api/export", { outputDir: els.outputDir.value, project });
    commitHistory();
    clearAutosave();
    const warningText = payload.warnings?.length ? `, ${payload.warnings.length} warnings` : "";
    setStatus(`Wrote ${payload.files.length} files to ${payload.outputDir}${warningText}`);
  } catch (error) {
    setStatus(error.message);
  }
}

async function importCase() {
  const caseDir = els.outputDir.value.trim();
  if (!caseDir) {
    setStatus("Set the case folder in Output");
    return;
  }
  setStatus("Importing SNaC case...");
  try {
    const payload = await postJson("/api/import", { caseDir });
    project = payload.project;
    selectedId = project.blocks[0]?.id ?? null;
    selectedIds = new Set(selectedId == null ? [] : [selectedId]);
    lastCheck = null;
    pendingRepair = null;
    decompositionResult = null;
    previewCache.clear();
    previewSequence += 1;
    renderAll();
    clearAutosave();
    resetHistory(false);
    fitView();
    const warningText = payload.warnings?.length
      ? `, ${payload.warnings.length} warning${payload.warnings.length === 1 ? "" : "s"}`
      : "";
    setStatus(`Imported ${project.blocks.length} blocks from ${payload.source}${warningText}`);
    await checkProject({ silent: true });
  } catch (error) {
    setStatus(error.message);
  }
}

function showLibraryDialog() {
  renderLibraryDialog();
  els.libraryDialog.showModal();
}

function renderLibraryDialog() {
  for (const button of els.libraryTabs.querySelectorAll("button")) {
    button.classList.toggle("active", button.dataset.libraryTab === currentLibraryTab);
  }
  const entries = libraryEntries();
  if (!entries.some((entry) => entry.key === selectedLibraryItem)) selectedLibraryItem = entries[0]?.key ?? null;
  els.libraryItems.innerHTML = entries
    .map((entry) => `<option value="${escapeHtml(entry.key)}">${escapeHtml(entry.item.name)}${entry.builtIn ? " (built-in)" : ""}</option>`)
    .join("");
  els.libraryItems.value = selectedLibraryItem ?? "";
  renderLibrarySelection();
  if (window.lucide) window.lucide.createIcons();
}

function renderLibrarySelection() {
  const entry = selectedLibraryEntry();
  const axisMode = currentLibraryTab === "axisPresets";
  if (!entry) {
    els.libraryName.value = "";
    els.librarySummary.textContent = "Empty library";
  } else {
    els.libraryName.value = entry.builtIn ? `${entry.item.name} copy` : entry.item.name;
    els.librarySummary.innerHTML = libraryItemSummary(entry);
  }
  els.applyLibraryItem.disabled = !entry || (axisMode && !selectedBlock());
  els.saveLibraryItem.disabled = axisMode && !selectedBlock();
  els.saveLibraryItem.title = axisMode ? "Save current axis to library" : "Save current project to library";
  els.renameLibraryItem.disabled = !entry || entry.builtIn;
  els.deleteLibraryItem.disabled = !entry || entry.builtIn;
}

function libraryEntries() {
  const builtIns = currentLibraryTab === "axisPresets" ? builtInAxisPresets() : builtInProjectTemplates();
  return [
    ...builtIns.map((item) => ({ key: `builtin:${item.id}`, item, builtIn: true })),
    ...userLibrary[currentLibraryTab].map((item) => ({ key: `user:${item.id}`, item, builtIn: false })),
  ];
}

function selectedLibraryEntry() {
  return libraryEntries().find((entry) => entry.key === selectedLibraryItem) ?? null;
}

function libraryItemSummary(entry) {
  const source = entry.builtIn ? "Built-in" : "User";
  if (currentLibraryTab === "axisPresets") {
    const axis = entry.item.axis;
    const rows = [
      `<strong>${escapeHtml(source)} axis preset</strong>`,
      `${formatInteger(entry.item.ng)} cells`,
      escapeHtml(axisKindLabel(axis.kind)),
    ];
    if (axis.profile) rows.push(escapeHtml(axis.profile));
    if (axis.kind === "multi") rows.push(`${axis.segments?.length ?? 0} regions`);
    if (entry.item.coordinateSpace === "relative") rows.push("Relative coordinates");
    return rows.map((row) => `<div>${row}</div>`).join("");
  }
  const candidate = entry.item.project;
  const periodic = AXES.filter((_, index) => candidate.periodicAxes?.[index]).map((axis) => axis.toUpperCase());
  return [
    `<strong>${escapeHtml(source)} case template</strong>`,
    `${formatInteger(candidate.blocks.length)} block${candidate.blocks.length === 1 ? "" : "s"}`,
    `${formatInteger(candidate.nscal ?? candidate.nScal ?? 0)} scalars`,
    `Periodic: ${periodic.length ? periodic.join(", ") : "none"}`,
  ].map((row) => `<div>${row}</div>`).join("");
}

function axisKindLabel(kind) {
  return {
    snac: "Native SNaC",
    simple_ratio: "Monotone grading",
    max_min: "Symmetric grading",
    multi: "Multi-region grading",
    explicit: "Explicit coordinates",
  }[kind] ?? String(kind ?? "Axis grid");
}

async function applySelectedLibraryItem() {
  const entry = selectedLibraryEntry();
  if (!entry) return;
  if (currentLibraryTab === "axisPresets") {
    const block = selectedBlock();
    if (!block) return;
    try {
      const axisIndex = AXES.indexOf(currentAxis);
      const applied = applyAxisPreset(entry.item, block.lmin[axisIndex], block.lmax[axisIndex]);
      block.axes[currentAxis] = applied.axis;
      block.ng[axisIndex] = applied.ng;
      markProjectDirty(block.id);
      renderAll();
      schedulePreview(block.id, 0);
      setStatus(`Applied ${entry.item.name} to ${currentAxis.toUpperCase()}`);
    } catch (error) {
      setStatus(error.message);
    }
    return;
  }
  setStatus(`Applying ${entry.item.name}...`);
  try {
    commitHistory();
    const payload = await postJson("/api/migrate", { project: entry.item.project });
    project = payload.project;
    selectedId = project.blocks[0]?.id ?? null;
    selectedIds = new Set(selectedId == null ? [] : [selectedId]);
    lastCheck = null;
    pendingRepair = null;
    decompositionResult = null;
    previewCache.clear();
    previewSequence += 1;
    renderAll();
    commitHistory();
    fitView();
    els.libraryDialog.close();
    const structured = project.inferConnectivity ? await updateStructure({ silent: true }) : true;
    if (structured) setStatus(`Applied ${entry.item.name}`);
  } catch (error) {
    setStatus(error.message);
  }
}

function saveCurrentLibraryItem() {
  const selected = selectedLibraryEntry();
  const existingId = selected && !selected.builtIn ? selected.item.id : undefined;
  try {
    assertLibraryNameAvailable(els.libraryName.value, existingId);
    let item;
    if (currentLibraryTab === "axisPresets") {
      const block = selectedBlock();
      if (!block) throw new Error("Select a block before saving an axis preset");
      const axisIndex = AXES.indexOf(currentAxis);
      item = captureAxisPreset(
        els.libraryName.value,
        ensureAxis(block, currentAxis),
        block.ng[axisIndex],
        block.lmin[axisIndex],
        block.lmax[axisIndex],
        existingId
      );
    } else {
      item = captureProjectTemplate(els.libraryName.value, project, existingId);
    }
    userLibrary = saveLibrary(replaceLibraryItem(userLibrary, currentLibraryTab, item));
    selectedLibraryItem = `user:${item.id}`;
    renderLibraryDialog();
    setStatus(`${existingId ? "Updated" : "Saved"} ${item.name}`);
  } catch (error) {
    setStatus(error.message);
  }
}

function renameSelectedLibraryItem() {
  const entry = selectedLibraryEntry();
  if (!entry || entry.builtIn) return;
  try {
    assertLibraryNameAvailable(els.libraryName.value, entry.item.id);
    userLibrary = saveLibrary(renameLibraryItem(userLibrary, currentLibraryTab, entry.item.id, els.libraryName.value));
    renderLibraryDialog();
    setStatus(`Renamed library item to ${els.libraryName.value}`);
  } catch (error) {
    setStatus(error.message);
  }
}

function deleteSelectedLibraryItem() {
  const entry = selectedLibraryEntry();
  if (!entry || entry.builtIn) return;
  const name = entry.item.name;
  try {
    userLibrary = saveLibrary(removeLibraryItem(userLibrary, currentLibraryTab, entry.item.id));
    selectedLibraryItem = null;
    renderLibraryDialog();
    setStatus(`Deleted ${name}`);
  } catch (error) {
    setStatus(error.message);
  }
}

async function importLibraryJson(event) {
  const file = event.target.files[0];
  if (!file) return;
  event.target.value = "";
  try {
    const incoming = normalizeLibrary(JSON.parse(await file.text()));
    const merged = mergeLibraries(userLibrary, incoming);
    userLibrary = saveLibrary(merged.library);
    selectedLibraryItem = null;
    renderLibraryDialog();
    const count = merged.counts.axisPresets + merged.counts.projectTemplates;
    setStatus(`Imported ${count} library item${count === 1 ? "" : "s"}`);
  } catch (error) {
    setStatus(error.message);
  }
}

function exportLibraryJson() {
  const content = `${JSON.stringify({ ...userLibrary, schemaVersion: LIBRARY_SCHEMA_VERSION }, null, 2)}\n`;
  downloadText(content, "snac-grid-library.json", "application/json");
  setStatus("Downloaded grid library");
}

function assertLibraryNameAvailable(name, exceptId) {
  const target = String(name ?? "").trim().toLocaleLowerCase();
  const duplicate = libraryEntries().some((entry) => entry.item.id !== exceptId && entry.item.name.toLocaleLowerCase() === target);
  if (duplicate) throw new Error("Library item names must be unique");
}

function downloadProjectJson() {
  const blob = new Blob([JSON.stringify(project, null, 2) + "\n"], { type: "application/json" });
  const link = document.createElement("a");
  link.href = URL.createObjectURL(blob);
  link.download = `${project.name || "snac-grid"}.json`;
  link.click();
  URL.revokeObjectURL(link.href);
  commitHistory();
  clearAutosave();
}

async function downloadCaseReport(format) {
  setStatus(`Building ${format === "json" ? "JSON" : "Markdown"} report...`);
  try {
    const payload = await postJson("/api/report", { project });
    const extension = format === "json" ? "json" : "md";
    const type = format === "json" ? "application/json" : "text/markdown";
    const content = format === "json"
      ? `${JSON.stringify(payload.report, null, 2)}\n`
      : payload.markdown;
    const base = (project.name || "snac-grid").replace(/[^A-Za-z0-9._-]+/g, "-");
    downloadText(content, `${base}-quality.${extension}`, type);
    setStatus(`Downloaded ${extension.toUpperCase()} quality report`);
  } catch (error) {
    setStatus(error.message);
  }
}

function downloadText(content, filename, type) {
  const blob = new Blob([content], { type });
  const link = document.createElement("a");
  link.href = URL.createObjectURL(blob);
  link.download = filename;
  link.click();
  URL.revokeObjectURL(link.href);
}

async function openProjectJson(event) {
  const file = event.target.files[0];
  if (!file) return;
  event.target.value = "";
  try {
    const loaded = JSON.parse(await file.text());
    const payload = await postJson("/api/migrate", { project: loaded });
    project = payload.project;
    selectedId = project.blocks[0]?.id ?? null;
    selectedIds = new Set(selectedId == null ? [] : [selectedId]);
    lastCheck = null;
    pendingRepair = null;
    decompositionResult = null;
    previewCache.clear();
    previewSequence += 1;
    renderAll();
    clearAutosave();
    resetHistory(false);
    setStatus(`Loaded ${file.name}`);
  } catch (error) {
    setStatus(error.message);
  }
}

function selectFromScene(event) {
  if (blockTool !== "select") return;
  const rect = renderer.domElement.getBoundingClientRect();
  pointer.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
  pointer.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;
  raycaster.setFromCamera(pointer, camera);
  const hit = raycaster.intersectObjects(blockGroup.children, false).find((item) => item.object.userData.blockId);
  if (hit && hit.object.userData.blockId !== selectedId) {
    selectedId = hit.object.userData.blockId;
    if (event.shiftKey) selectedIds.add(selectedId);
    renderAll();
  }
}

function animate() {
  requestAnimationFrame(animate);
  controls.update();
  renderer.render(scene, camera);
}

function resizeRenderer() {
  const rect = els.scene.parentElement.getBoundingClientRect();
  const width = Math.max(1, Math.floor(rect.width));
  const height = Math.max(1, Math.floor(rect.height));
  renderer.setSize(width, height, false);
  camera.aspect = width / height;
  camera.updateProjectionMatrix();
}

function fitView() {
  const box = projectBox(project.blocks);
  if (box.isEmpty()) {
    controls.target.set(0, 0, 0);
    camera.position.set(3.0, -4.5, 3.0);
    controls.update();
    return;
  }
  const center = box.getCenter(new THREE.Vector3());
  const size = box.getSize(new THREE.Vector3());
  const radius = Math.max(size.x, size.y, size.z, 1);
  camera.position.copy(center.clone().add(new THREE.Vector3(radius * 1.6, -radius * 2.2, radius * 1.45)));
  controls.target.copy(center);
  controls.update();
  gridHelper.position.copy(center);
  gridHelper.scale.setScalar(Math.max(1, radius / 5));
}

function selectedBlock() {
  return project?.blocks?.find((block) => block.id === selectedId);
}

function selectedBlockIds() {
  const ids = [...selectedIds].filter((id) => project.blocks.some((block) => block.id === id));
  return ids.length ? ids : selectedId == null ? [] : [selectedId];
}

function blockLabel(block) {
  return block.name ? `${block.id}: ${block.name}` : `Block ${block.id}`;
}

function numberInput(field, index, value, min = null, disabled = false) {
  const minAttr = min == null ? "" : `min="${min}"`;
  const disabledAttr = disabled ? 'disabled title="Derived from min/max spacing"' : "";
  return `<label class="field"><input data-field="${field}" data-index="${index}" type="number" step="any" ${minAttr} ${disabledAttr} value="${value}"></label>`;
}

function boundaryComponent(name, face, code, value, component = null, forcedKind = null) {
  const kind = forcedKind ?? (component == null ? "pre" : "vel");
  const componentAttr = component == null ? "" : `data-component="${component}"`;
  return `
    <div class="boundary-component">
      <span class="component-name">${name}</span>
      ${bcSelect(kind, face, code, componentAttr)}
      ${bcValueInput(kind, face, value, componentAttr)}
    </div>
  `;
}

function bcSelect(kind, face, selected, componentAttr = "") {
  return `<select data-kind="${kind}" data-face="${face}" ${componentAttr}>
    ${["D", "N", "F"].map((value) => `<option value="${value}" ${value === selected ? "selected" : ""}>${value}</option>`).join("")}
  </select>`;
}

function bcValueInput(kind, face, value, componentAttr = "") {
  return `
    <input data-kind="${kind}" data-face="${face}" ${componentAttr} type="number" step="any" value="${value ?? 0}">
  `;
}

function positiveNumber(value, fallback) {
  const number = Number(value);
  return Number.isFinite(number) && number > 0 ? number : fallback;
}

function formatNumber(value) {
  if (!Number.isFinite(value)) return "n/a";
  if (Math.abs(value) >= 1000 || Math.abs(value) < 0.001) return value.toExponential(3);
  return value.toPrecision(5);
}

function formatInteger(value) {
  return Number.isFinite(value) ? Math.round(value).toLocaleString() : "n/a";
}

function formatPercent(value) {
  if (!Number.isFinite(value)) return "n/a";
  return `${(100 * value).toPrecision(value < 0.001 ? 2 : 3)}%`;
}

function formatBytes(value) {
  let size = Number(value) || 0;
  const units = ["B", "KiB", "MiB", "GiB"];
  let index = 0;
  while (size >= 1024 && index < units.length - 1) {
    size /= 1024;
    index += 1;
  }
  return `${formatNumber(size)} ${units[index]}`;
}

function formatInputNumber(value) {
  return Number(value.toPrecision(12)).toString();
}

function setStatus(message) {
  els.status.textContent = message;
}

function clearCheckState() {
  lastCheck = null;
  renderCheckResults();
  renderQualitySummary();
  setStatus("");
}

function decompositionSettingEdited() {
  pendingRepair = null;
  decompositionResult = null;
  lastCheck = null;
}

function markProjectDirty(blockId = null) {
  if (blockId != null) previewCache.delete(blockId);
  pendingRepair = null;
  decompositionResult = null;
  clearCheckState();
  renderDecompositionSummary();
  scheduleHistoryCommit();
}

function projectSnapshot() {
  return {
    project: JSON.stringify(project),
    decompositionResult: decompositionResult ? JSON.stringify(decompositionResult) : null,
    selectedId,
    selectedIds: [...selectedIds],
    currentAxis,
  };
}

function resetHistory(save = true) {
  window.clearTimeout(historyTimer);
  historyTimer = null;
  projectHistory.reset(projectSnapshot());
  updateHistoryButtons();
  if (save) saveAutosave();
}

function scheduleHistoryCommit(delay = 250) {
  if (isTransformDragging) return;
  saveAutosave();
  window.clearTimeout(historyTimer);
  historyTimer = window.setTimeout(commitHistory, delay);
}

function commitHistory() {
  window.clearTimeout(historyTimer);
  historyTimer = null;
  const snapshot = projectSnapshot();
  saveAutosave();
  projectHistory.commit(snapshot);
  updateHistoryButtons();
}

function undoProject() {
  commitHistory();
  const state = projectHistory.undo();
  if (state) restoreHistory(state);
}

function redoProject() {
  const state = projectHistory.redo();
  if (state) restoreHistory(state);
}

function restoreHistory(state) {
  window.clearTimeout(historyTimer);
  historyTimer = null;
  project = JSON.parse(state.project);
  selectedId = project.blocks.some((block) => block.id === state.selectedId) ? state.selectedId : project.blocks[0]?.id ?? null;
  selectedIds = new Set((state.selectedIds ?? []).filter((id) => project.blocks.some((block) => block.id === id)));
  currentAxis = state.currentAxis;
  lastCheck = null;
  pendingRepair = null;
  decompositionResult = state.decompositionResult ? JSON.parse(state.decompositionResult) : null;
  previewCache.clear();
  previewSequence += 1;
  renderAll();
  setStatus("");
}

function saveAutosave() {
  saveStoredProject(AUTOSAVE_KEY, project);
}

function clearAutosave() {
  clearStoredProject(AUTOSAVE_KEY);
}

function loadAutosave() {
  return loadStoredProject(AUTOSAVE_KEY);
}

function showRecoveryDialog(recovery) {
  pendingRecovery = recovery;
  const savedAt = new Date(recovery.savedAt);
  const time = Number.isNaN(savedAt.getTime()) ? "an earlier session" : savedAt.toLocaleString();
  els.recoverySummary.textContent = `${recovery.project.name || "Untitled project"}, saved ${time}`;
  els.recoveryDialog.showModal();
}

async function restoreRecovery() {
  if (!pendingRecovery) return;
  try {
    const payload = await postJson("/api/migrate", { project: pendingRecovery.project });
    project = payload.project;
    selectedId = project.blocks[0]?.id ?? null;
    selectedIds = new Set(selectedId == null ? [] : [selectedId]);
    pendingRecovery = null;
    els.recoveryDialog.close();
    previewCache.clear();
    previewSequence += 1;
    renderAll();
    resetHistory();
    fitView();
    setStatus("Recovered autosaved project");
  } catch (error) {
    setStatus(error.message);
  }
}

function discardRecovery() {
  pendingRecovery = null;
  clearAutosave();
  els.recoveryDialog.close();
  resetHistory(false);
  setStatus("Started with the sample project");
}

function updateHistoryButtons() {
  if (!els.undoProject || !els.redoProject) return;
  els.undoProject.disabled = !projectHistory.canUndo;
  els.redoProject.disabled = !projectHistory.canRedo;
}

function escapeHtml(value) {
  return String(value).replace(/[&<>"']/g, (char) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#039;" })[char]);
}

function toCamel(id) {
  return id.replace(/-([a-z])/g, (_, char) => char.toUpperCase());
}
