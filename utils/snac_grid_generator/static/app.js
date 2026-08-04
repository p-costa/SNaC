import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { TransformControls } from "three/addons/controls/TransformControls.js";

const AXES = ["x", "y", "z"];
const FACE_ORDER = ["x-", "x+", "y-", "y+", "z-", "z+"];
const INIT_FIELDS = ["zer", "uni", "cou", "poi", "log", "hcp", "hcl", "tgv"];
const GRID_FUNCTIONS = [
  "Cluster both ends",
  "Cluster lower end",
  "Cluster middle",
  "Cluster upper end",
  "Geometric lower end",
  "Geometric upper end",
  "Geometric both ends",
  "Geometric middle",
];
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
let currentAxis = "x";
let groundGridVisible = true;
let boundaryColorsVisible = true;
let partitionPlanesVisible = true;
let blockTool = "select";
let lastCheck = null;
let previewTimer = null;
let previewSequence = 0;
const previewCache = new Map();
let historyTimer = null;
let historyIndex = -1;
let projectHistory = [];
const HISTORY_LIMIT = 100;

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
let transformControls;
let selectedBlockMesh = null;
let transformDrag = null;
let isTransformDragging = false;
const axisLabelMaterials = new Map();

bootstrap();

async function bootstrap() {
  bindElements();
  initScene();
  bindStaticEvents();
  const response = await fetch("/api/sample");
  project = await response.json();
  selectedId = project.blocks[0]?.id ?? null;
  renderAll();
  resetHistory();
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
    "download-json",
    "open-json",
    "undo-project",
    "redo-project",
    "reset-project",
    "status",
    "block-list",
    "add-block",
    "add-adjacent",
    "adjacent-face",
    "adjacent-size",
    "scene",
    "fit-view",
    "toggle-grid",
    "tool-select",
    "tool-move",
    "tool-scale",
    "toggle-boundaries",
    "toggle-partitions",
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
    "check-results",
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
    for (const block of project.blocks) ensureBoundaryArrays(block);
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
  els.exportCase.addEventListener("click", exportCase);
  els.checkProject.addEventListener("click", () => {
    if (project.inferConnectivity) updateStructure({ successMessage: "Grid checks passed" });
    else checkProject();
  });
  els.updateStructure.addEventListener("click", updateStructure);
  els.downloadJson.addEventListener("click", downloadProjectJson);
  els.openJson.addEventListener("change", openProjectJson);
  els.undoProject.addEventListener("click", undoProject);
  els.redoProject.addEventListener("click", redoProject);
  els.resetProject.addEventListener("click", resetProject);
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
  partitionGroup = new THREE.Group();
  axesGroup = new THREE.Group();
  gridHelper = new THREE.GridHelper(10, 20, 0x46515c, 0x252b32);
  gridHelper.rotation.x = Math.PI / 2;
  scene.add(gridHelper);
  scene.add(axesGroup);
  scene.add(blockGroup);
  scene.add(boundaryFaceGroup);
  scene.add(selectedGridGroup);
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
  project.schemaVersion = project.schemaVersion ?? 1;
  project.blocks = project.blocks || [];
  project.nscal = normalizedScalarCount(project.nscal ?? project.nScal);
  project.periodicAxes = normalizedPeriodicAxes();
  for (const block of project.blocks) ensureBoundaryArrays(block);
  project.blocks.sort((a, b) => a.id - b.id);
  if (!selectedBlock()) selectedId = project.blocks[0]?.id ?? null;
  els.projectName.value = project.name ?? "snac-grid";
  els.scalarCount.value = project.nscal;
  els.inferConnectivity.checked = project.inferConnectivity !== false;
  project.writeExternalGrid = project.writeExternalGrid !== false;
  els.writeExternalGrid.checked = project.writeExternalGrid;
  if (!["grid", "both"].includes(project.externalGridSource)) project.externalGridSource = "grid";
  els.gridSource.value = project.externalGridSource;
  els.gridSource.disabled = !project.writeExternalGrid;
  els.addAdjacent.disabled = !selectedBlock();
  renderBlockList();
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
    const button = document.createElement("button");
    button.className = `block-item${block.id === selectedId ? " selected" : ""}`;
    button.innerHTML = `
      <span class="swatch" style="background:#${COLORS[(block.id - 1) % COLORS.length].toString(16).padStart(6, "0")}"></span>
      <span><strong>${escapeHtml(blockLabel(block))}</strong><span>${block.lmin.join(", ")} -> ${block.lmax.join(", ")}</span></span>
    `;
    button.addEventListener("click", () => {
      selectedId = block.id;
      renderAll();
    });
    els.blockList.appendChild(button);
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
      block[field][index] = value;
      if (field === "dims") propagateMpiPartition(block.id, index, value);
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
  const selectedKind = axis.kind ?? "snac";
  els.axisEditor.innerHTML = `
    <label class="field">
      <span>Mode</span>
      <select id="grid-kind">
        <option value="snac">Native SNaC</option>
        <option value="simple_ratio">Monotone grading</option>
        <option value="max_min">Symmetric grading</option>
        <option value="multi">Multi-region grading</option>
      </select>
    </label>
    <div id="grid-kind-body"></div>
    <button id="apply-axis" class="button" type="button"><i data-lucide="copy-check"></i><span>Apply to aligned blocks</span></button>
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
  if (axis.kind === "multi") {
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
    renderAll();
    scheduleHistoryCommit();
    const count = payload.changedBlockIds?.length ?? 0;
    setStatus(count ? `Applied ${currentAxis.toUpperCase()} grid to ${count} aligned block${count === 1 ? "" : "s"}` : "No other aligned blocks");
  } catch (error) {
    setStatus(error.message);
  }
}

function renderBoundaryTable(block) {
  ensureBoundaryArrays(block);
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
  els.exportCase.disabled = !project.blocks.length || errors.length > 0;
  if (!lastCheck) {
    els.checkResults.innerHTML = "";
    return;
  }
  if (!errors.length && !warnings.length) {
    els.checkResults.innerHTML = `<div class="check-line ok">Structured grid checks passed</div>`;
    return;
  }
  els.checkResults.innerHTML = `
    ${errors.map((message) => `<div class="check-line error">${escapeHtml(message)}</div>`).join("")}
    ${warnings.map((message) => `<div class="check-line warning">${escapeHtml(message)}</div>`).join("")}
  `;
}

function renderScene() {
  transformControls.detach();
  selectedBlockMesh = null;
  boundaryFaceGroup.visible = true;
  selectedGridGroup.visible = true;
  partitionGroup.visible = true;
  clearDisposableGroup(blockGroup);
  clearDisposableGroup(boundaryFaceGroup);
  clearDisposableGroup(selectedGridGroup);
  clearDisposableGroup(partitionGroup);
  for (const block of project.blocks) {
    const size = blockSize(block);
    const center = blockCenter(block);
    const color = COLORS[(block.id - 1) % COLORS.length];
    const material = new THREE.MeshStandardMaterial({
      color,
      transparent: true,
      opacity: block.id === selectedId ? 0.43 : 0.24,
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
      new THREE.LineBasicMaterial({ color: block.id === selectedId ? 0xffffff : 0x9aa6b2, transparent: true, opacity: block.id === selectedId ? 0.95 : 0.48 })
    );
    mesh.add(edges);
  }
  drawBoundaryFaces();
  drawSelectedGrid();
  drawPartitionPlanes();
  updateAxesHelper();
  attachSelectedTransform();
  renderViewControls();
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

function partitionFaceIndices(n, partitions) {
  if (partitions <= 1) return [];
  const base = Math.floor(n / partitions);
  const remainder = n % partitions;
  return Array.from({ length: partitions - 1 }, (_, rank) => {
    const count = rank + 1;
    return count * base + Math.min(count, remainder);
  });
}

function planeOutlinePoints(block, axisIndex, value) {
  const [xmin, ymin, zmin] = block.lmin;
  const [xmax, ymax, zmax] = block.lmax;
  if (axisIndex === 0) {
    return [
      new THREE.Vector3(value, ymin, zmin),
      new THREE.Vector3(value, ymax, zmin),
      new THREE.Vector3(value, ymax, zmax),
      new THREE.Vector3(value, ymin, zmax),
      new THREE.Vector3(value, ymin, zmin),
    ];
  }
  if (axisIndex === 1) {
    return [
      new THREE.Vector3(xmin, value, zmin),
      new THREE.Vector3(xmax, value, zmin),
      new THREE.Vector3(xmax, value, zmax),
      new THREE.Vector3(xmin, value, zmax),
      new THREE.Vector3(xmin, value, zmin),
    ];
  }
  return [
    new THREE.Vector3(xmin, ymin, value),
    new THREE.Vector3(xmax, ymin, value),
    new THREE.Vector3(xmax, ymax, value),
    new THREE.Vector3(xmin, ymax, value),
    new THREE.Vector3(xmin, ymin, value),
  ];
}

function updateAxesHelper() {
  clearDisposableGroup(axesGroup, false);
  const box = projectBox();
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

function axisLabel(label, color) {
  const key = `${label}-${color}`;
  if (axisLabelMaterials.has(key)) {
    const sprite = new THREE.Sprite(axisLabelMaterials.get(key));
    sprite.userData.sharedMaterial = true;
    return sprite;
  }
  const canvas = document.createElement("canvas");
  canvas.width = 64;
  canvas.height = 64;
  const ctx = canvas.getContext("2d");
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  ctx.fillStyle = `#${color.toString(16).padStart(6, "0")}`;
  ctx.font = "700 42px Inter, Arial, sans-serif";
  ctx.textAlign = "center";
  ctx.textBaseline = "middle";
  ctx.fillText(label, 32, 34);
  const texture = new THREE.CanvasTexture(canvas);
  texture.colorSpace = THREE.SRGBColorSpace;
  const material = new THREE.SpriteMaterial({ map: texture, transparent: true, depthTest: false });
  axisLabelMaterials.set(key, material);
  const sprite = new THREE.Sprite(material);
  sprite.userData.sharedMaterial = true;
  return sprite;
}

function drawSelectedGrid() {
  const block = selectedBlock();
  if (!block) return;
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
    block.ng = payload.ng;
    scheduleHistoryCommit();
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

function clearDisposableGroup(group, disposeGeometries = true) {
  const geometries = new Set();
  const materials = new Set();
  for (const child of [...group.children]) {
    child.traverse((object) => {
      if (object.userData.sharedMaterial) return;
      if (disposeGeometries && object.geometry) geometries.add(object.geometry);
      if (!object.material) return;
      for (const material of Array.isArray(object.material) ? object.material : [object.material]) materials.add(material);
    });
    group.remove(child);
  }
  for (const geometry of geometries) geometry.dispose();
  for (const material of materials) material.dispose();
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
}

function addFreeBlock() {
  const source = selectedBlock() ?? project.blocks[0];
  const copy = source ? structuredClone(source) : defaultBlock(nextBlockId());
  copy.id = nextBlockId();
  copy.name = `block-${copy.id}`;
  if (source) {
    const width = source.lmax[0] - source.lmin[0];
    copy.lmin[0] += width;
    copy.lmax[0] += width;
  }
  resetFriendBoundaries(copy);
  project.blocks.push(copy);
  selectedId = copy.id;
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
  copy.id = nextBlockId();
  copy.name = `block-${copy.id}`;
  if (sign > 0) {
    copy.lmin[axis] = source.lmax[axis];
    copy.lmax[axis] = source.lmax[axis] + thickness;
  } else {
    copy.lmax[axis] = source.lmin[axis];
    copy.lmin[axis] = source.lmin[axis] - thickness;
  }
  resetFriendBoundaries(copy);
  project.blocks.push(copy);
  selectedId = copy.id;
  markProjectDirty(copy.id);
  renderAll();
}

function removeSelectedBlock() {
  if (!selectedBlock()) return;
  project.blocks = project.blocks.filter((block) => block.id !== selectedId);
  previewCache.delete(selectedId);
  selectedId = project.blocks[0]?.id ?? null;
  markProjectDirty();
  renderAll();
}

function resetProject() {
  project.blocks = [];
  selectedId = null;
  lastCheck = null;
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
    const warningText = payload.warnings?.length ? `, ${payload.warnings.length} warnings` : "";
    setStatus(`Wrote ${payload.files.length} files to ${payload.outputDir}${warningText}`);
  } catch (error) {
    setStatus(error.message);
  }
}

async function postJson(url, body, options = {}) {
  const response = await fetch(url, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
  const payload = await response.json();
  if (!response.ok || (payload.ok === false && !options.allowInvalid)) {
    throw new Error(payload.error || payload.errors?.join("; ") || "request failed");
  }
  return payload;
}

function downloadProjectJson() {
  const blob = new Blob([JSON.stringify(project, null, 2) + "\n"], { type: "application/json" });
  const link = document.createElement("a");
  link.href = URL.createObjectURL(blob);
  link.download = `${project.name || "snac-grid"}.json`;
  link.click();
  URL.revokeObjectURL(link.href);
}

function openProjectJson(event) {
  const file = event.target.files[0];
  if (!file) return;
  const reader = new FileReader();
  reader.onload = () => {
    try {
      project = JSON.parse(reader.result);
      project.blocks = project.blocks || [];
      selectedId = project.blocks[0]?.id ?? null;
      lastCheck = null;
      previewCache.clear();
      previewSequence += 1;
      renderAll();
      resetHistory();
      setStatus(`Loaded ${file.name}`);
    } catch (error) {
      setStatus(error.message);
    }
  };
  reader.readAsText(file);
  event.target.value = "";
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

function projectBox() {
  const box = new THREE.Box3();
  for (const block of project.blocks) {
    box.expandByPoint(new THREE.Vector3(...block.lmin));
    box.expandByPoint(new THREE.Vector3(...block.lmax));
  }
  if (box.isEmpty()) {
    box.expandByPoint(new THREE.Vector3(0, 0, 0));
    box.expandByPoint(new THREE.Vector3(1, 1, 1));
  }
  return box;
}

function fitView() {
  const box = projectBox();
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

function blockLabel(block) {
  return block.name ? `${block.id}: ${block.name}` : `Block ${block.id}`;
}

function blockSize(block) {
  return block.lmax.map((value, index) => Math.max(1e-9, value - block.lmin[index]));
}

function blockCenter(block) {
  return block.lmin.map((value, index) => 0.5 * (value + block.lmax[index]));
}

function ensureAxis(block, axis) {
  block.axes = block.axes || {};
  block.axes[axis] = block.axes[axis] || defaultAxis();
  block.axes[axis].segments = block.axes[axis].segments || [{ length: 1, cells: 1, ratio: 1, continuous: false }];
  for (const segment of block.axes[axis].segments) segment.continuous = Boolean(segment.continuous);
  block.axes[axis].profile = block.axes[axis].profile || "geometric";
  block.axes[axis].controls = block.axes[axis].controls || ["n", "ratio"];
  return block.axes[axis];
}

function isCellCountDerived(axis) {
  return axis.kind === "max_min" || (axis.kind === "simple_ratio" && !axis.controls?.includes("n"));
}

function defaultAxis() {
  return {
    kind: "snac",
    gt: 0,
    gr: 0,
    ratio: 1,
    cell_ratio: 1,
    width_start: 0.01,
    width_end: 0.02,
    controls: ["n", "ratio"],
    profile: "geometric",
    min: 0.01,
    max: 0.02,
    side: "end",
    segments: [{ length: 1, cells: 1, ratio: 1, continuous: false }],
  };
}

function defaultBlock(id) {
  return {
    id,
    name: "",
    dims: [1, 1, 1],
    ng: [32, 32, 2],
    lmin: [0, 0, 0],
    lmax: [1, 1, 0.05],
    axes: { x: defaultAxis(), y: defaultAxis(), z: defaultAxis() },
    cbcvel: Array.from({ length: 3 }, () => ["D", "D", "D", "D", "D", "D"]),
    cbcpre: ["N", "N", "N", "N", "N", "N"],
    bcvel: Array.from({ length: 3 }, () => [0, 0, 0, 0, 0, 0]),
    bcpre: [0, 0, 0, 0, 0, 0],
    cbcscal: [],
    bcscal: [],
    inflow: [0, 0, 0, 0, 0, 0],
    inivel: "zer",
  };
}

function ensureBoundaryArrays(block) {
  const nscal = normalizedScalarCount(project.nscal);
  block.cbcvel = Array.from({ length: 3 }, (_, component) => completeArray(block.cbcvel?.[component], "D"));
  block.cbcpre = completeArray(block.cbcpre, "N");
  block.bcvel = Array.from({ length: 3 }, (_, component) => completeArray(block.bcvel?.[component], 0));
  block.bcpre = completeArray(block.bcpre, 0);
  block.cbcscal = Array.from({ length: nscal }, (_, iscal) => completeArray(block.cbcscal?.[iscal], "N"));
  block.bcscal = Array.from({ length: nscal }, (_, iscal) => completeArray(block.bcscal?.[iscal], 0));
  block.inflow = completeArray(block.inflow, 0);
}

function completeArray(values, fallback) {
  const result = Array.isArray(values) ? values.slice(0, 6) : [];
  while (result.length < 6) result.push(fallback);
  return result;
}

function normalizedPeriodicAxes() {
  const values = Array.isArray(project.periodicAxes) ? project.periodicAxes : [false, false, false];
  return AXES.map((_, index) => Boolean(values[index]));
}

function normalizedScalarCount(value = project?.nscal ?? 0) {
  return Math.max(0, Math.round(Number(value) || 0));
}

function resetFriendBoundaries(block) {
  block.cbcvel = Array.from({ length: 3 }, () => ["D", "D", "D", "D", "D", "D"]);
  block.cbcpre = ["N", "N", "N", "N", "N", "N"];
  block.bcvel = Array.from({ length: 3 }, () => [0, 0, 0, 0, 0, 0]);
  block.bcpre = [0, 0, 0, 0, 0, 0];
  block.cbcscal = Array.from({ length: normalizedScalarCount(project.nscal) }, () => ["N", "N", "N", "N", "N", "N"]);
  block.bcscal = Array.from({ length: normalizedScalarCount(project.nscal) }, () => [0, 0, 0, 0, 0, 0]);
}

function propagateMpiPartition(sourceId, axisIndex, value) {
  const graph = new Map(project.blocks.map((block) => [block.id, new Set()]));
  for (const connection of fullFaceConnections(project.blocks)) {
    if (connection.axisIndex === axisIndex) continue;
    graph.get(connection.a.id).add(connection.b.id);
    graph.get(connection.b.id).add(connection.a.id);
  }
  const queue = [sourceId];
  const seen = new Set(queue);
  while (queue.length) {
    const blockId = queue.shift();
    const block = project.blocks.find((item) => item.id === blockId);
    if (block) block.dims[axisIndex] = value;
    for (const neighborId of graph.get(blockId) ?? []) {
      if (!seen.has(neighborId)) {
        seen.add(neighborId);
        queue.push(neighborId);
      }
    }
  }
}

function fullFaceConnections(blocks) {
  const connections = [];
  for (let i = 0; i < blocks.length; i += 1) {
    for (let j = i + 1; j < blocks.length; j += 1) {
      const a = blocks[i];
      const b = blocks[j];
      for (let axisIndex = 0; axisIndex < 3; axisIndex += 1) {
        if (touches(a.lmax[axisIndex], b.lmin[axisIndex]) && sameCrossSection(a, b, axisIndex)) {
          connections.push({ a, b, axisIndex });
        }
        if (touches(b.lmax[axisIndex], a.lmin[axisIndex]) && sameCrossSection(a, b, axisIndex)) {
          connections.push({ a: b, b: a, axisIndex });
        }
      }
    }
  }
  return connections;
}

function sameCrossSection(a, b, axisIndex) {
  for (let index = 0; index < 3; index += 1) {
    if (index === axisIndex) continue;
    if (!touches(a.lmin[index], b.lmin[index]) || !touches(a.lmax[index], b.lmax[index])) return false;
  }
  return true;
}

function touches(a, b) {
  return Math.abs(a - b) <= 1e-10;
}

function nextBlockId() {
  return Math.max(0, ...project.blocks.map((block) => block.id)) + 1;
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

function formatPercent(value) {
  if (!Number.isFinite(value)) return "n/a";
  return `${(100 * value).toPrecision(value < 0.001 ? 2 : 3)}%`;
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
  setStatus("");
}

function markProjectDirty(blockId = null) {
  if (blockId != null) previewCache.delete(blockId);
  clearCheckState();
  scheduleHistoryCommit();
}

function projectSnapshot() {
  return {
    project: JSON.stringify(project),
    selectedId,
    currentAxis,
  };
}

function resetHistory() {
  window.clearTimeout(historyTimer);
  historyTimer = null;
  projectHistory = [projectSnapshot()];
  historyIndex = 0;
  updateHistoryButtons();
}

function scheduleHistoryCommit(delay = 250) {
  if (isTransformDragging) return;
  window.clearTimeout(historyTimer);
  historyTimer = window.setTimeout(commitHistory, delay);
}

function commitHistory() {
  window.clearTimeout(historyTimer);
  historyTimer = null;
  const snapshot = projectSnapshot();
  if (projectHistory[historyIndex]?.project === snapshot.project) return;
  projectHistory = projectHistory.slice(0, historyIndex + 1);
  projectHistory.push(snapshot);
  if (projectHistory.length > HISTORY_LIMIT) projectHistory.shift();
  historyIndex = projectHistory.length - 1;
  updateHistoryButtons();
}

function undoProject() {
  commitHistory();
  if (historyIndex <= 0) return;
  restoreHistory(historyIndex - 1);
}

function redoProject() {
  if (historyIndex >= projectHistory.length - 1) return;
  restoreHistory(historyIndex + 1);
}

function restoreHistory(index) {
  window.clearTimeout(historyTimer);
  historyTimer = null;
  historyIndex = index;
  const state = projectHistory[historyIndex];
  project = JSON.parse(state.project);
  selectedId = project.blocks.some((block) => block.id === state.selectedId) ? state.selectedId : project.blocks[0]?.id ?? null;
  currentAxis = state.currentAxis;
  lastCheck = null;
  previewCache.clear();
  previewSequence += 1;
  renderAll();
  setStatus("");
}

function updateHistoryButtons() {
  if (!els.undoProject || !els.redoProject) return;
  els.undoProject.disabled = historyIndex <= 0;
  els.redoProject.disabled = historyIndex < 0 || historyIndex >= projectHistory.length - 1;
}

function escapeHtml(value) {
  return String(value).replace(/[&<>"']/g, (char) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#039;" })[char]);
}

function toCamel(id) {
  return id.replace(/-([a-z])/g, (_, char) => char.toUpperCase());
}
