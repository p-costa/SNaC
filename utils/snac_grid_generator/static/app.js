import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";

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
const COLORS = [0x43c989, 0x76a9ff, 0xd28a3e, 0xc789e8, 0xe46868, 0x82c6c2, 0xd6bf5b, 0x9ba7ff];

let project = null;
let selectedId = 1;
let currentAxis = "x";
let groundGridVisible = true;
let lastCheck = null;

const els = {};
let scene;
let camera;
let renderer;
let controls;
let raycaster;
let pointer;
let blockGroup;
let gridHelper;
let axesGroup;
let selectedGridGroup;
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
  });
  els.scalarCount.addEventListener("input", () => {
    project.nscal = normalizedScalarCount(els.scalarCount.value);
    for (const block of project.blocks) ensureBoundaryArrays(block);
    renderInspector();
    clearCheckState();
  });
  els.inferConnectivity.addEventListener("change", () => {
    project.inferConnectivity = els.inferConnectivity.checked;
  });
  els.gridSource.addEventListener("change", () => {
    project.externalGridSource = els.gridSource.value;
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
  els.exportCase.addEventListener("click", exportCase);
  els.checkProject.addEventListener("click", () => checkProject());
  els.updateStructure.addEventListener("click", updateStructure);
  els.downloadJson.addEventListener("click", downloadProjectJson);
  els.openJson.addEventListener("change", openProjectJson);
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
  selectedGridGroup = new THREE.Group();
  axesGroup = new THREE.Group();
  gridHelper = new THREE.GridHelper(10, 20, 0x46515c, 0x252b32);
  gridHelper.rotation.x = Math.PI / 2;
  scene.add(gridHelper);
  scene.add(axesGroup);
  scene.add(blockGroup);
  scene.add(selectedGridGroup);

  const ambient = new THREE.HemisphereLight(0xffffff, 0x1a2028, 1.9);
  const key = new THREE.DirectionalLight(0xffffff, 1.4);
  key.position.set(4, -5, 6);
  scene.add(ambient, key);
  resizeRenderer();
}

function renderAll() {
  if (!project) return;
  project.blocks = project.blocks || [];
  project.nscal = normalizedScalarCount(project.nscal ?? project.nScal);
  project.periodicAxes = normalizedPeriodicAxes();
  for (const block of project.blocks) ensureBoundaryArrays(block);
  project.blocks.sort((a, b) => a.id - b.id);
  if (!selectedBlock()) selectedId = project.blocks[0]?.id ?? null;
  els.projectName.value = project.name ?? "snac-grid";
  els.scalarCount.value = project.nscal;
  els.inferConnectivity.checked = project.inferConnectivity !== false;
  if (!["grid", "both"].includes(project.externalGridSource)) project.externalGridSource = "grid";
  els.gridSource.value = project.externalGridSource;
  els.addAdjacent.disabled = !selectedBlock();
  renderBlockList();
  renderInspector();
  renderScene();
  renderSpacingPreview();
  renderCheckResults();
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
    renderBlockList();
    renderScene();
  };

  els.inivel.innerHTML = INIT_FIELDS.map((value) => `<option value="${value}">${value}</option>`).join("");
  els.inivel.value = block.inivel ?? "zer";
  els.inivel.onchange = () => {
    block.inivel = els.inivel.value;
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
    <span class="axis-label">ng</span>${AXES.map((axis, index) => numberInput("ng", index, block.ng[index], 1)).join("")}
    <span class="axis-label">mpi</span>${AXES.map((axis, index) => numberInput("dims", index, block.dims[index], 1)).join("")}
  `;
  for (const input of [...els.geometryGrid.querySelectorAll("input"), ...els.resolutionGrid.querySelectorAll("input")]) {
    input.addEventListener("input", () => {
      const field = input.dataset.field;
      const index = Number(input.dataset.index);
      const value = field === "ng" || field === "dims" ? Math.max(1, Math.round(Number(input.value) || 1)) : Number(input.value);
      block[field][index] = value;
      if (field === "dims") propagateMpiPartition(block.id, index, value);
      if (field === "ng" || field === "lmin" || field === "lmax") renderSpacingPreview();
      renderScene();
      renderBlockList();
    });
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
        <option value="snac">SNaC mapping</option>
        <option value="multi">OpenFOAM multi-grading</option>
        <option value="simple_ratio">Single ratio</option>
        <option value="max_min">Min/max spacing</option>
      </select>
    </label>
    <div id="grid-kind-body"></div>
  `;
  const kindSelect = document.getElementById("grid-kind");
  kindSelect.value = selectedKind;
  kindSelect.onchange = () => {
    axis.kind = kindSelect.value;
    renderAxisEditor(block);
    renderSpacingPreview();
    renderScene();
  };
  renderAxisKindBody(block, axis);
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

function bindProfileSelect(axis) {
  const select = document.getElementById("grid-profile");
  if (!select) return;
  select.onchange = () => {
    axis.profile = select.value;
    renderSpacingPreview();
    renderScene();
  };
}

function renderAxisKindBody(block, axis) {
  const body = document.getElementById("grid-kind-body");
  if (axis.kind === "multi") {
    axis.segments = axis.segments?.length ? axis.segments : [{ length: 1, cells: 1, ratio: 1 }];
    body.innerHTML = `
      ${profileSelect(axis)}
      <div class="segment-table">
        ${axis.segments
          .map(
            (segment, index) => `
          <div class="segment-row" data-segment="${index}">
            <label class="field"><span>Length</span><input data-key="length" type="number" step="any" min="1e-12" value="${segment.length}"></label>
            <label class="field"><span>Cells</span><input data-key="cells" type="number" step="any" min="1e-12" value="${segment.cells}"></label>
            <label class="field"><span>Ratio</span><input data-key="ratio" type="number" step="any" min="1e-12" value="${segment.ratio}"></label>
            <button class="button icon segment-remove" title="Remove segment"><i data-lucide="minus"></i></button>
          </div>
        `
          )
          .join("")}
      </div>
      <button id="add-segment" class="button"><i data-lucide="list-plus"></i><span>Add Segment</span></button>
    `;
    bindProfileSelect(axis);
    for (const row of body.querySelectorAll(".segment-row")) {
      const index = Number(row.dataset.segment);
      for (const input of row.querySelectorAll("input")) {
        input.oninput = () => {
          axis.segments[index][input.dataset.key] = positiveNumber(input.value, 1);
          renderSpacingPreview();
          renderScene();
        };
      }
      row.querySelector(".segment-remove").onclick = () => {
        if (axis.segments.length > 1) axis.segments.splice(index, 1);
        renderAxisKindBody(block, axis);
        renderSpacingPreview();
        renderScene();
      };
    }
    document.getElementById("add-segment").onclick = () => {
      axis.segments.push({ length: 1, cells: 1, ratio: 1 });
      renderAxisKindBody(block, axis);
      renderSpacingPreview();
      renderScene();
    };
  } else if (axis.kind === "simple_ratio") {
    body.innerHTML = `
      ${profileSelect(axis)}
      <div class="form-grid two">
        <label class="field"><span>Expansion</span><input id="simple-ratio" type="number" step="any" min="1e-12" value="${axis.ratio ?? 1}"></label>
        <label class="field"><span>Bias</span>
          <select id="simple-side">
            <option value="end">Lower end</option>
            <option value="start">Upper end</option>
          </select>
        </label>
      </div>
    `;
    bindProfileSelect(axis);
    document.getElementById("simple-side").value = axis.side ?? "end";
    document.getElementById("simple-ratio").oninput = (event) => {
      axis.ratio = positiveNumber(event.target.value, 1);
      renderSpacingPreview();
      renderScene();
    };
    document.getElementById("simple-side").onchange = (event) => {
      axis.side = event.target.value;
      renderSpacingPreview();
      renderScene();
    };
  } else if (axis.kind === "max_min") {
    body.innerHTML = `
      ${profileSelect(axis)}
      <div class="form-grid two">
        <label class="field"><span>Min dx</span><input id="min-dx" type="number" step="any" min="1e-12" value="${axis.min ?? 0.01}"></label>
        <label class="field"><span>Max dx</span><input id="max-dx" type="number" step="any" min="1e-12" value="${axis.max ?? 0.02}"></label>
      </div>
      <label class="field"><span>Cluster</span>
        <select id="max-min-side">
          <option value="end">Lower end</option>
          <option value="start">Upper end</option>
          <option value="both">Both ends</option>
          <option value="middle">Middle</option>
        </select>
      </label>
    `;
    bindProfileSelect(axis);
    document.getElementById("max-min-side").value = axis.side ?? "end";
    for (const id of ["min-dx", "max-dx"]) {
      document.getElementById(id).oninput = () => {
        axis[id === "min-dx" ? "min" : "max"] = positiveNumber(document.getElementById(id).value, 0.01);
        renderSpacingPreview();
        renderScene();
      };
    }
    document.getElementById("max-min-side").onchange = (event) => {
      axis.side = event.target.value;
      renderSpacingPreview();
      renderScene();
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
      renderSpacingPreview();
      renderScene();
    };
    document.getElementById("snac-gr").oninput = (event) => {
      axis.gr = Math.max(0, Number(event.target.value) || 0);
      renderSpacingPreview();
      renderScene();
    };
  }
  if (window.lucide) window.lucide.createIcons();
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
      clearCheckState();
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
      clearCheckState();
    };
  }
}

function renderCheckResults() {
  if (!els.checkResults) return;
  const errors = lastCheck?.errors ?? [];
  const warnings = lastCheck?.warnings ?? [];
  els.exportCase.disabled = false;
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
  blockGroup.clear();
  selectedGridGroup.clear();
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
    blockGroup.add(mesh);

    const edges = new THREE.LineSegments(
      new THREE.EdgesGeometry(mesh.geometry),
      new THREE.LineBasicMaterial({ color: block.id === selectedId ? 0xffffff : 0x9aa6b2, transparent: true, opacity: block.id === selectedId ? 0.95 : 0.48 })
    );
    edges.position.copy(mesh.position);
    blockGroup.add(edges);
  }
  drawSelectedGrid();
  updateAxesHelper();
}

function updateAxesHelper() {
  axesGroup.clear();
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
  if (axisLabelMaterials.has(key)) return new THREE.Sprite(axisLabelMaterials.get(key));
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
  return new THREE.Sprite(material);
}

function drawSelectedGrid() {
  selectedGridGroup.clear();
  const block = selectedBlock();
  if (!block) return;
  const material = new THREE.LineBasicMaterial({ color: 0xffffff, transparent: true, opacity: 0.32 });
  const [xmin, ymin, zmin] = block.lmin;
  const [xmax, ymax, zmax] = block.lmax;
  for (const axis of AXES) {
    const idx = AXES.indexOf(axis);
    const faces = axisFaces(block, axis);
    const step = Math.max(1, Math.floor(faces.length / 34));
    for (let i = 0; i < faces.length; i += step) {
      const value = faces[i];
      const points = [];
      if (idx === 0) {
        points.push(new THREE.Vector3(value, ymin, zmin), new THREE.Vector3(value, ymax, zmin), new THREE.Vector3(value, ymax, zmax), new THREE.Vector3(value, ymin, zmax), new THREE.Vector3(value, ymin, zmin));
      } else if (idx === 1) {
        points.push(new THREE.Vector3(xmin, value, zmin), new THREE.Vector3(xmax, value, zmin), new THREE.Vector3(xmax, value, zmax), new THREE.Vector3(xmin, value, zmax), new THREE.Vector3(xmin, value, zmin));
      } else {
        points.push(new THREE.Vector3(xmin, ymin, value), new THREE.Vector3(xmax, ymin, value), new THREE.Vector3(xmax, ymax, value), new THREE.Vector3(xmin, ymax, value), new THREE.Vector3(xmin, ymin, value));
      }
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
  const faces = axisFaces(block, currentAxis);
  const widths = [];
  for (let i = 1; i < faces.length; i += 1) widths.push(faces[i] - faces[i - 1]);
  const min = Math.min(...widths);
  const max = Math.max(...widths);
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
  els.spacingStats.innerHTML = `
    <span>min<strong>${formatNumber(min)}</strong></span>
    <span>max<strong>${formatNumber(max)}</strong></span>
    <span>ratio<strong>${formatNumber(max / min)}</strong></span>
  `;
}

function axisFaces(block, axisName) {
  const idx = AXES.indexOf(axisName);
  const n = Math.max(1, Math.round(block.ng[idx]));
  const lmin = block.lmin[idx];
  const lmax = block.lmax[idx];
  const axis = ensureAxis(block, axisName);
  if (axis.kind === "multi" || axis.kind === "simple_ratio") {
    const segments = segmentsFromAxis(axis, lmax - lmin);
    const widths = multiWidths(lmax - lmin, n, segments, axis.profile ?? "geometric");
    const faces = [lmin];
    for (const width of widths) faces.push(faces[faces.length - 1] + width);
    faces[faces.length - 1] = lmax;
    return faces;
  }
  if (axis.kind === "max_min") {
    const segments = segmentsFromAxis(axis, lmax - lmin);
    const widths = multiWidths(lmax - lmin, n, segments, axis.profile ?? "geometric");
    const faces = [lmin];
    for (const width of widths) faces.push(faces[faces.length - 1] + width);
    faces[faces.length - 1] = lmax;
    return faces;
  }
  const mapper = snacMapper(axis.gt ?? 0);
  const faces = [lmin];
  for (let q = 1; q <= n; q += 1) {
    const r0 = q / n;
    faces.push(lmin + mapper(axis.gr ?? 0, r0) * (lmax - lmin));
  }
  faces[faces.length - 1] = lmax;
  return faces;
}

function segmentsFromAxis(axis, length) {
  if (axis.kind === "multi") return axis.segments?.length ? axis.segments : [{ length: 1, cells: 1, ratio: 1 }];
  if (axis.kind === "simple_ratio") {
    const ratio = positiveNumber(axis.ratio, 1);
    return [{ length: 1, cells: 1, ratio: axis.side === "start" ? 1 / ratio : ratio }];
  }
  if (axis.kind === "max_min") {
    const min = positiveNumber(axis.min, length);
    const max = positiveNumber(axis.max, length);
    const ratio = Math.max(min, max) / Math.min(min, max);
    if (axis.side === "start") return [{ length: 1, cells: 1, ratio: 1 / ratio }];
    if (axis.side === "both") return [{ length: 0.5, cells: 0.5, ratio }, { length: 0.5, cells: 0.5, ratio: 1 / ratio }];
    if (axis.side === "middle") return [{ length: 0.5, cells: 0.5, ratio: 1 / ratio }, { length: 0.5, cells: 0.5, ratio }];
    return [{ length: 1, cells: 1, ratio }];
  }
  return [{ length: 1, cells: 1, ratio: 1 }];
}

function multiWidths(length, n, segments, profile = "geometric") {
  const lengthSum = segments.reduce((sum, item) => sum + positiveNumber(item.length, 1), 0);
  const cellWeights = segments.map((item) => positiveNumber(item.cells, 1));
  const cellSum = cellWeights.reduce((sum, value) => sum + value, 0);
  let cells = cellWeights.map((value) => Math.max(1, Math.floor((value / cellSum) * n)));
  while (cells.reduce((sum, value) => sum + value, 0) > n) {
    const index = cells.findIndex((value) => value > 1);
    if (index < 0) break;
    cells[index] -= 1;
  }
  while (cells.reduce((sum, value) => sum + value, 0) < n) {
    let best = 0;
    let bestRem = -1;
    for (let i = 0; i < cellWeights.length; i += 1) {
      const rem = (cellWeights[i] / cellSum) * n - Math.floor((cellWeights[i] / cellSum) * n);
      if (rem > bestRem) {
        best = i;
        bestRem = rem;
      }
    }
    cells[best] += 1;
  }
  const widths = [];
  segments.forEach((segment, index) => {
    const segLength = (positiveNumber(segment.length, 1) / lengthSum) * length;
    widths.push(...profileWidths(segLength, cells[index], positiveNumber(segment.ratio, 1), profile));
  });
  widths[widths.length - 1] += length - widths.reduce((sum, value) => sum + value, 0);
  return widths;
}

function profileWidths(length, n, ratio, profile) {
  if (profile === "geometric") return geometricWidths(length, n, ratio);
  if (profile === "tanh" || profile === "erf") return mappedRatioWidths(length, n, ratio, profile);
  return geometricWidths(length, n, ratio);
}

function geometricWidths(length, n, ratio) {
  if (n <= 1) return [length];
  const cellRatio = Math.pow(ratio, 1 / (n - 1));
  const first = Math.abs(cellRatio - 1) < 1e-12 ? length / n : (length * (1 - cellRatio)) / (1 - Math.pow(cellRatio, n));
  const widths = Array.from({ length: n }, (_, index) => first * Math.pow(cellRatio, index));
  widths[widths.length - 1] += length - widths.reduce((sum, value) => sum + value, 0);
  return widths;
}

function mappedRatioWidths(length, n, ratio, profile) {
  if (n <= 1 || Math.abs(ratio - 1) < 1e-12) return Array.from({ length: n }, () => length / n);
  const target = Math.max(ratio, 1 / ratio);
  const alpha = solveMappingAlpha(n, target, profile);
  const faces = Array.from({ length: n + 1 }, (_, index) => mapOneEnd(alpha, index / n, profile));
  let widths = Array.from({ length: n }, (_, index) => (faces[index + 1] - faces[index]) * length);
  if (ratio < 1) widths = widths.reverse();
  widths[widths.length - 1] += length - widths.reduce((sum, value) => sum + value, 0);
  return widths;
}

function solveMappingAlpha(n, target, profile) {
  let high = 1;
  while (mappingRatio(n, high, profile) < target && high < 80) high *= 2;
  if (mappingRatio(n, high, profile) < target) return high;
  let low = 0;
  for (let step = 0; step < 80; step += 1) {
    const mid = 0.5 * (low + high);
    if (mappingRatio(n, mid, profile) < target) low = mid;
    else high = mid;
  }
  return 0.5 * (low + high);
}

function mappingRatio(n, alpha, profile) {
  const faces = Array.from({ length: n + 1 }, (_, index) => mapOneEnd(alpha, index / n, profile));
  const start = faces[1] - faces[0];
  const end = faces[n] - faces[n - 1];
  return end / start;
}

function mapOneEnd(alpha, r0, profile) {
  if (!alpha) return r0;
  if (profile === "erf") return 1 + erf((r0 - 1) * alpha) / erf(alpha);
  return 1 + Math.tanh((r0 - 1) * alpha) / Math.tanh(alpha);
}

function erf(value) {
  const sign = value < 0 ? -1 : 1;
  const x = Math.abs(value);
  const a1 = 0.254829592;
  const a2 = -0.284496736;
  const a3 = 1.421413741;
  const a4 = -1.453152027;
  const a5 = 1.061405429;
  const p = 0.3275911;
  const t = 1 / (1 + p * x);
  const y = 1 - (((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t) * Math.exp(-x * x);
  return sign * y;
}

function snacMapper(gt) {
  const functions = [
    (alpha, r0) => (alpha ? 0.5 * (1 + Math.tanh((r0 - 0.5) * alpha) / Math.tanh(alpha / 2)) : r0),
    (alpha, r0) => (alpha ? 1 + Math.tanh((r0 - 1) * alpha) / Math.tanh(alpha) : r0),
    (alpha, r0) => {
      if (!alpha) return r0;
      if (r0 <= 0.5) return 0.5 * Math.tanh(2 * alpha * r0) / Math.tanh(alpha);
      return 0.5 * (2 + Math.tanh(2 * alpha * (r0 - 1)) / Math.tanh(alpha));
    },
    (alpha, r0) => (alpha ? 1 - (1 + Math.tanh(((1 - r0) - 1) * alpha) / Math.tanh(alpha)) : r0),
    geometricOneEnd,
    (alpha, r0) => 1 - geometricOneEnd(alpha, 1 - r0),
    (alpha, r0) => (r0 <= 0.5 ? 0.5 * geometricOneEnd(alpha / 2, 2 * r0) : 1 - 0.5 * geometricOneEnd(alpha / 2, 2 * (1 - r0))),
    (alpha, r0) => (r0 <= 0.5 ? 0.5 * (1 - geometricOneEnd(alpha / 2, (0.5 - r0) * 2)) : 0.5 * (1 + geometricOneEnd(alpha / 2, (r0 - 0.5) * 2))),
  ];
  return functions[gt] ?? functions[0];
}

function geometricOneEnd(alpha, r0) {
  if (r0 === 1 || !alpha) return r0;
  return r0 * Math.pow((1 - Math.pow(r0, alpha)) / (1 - r0) / alpha, 1 / 3);
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
  renderAll();
}

function removeSelectedBlock() {
  if (!selectedBlock()) return;
  project.blocks = project.blocks.filter((block) => block.id !== selectedId);
  selectedId = project.blocks[0]?.id ?? null;
  renderAll();
}

function resetProject() {
  project.blocks = [];
  selectedId = null;
  lastCheck = null;
  renderAll();
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

async function updateStructure() {
  setStatus("Updating structure...");
  try {
    const payload = await postJson("/api/update", { project, sourceBlockId: selectedId }, { allowInvalid: true });
    if (payload.project) {
      const oldSelectedId = selectedId;
      project = payload.project;
      selectedId = project.blocks.find((block) => block.id === oldSelectedId)?.id ?? project.blocks[0]?.id ?? null;
    }
    lastCheck = payload;
    renderAll();
    const warningText = payload.warnings?.length ? `, ${payload.warnings.length} warnings` : "";
    setStatus(payload.ok ? `Structure updated${warningText}` : `${payload.errors.length} check issue${payload.errors.length === 1 ? "" : "s"}`);
  } catch (error) {
    setStatus(error.message);
  }
}

async function exportCase() {
  const valid = await checkProject({ silent: true });
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
      renderAll();
      setStatus(`Loaded ${file.name}`);
    } catch (error) {
      setStatus(error.message);
    }
  };
  reader.readAsText(file);
  event.target.value = "";
}

function selectFromScene(event) {
  const rect = renderer.domElement.getBoundingClientRect();
  pointer.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
  pointer.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;
  raycaster.setFromCamera(pointer, camera);
  const hit = raycaster.intersectObjects(blockGroup.children, false).find((item) => item.object.userData.blockId);
  if (hit) {
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
  block.axes[axis].segments = block.axes[axis].segments || [{ length: 1, cells: 1, ratio: 1 }];
  block.axes[axis].profile = block.axes[axis].profile || "geometric";
  return block.axes[axis];
}

function defaultAxis() {
  return {
    kind: "snac",
    gt: 0,
    gr: 0,
    ratio: 1,
    profile: "geometric",
    min: 0.01,
    max: 0.02,
    side: "end",
    segments: [{ length: 1, cells: 1, ratio: 1 }],
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

function numberInput(field, index, value, min = null) {
  const minAttr = min == null ? "" : `min="${min}"`;
  return `<label class="field"><input data-field="${field}" data-index="${index}" type="number" step="any" ${minAttr} value="${value}"></label>`;
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

function setStatus(message) {
  els.status.textContent = message;
}

function clearCheckState() {
  lastCheck = null;
  renderCheckResults();
  setStatus("");
}

function escapeHtml(value) {
  return String(value).replace(/[&<>"']/g, (char) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#039;" })[char]);
}

function toCamel(id) {
  return id.replace(/-([a-z])/g, (_, char) => char.toUpperCase());
}
