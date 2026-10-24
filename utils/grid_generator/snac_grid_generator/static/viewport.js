// The three.js scene: drawing, inspection tools, clipping, transforms, and picking.

import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { TransformControls } from "three/addons/controls/TransformControls.js";
import { renderBlockList, selectBlock } from "./block-list.js";
import { AXES, FACE_ORDER, blockCenter, blockSize, ensureAxis, rescaleExplicitFaces } from "./block-model.js";
import { handleCameraChange, resizeRenderer } from "./camera.js";
import { markProjectDirty, renderAll, setStatus } from "./core.js";
import { commitHistory } from "./history.js";
import { updateGeometryInputValues } from "./inspector.js";
import { scheduleGridPreviews, schedulePreview } from "./previews.js";
import {
  axisLabel,
  clearDisposableGroup,
  partitionFaceIndices,
  planeOutlinePoints,
  projectBox,
} from "./scene-helpers.js";
import {
  COLORS,
  els,
  previewCache,
  selectedBlock,
  selectedBlockIds,
  selectedVisibleBlocks,
  state,
  ui,
  visibleBlocks,
} from "./state.js";

export const NOTICEABLE_SPACING_JUMP = 1.5;
export const BOUNDARY_COLORS = { D: 0xe46868, F: 0x43c989, N: 0x76a9ff, mixed: 0xd6bf5b };

export const clippingPlane = new THREE.Plane(new THREE.Vector3(1, 0, 0), 0);

export function initScene() {
  ui.scene = new THREE.Scene();
  ui.scene.background = new THREE.Color(0x101216);
  ui.perspectiveCamera = new THREE.PerspectiveCamera(45, 1, 0.01, 10000);
  ui.perspectiveCamera.position.set(3.0, -4.5, 3.0);
  ui.orthographicCamera = new THREE.OrthographicCamera(-1, 1, 1, -1, 0.01, 10000);
  ui.camera = ui.perspectiveCamera;
  els.scene.dataset.projection = "perspective";
  els.scene.dataset.view = "3d";
  ui.renderer = new THREE.WebGLRenderer({ canvas: els.scene, antialias: true, alpha: false });
  ui.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  ui.renderer.localClippingEnabled = true;
  ui.controls = new OrbitControls(ui.camera, ui.renderer.domElement);
  ui.controls.enableDamping = true;
  ui.controls.dampingFactor = 0.08;
  ui.controls.screenSpacePanning = true;
  ui.controls.addEventListener("change", handleCameraChange);
  ui.raycaster = new THREE.Raycaster();
  ui.pointer = new THREE.Vector2();
  ui.blockGroup = new THREE.Group();
  ui.boundaryFaceGroup = new THREE.Group();
  ui.gridLineGroup = new THREE.Group();
  ui.spacingJumpGroup = new THREE.Group();
  ui.diagnosticGroup = new THREE.Group();
  ui.partitionGroup = new THREE.Group();
  ui.axesGroup = new THREE.Group();
  ui.gridHelper = new THREE.GridHelper(10, 20, 0x46515c, 0x252b32);
  ui.gridHelper.rotation.x = Math.PI / 2;
  ui.scene.add(ui.gridHelper);
  ui.scene.add(ui.axesGroup);
  ui.scene.add(ui.blockGroup);
  ui.scene.add(ui.boundaryFaceGroup);
  ui.scene.add(ui.gridLineGroup);
  ui.scene.add(ui.spacingJumpGroup);
  ui.scene.add(ui.diagnosticGroup);
  ui.scene.add(ui.partitionGroup);

  ui.transformControls = new TransformControls(ui.camera, ui.renderer.domElement);
  ui.transformControls.setSpace("world");
  ui.transformControls.setSize(0.82);
  ui.transformControls.addEventListener("dragging-changed", handleTransformDragging);
  ui.transformControls.addEventListener("objectChange", handleTransformObjectChange);
  ui.scene.add(ui.transformControls.getHelper ? ui.transformControls.getHelper() : ui.transformControls);

  const ambient = new THREE.HemisphereLight(0xffffff, 0x1a2028, 1.9);
  const key = new THREE.DirectionalLight(0xffffff, 1.4);
  key.position.set(4, -5, 6);
  ui.scene.add(ambient, key);
  resizeRenderer();
}

export function renderScene() {
  ui.transformControls.detach();
  ui.selectedBlockMesh = null;
  ui.boundaryFaceGroup.visible = true;
  ui.gridLineGroup.visible = ui.gridLinesMode !== "off";
  ui.partitionGroup.visible = true;
  ui.spacingJumpGroup.visible = ui.spacingJumpsVisible;
  clearDisposableGroup(ui.blockGroup);
  clearDisposableGroup(ui.boundaryFaceGroup);
  clearDisposableGroup(ui.gridLineGroup);
  clearDisposableGroup(ui.partitionGroup);
  clearDisposableGroup(ui.spacingJumpGroup);
  clearDisposableGroup(ui.diagnosticGroup);
  const diagnosticIds = new Set(
    (state.activeDiagnostic?.locations ?? []).map((location) => location.blockId)
  );
  const selectedSet = new Set(selectedBlockIds());
  for (const block of state.project.blocks) {
    if (ui.hiddenBlockIds.has(block.id)) continue;
    const size = blockSize(block);
    const center = blockCenter(block);
    const color = COLORS[(block.id - 1) % COLORS.length];
    const material = new THREE.MeshStandardMaterial({
      color,
      transparent: true,
      opacity: diagnosticIds.has(block.id)
        ? 0.52
        : block.id === state.selectedId ? 0.43 : selectedSet.has(block.id) ? 0.34 : 0.24,
      roughness: 0.62,
      metalness: 0.0,
      side: THREE.DoubleSide,
    });
    const mesh = new THREE.Mesh(new THREE.BoxGeometry(size[0], size[1], size[2]), material);
    mesh.position.set(center[0], center[1], center[2]);
    mesh.userData.blockId = block.id;
    mesh.userData.baseSize = size;
    if (block.id === state.selectedId) ui.selectedBlockMesh = mesh;
    ui.blockGroup.add(mesh);

    const edges = new THREE.LineSegments(
      new THREE.EdgesGeometry(mesh.geometry),
      new THREE.LineBasicMaterial({
        color: diagnosticIds.has(block.id)
          ? state.activeDiagnostic.severity === "error" ? 0xff6b6b : 0xf0c36a
          : block.id === state.selectedId ? 0xffffff : selectedSet.has(block.id) ? color : 0x9aa6b2,
        transparent: true,
        opacity: diagnosticIds.has(block.id) || selectedSet.has(block.id) ? 0.95 : 0.48,
      })
    );
    mesh.add(edges);
  }
  drawBoundaryFaces();
  drawGridLines();
  drawPartitionPlanes();
  drawInterfaceSpacing();
  drawDiagnosticHighlight();
  updateAxesHelper();
  updateClippingPlane();
  applyClipping();
  attachSelectedTransform();
  renderViewControls();
}

export function drawDiagnosticHighlight() {
  if (!state.activeDiagnostic) return;
  const color = state.activeDiagnostic.severity === "error" ? 0xff6b6b : 0xf0c36a;
  for (const location of state.activeDiagnostic.locations ?? []) {
    if (!location.face) continue;
    if (ui.hiddenBlockIds.has(location.blockId)) continue;
    const block = state.project.blocks.find((item) => item.id === location.blockId);
    const face = FACE_ORDER.indexOf(location.face);
    if (!block || face < 0) continue;
    const axisIndex = Math.floor(face / 2);
    const value = face % 2 === 0 ? block.lmin[axisIndex] : block.lmax[axisIndex];
    const line = new THREE.Line(
      new THREE.BufferGeometry().setFromPoints(planeOutlinePoints(block, axisIndex, value)),
      new THREE.LineBasicMaterial({ color, depthTest: false })
    );
    line.renderOrder = 8;
    ui.diagnosticGroup.add(line);
  }
}

export function drawInterfaceSpacing() {
  if (!ui.spacingJumpsVisible) return;
  for (const metric of state.lastCheck?.interfaces ?? []) {
    if (ui.hiddenBlockIds.has(metric.blockId) || ui.hiddenBlockIds.has(metric.neighborId)) continue;
    const block = state.project.blocks.find((item) => item.id === metric.blockId);
    if (!block) continue;
    const axisIndex = AXES.indexOf(metric.axis);
    const value = metric.face.endsWith("+") ? block.lmax[axisIndex] : block.lmin[axisIndex];
    // Red follows the server's warning threshold; amber marks smaller visible jumps.
    const color = metric.exceedsWarning ? 0xe46868 : metric.ratio > NOTICEABLE_SPACING_JUMP ? 0xd6bf5b : 0x43c989;
    const line = new THREE.Line(
      new THREE.BufferGeometry().setFromPoints(planeOutlinePoints(block, axisIndex, value)),
      new THREE.LineBasicMaterial({ color, transparent: true, opacity: 0.95, depthTest: false })
    );
    line.renderOrder = 5;
    ui.spacingJumpGroup.add(line);
  }
}

export function drawBoundaryFaces() {
  const drawnBlockIds = [];
  if (!ui.boundaryColorsVisible) {
    els.scene.dataset.boundaryBlocks = "";
    return;
  }
  for (const block of selectedVisibleBlocks()) {
    drawnBlockIds.push(block.id);
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
      ui.boundaryFaceGroup.add(mesh);
    }
  }
  els.scene.dataset.boundaryBlocks = drawnBlockIds.join(",");
}

export function boundaryFaceColor(block, face) {
  const velocityCodes = block.cbcvel.map((component) => component[face]);
  const allCodes = [block.cbcpre[face], ...velocityCodes, ...block.cbcscal.map((scalar) => scalar[face])];
  if (allCodes.every((code) => code === "F")) return BOUNDARY_COLORS.F;
  const velocityCode = velocityCodes[0];
  if (["D", "N"].includes(velocityCode) && velocityCodes.every((code) => code === velocityCode)) {
    return BOUNDARY_COLORS[velocityCode];
  }
  return BOUNDARY_COLORS.mixed;
}

export function drawPartitionPlanes() {
  const drawnBlockIds = [];
  if (!ui.partitionPlanesVisible) {
    els.scene.dataset.partitionBlocks = "";
    return;
  }
  for (const block of selectedVisibleBlocks()) {
    const preview = previewCache.get(block.id);
    if (!preview) continue;
    drawnBlockIds.push(block.id);
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
        ui.partitionGroup.add(line);
      }
    }
  }
  els.scene.dataset.partitionBlocks = drawnBlockIds.join(",");
}

export function updateAxesHelper() {
  clearDisposableGroup(ui.axesGroup, false);
  const box = projectBox(visibleBlocks());
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
    ui.axesGroup.add(arrow, sprite);
  }
}

export function drawGridLines() {
  const drawnBlockIds = [];
  if (ui.gridLinesMode === "off") {
    els.scene.dataset.gridBlocks = "";
    return;
  }
  const blocks = ui.gridLinesMode === "all" ? visibleBlocks() : selectedVisibleBlocks();
  for (const block of blocks) {
    const preview = previewCache.get(block.id);
    if (!preview) continue;
    drawnBlockIds.push(block.id);
    const material = new THREE.LineBasicMaterial({
      color: block.id === state.selectedId ? 0xffffff : COLORS[(block.id - 1) % COLORS.length],
      transparent: true,
      opacity: block.id === state.selectedId ? 0.36 : 0.28,
    });
    for (const axis of AXES) {
      const idx = AXES.indexOf(axis);
      const faces = preview[axis]?.faces ?? [];
      const step = Math.max(1, Math.floor(faces.length / 34));
      for (let i = 0; i < faces.length; i += step) {
        const points = planeOutlinePoints(block, idx, faces[i]);
        ui.gridLineGroup.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(points), material));
      }
    }
  }
  els.scene.dataset.gridBlocks = drawnBlockIds.join(",");
}

export function hideSelectedBlocks() {
  const blockIds = selectedBlockIds();
  if (!blockIds.length) return;
  for (const blockId of blockIds) ui.hiddenBlockIds.add(blockId);
  refreshInspectionView(`Hidden ${blockIds.length} block${blockIds.length === 1 ? "" : "s"}`);
}

export function isolateSelectedBlocks() {
  const blockIds = new Set(selectedBlockIds());
  if (!blockIds.size) return;
  ui.hiddenBlockIds = new Set(
    state.project.blocks.map((block) => block.id).filter((blockId) => !blockIds.has(blockId))
  );
  refreshInspectionView(`Isolated ${blockIds.size} block${blockIds.size === 1 ? "" : "s"}`);
}

export function showAllBlocks() {
  if (!ui.hiddenBlockIds.size) return;
  ui.hiddenBlockIds.clear();
  refreshInspectionView("Showing all blocks");
}

export function refreshInspectionView(message) {
  renderBlockList();
  renderScene();
  scheduleGridPreviews(0);
  setStatus(message);
}

export function resetInspectionState() {
  ui.hiddenBlockIds.clear();
  ui.clippingEnabled = false;
  ui.clippingAxisIndex = 0;
  ui.clippingFraction = 0.5;
  ui.clippingDirection = 1;
}

export function toggleClipping() {
  ui.clippingEnabled = !ui.clippingEnabled;
  els.clipTools.hidden = !ui.clippingEnabled;
  updateClippingPlane();
  applyClipping();
  renderViewControls();
}

export function setClippingAxis(axis) {
  const axisIndex = AXES.indexOf(axis);
  if (axisIndex < 0) return;
  ui.clippingAxisIndex = axisIndex;
  updateClippingPlane();
  renderViewControls();
}

export function updateClippingPlane() {
  const box = projectBox(visibleBlocks());
  const coordinate = box.min.getComponent(ui.clippingAxisIndex)
    + ui.clippingFraction * (box.max.getComponent(ui.clippingAxisIndex) - box.min.getComponent(ui.clippingAxisIndex));
  const normal = new THREE.Vector3().setComponent(ui.clippingAxisIndex, ui.clippingDirection);
  clippingPlane.set(normal, -ui.clippingDirection * coordinate);
  updateInspectionDataset();
}

export function applyClipping() {
  const planes = ui.clippingEnabled ? [clippingPlane] : null;
  for (const group of [
    ui.blockGroup,
    ui.boundaryFaceGroup,
    ui.gridLineGroup,
    ui.partitionGroup,
    ui.spacingJumpGroup,
    ui.diagnosticGroup,
  ]) {
    group.traverse((object) => {
      if (!object.material) return;
      for (const material of Array.isArray(object.material) ? object.material : [object.material]) {
        material.clippingPlanes = planes;
        material.needsUpdate = true;
      }
    });
  }
  updateInspectionDataset();
}

export function updateInspectionDataset() {
  els.scene.dataset.clipping = String(ui.clippingEnabled);
  els.scene.dataset.clipAxis = AXES[ui.clippingAxisIndex];
  els.scene.dataset.clipDirection = String(ui.clippingDirection);
  els.scene.dataset.hiddenBlocks = [...ui.hiddenBlockIds].sort((a, b) => a - b).join(",");
  els.scene.dataset.selectedBlocks = selectedBlockIds().join(",");
  els.scene.dataset.primaryBlock = state.selectedId == null ? "" : String(state.selectedId);
  els.scene.dataset.gridMode = ui.gridLinesMode;
}

export function setBlockTool(tool) {
  ui.blockTool = tool === ui.blockTool && tool !== "select" ? "select" : tool;
  renderScene();
}

export function attachSelectedTransform() {
  if (!ui.selectedBlockMesh || ui.blockTool === "select") return;
  ui.transformControls.setMode(ui.blockTool);
  ui.transformControls.attach(ui.selectedBlockMesh);
}

export function handleTransformDragging(event) {
  if (event.value) {
    if (!ui.selectedBlockMesh) return;
    commitHistory();
    ui.isTransformDragging = true;
    ui.transformDrag = {
      blockId: ui.selectedBlockMesh.userData.blockId,
      baseSize: [...ui.selectedBlockMesh.userData.baseSize],
      baseLmin: [...selectedBlock().lmin],
      baseLmax: [...selectedBlock().lmax],
      baseFaces: AXES.map((axis) => {
        const spec = ensureAxis(selectedBlock(), axis);
        return spec.kind === "explicit" ? [...spec.faces] : null;
      }),
    };
    ui.controls.enabled = false;
    ui.boundaryFaceGroup.visible = false;
    ui.gridLineGroup.visible = false;
    ui.partitionGroup.visible = false;
    return;
  }

  ui.controls.enabled = true;
  if (!ui.isTransformDragging || !ui.transformDrag) return;
  const blockId = ui.transformDrag.blockId;
  ui.isTransformDragging = false;
  ui.transformDrag = null;
  commitHistory();
  window.setTimeout(() => {
    renderAll();
    schedulePreview(blockId, 0);
  }, 0);
}

export function handleTransformObjectChange() {
  const mesh = ui.transformControls.object;
  if (!ui.isTransformDragging || !ui.transformDrag || !mesh) return;
  const block = state.project.blocks.find((item) => item.id === ui.transformDrag.blockId);
  if (!block) return;
  const center = mesh.position.toArray();
  const size = ui.transformDrag.baseSize.map((value, index) => Math.max(1.0e-9, value * Math.abs(mesh.scale.getComponent(index))));
  block.lmin = center.map((value, index) => value - 0.5 * size[index]);
  block.lmax = center.map((value, index) => value + 0.5 * size[index]);
  for (let index = 0; index < 3; index += 1) {
    rescaleExplicitFaces(
      block,
      index,
      ui.transformDrag.baseLmin[index],
      ui.transformDrag.baseLmax[index],
      block.lmin[index],
      block.lmax[index],
      ui.transformDrag.baseFaces[index]
    );
  }
  updateGeometryInputValues(block);
  if (!ui.transformDrag.dirty) {
    // Clearing check, repair, and decomposition state once per drag is enough.
    markProjectDirty(block.id);
    ui.transformDrag.dirty = true;
  } else {
    previewCache.delete(block.id);
  }
  renderBlockList();
}

export function renderViewControls() {
  const selectedVisible = selectedVisibleBlocks();
  if (!selectedVisible.length && ui.blockTool !== "select") ui.blockTool = "select";
  const toolButtons = {
    select: els.toolSelect,
    translate: els.toolMove,
    scale: els.toolScale,
  };
  for (const [tool, button] of Object.entries(toolButtons)) {
    button.classList.toggle("active", ui.blockTool === tool);
    button.setAttribute("aria-pressed", String(ui.blockTool === tool));
  }
  els.toolMove.disabled = !selectedVisible.length;
  els.toolScale.disabled = !selectedVisible.length;
  els.focusSelection.disabled = !selectedVisible.length;
  els.hideSelected.disabled = !selectedVisible.length;
  els.isolateSelected.disabled = !selectedBlockIds().length;
  els.showAllBlocks.disabled = !ui.hiddenBlockIds.size;
  els.toggleClipping.classList.toggle("active", ui.clippingEnabled);
  els.toggleClipping.setAttribute("aria-pressed", String(ui.clippingEnabled));
  els.clipTools.hidden = !ui.clippingEnabled;
  els.clipPosition.value = String(ui.clippingFraction);
  for (const button of els.clipAxis.querySelectorAll("button")) {
    const active = AXES[ui.clippingAxisIndex] === button.dataset.clipAxis;
    button.classList.toggle("active", active);
    button.setAttribute("aria-pressed", String(active));
  }
  els.flipClipping.classList.toggle("active", ui.clippingDirection < 0);
  els.flipClipping.setAttribute("aria-pressed", String(ui.clippingDirection < 0));
  els.toggleBoundaries.classList.toggle("active", ui.boundaryColorsVisible);
  els.toggleBoundaries.setAttribute("aria-pressed", String(ui.boundaryColorsVisible));
  els.togglePartitions.classList.toggle("active", ui.partitionPlanesVisible);
  els.togglePartitions.setAttribute("aria-pressed", String(ui.partitionPlanesVisible));
  els.gridLinesMode.value = ui.gridLinesMode;
  els.toggleSpacingJumps.classList.toggle("active", ui.spacingJumpsVisible);
  els.toggleSpacingJumps.setAttribute("aria-pressed", String(ui.spacingJumpsVisible));
  els.toggleGrid.classList.toggle("active", ui.groundGridVisible);
  els.toggleGrid.setAttribute("aria-pressed", String(ui.groundGridVisible));
}

export function beginSceneSelection(event) {
  if (ui.blockTool !== "select" || event.button !== 0) return;
  ui.scenePointerDown = {
    pointerId: event.pointerId,
    x: event.clientX,
    y: event.clientY,
    shiftKey: event.shiftKey,
  };
}

export function completeSceneSelection(event) {
  if (!ui.scenePointerDown || event.pointerId !== ui.scenePointerDown.pointerId) return;
  const pointerDown = ui.scenePointerDown;
  ui.scenePointerDown = null;
  if (Math.hypot(event.clientX - pointerDown.x, event.clientY - pointerDown.y) > 4) return;
  selectFromScene(event, pointerDown.shiftKey || event.shiftKey);
}

export function cancelSceneSelection() {
  ui.scenePointerDown = null;
}

export function selectFromScene(event, additive = false) {
  if (ui.blockTool !== "select") return;
  const rect = ui.renderer.domElement.getBoundingClientRect();
  ui.pointer.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
  ui.pointer.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;
  ui.raycaster.setFromCamera(ui.pointer, ui.camera);
  const hit = ui.raycaster.intersectObjects(ui.blockGroup.children, false).find((item) => item.object.userData.blockId);
  if (hit) selectBlock(hit.object.userData.blockId, additive);
}

export function animate() {
  requestAnimationFrame(animate);
  ui.controls.update();
  ui.renderer.render(ui.scene, ui.camera);
}
