// The three.js scene: drawing, inspection tools, clipping, transforms, and picking.

import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { TransformControls } from "three/addons/controls/TransformControls.js";
import { renderBlockList, selectBlock } from "./block-list.js";
import { AXES, FACE_ORDER, blockCenter, blockSize, copyBlock, ensureAxis, rescaleExplicitFaces } from "./block-model.js";
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
  blockLabel,
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

// Shift snaps a dragged face to another block's face within this many pixels.
const SNAP_PIXELS = 10;
const SNAP_COLOR = 0x5ce1e6;
// A drag never shrinks a block below this fraction of its size when the drag started.
const MIN_DRAG_FRACTION = 1e-3;

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
  ui.snapGroup = new THREE.Group();
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
  ui.scene.add(ui.snapGroup);

  ui.transformControls = new TransformControls(ui.camera, ui.renderer.domElement);
  ui.transformControls.setSpace("world");
  ui.transformControls.setSize(0.82);
  ui.transformControls.addEventListener("dragging-changed", handleTransformDragging);
  ui.transformControls.addEventListener("objectChange", handleTransformObjectChange);
  ui.scene.add(ui.transformControls.getHelper ? ui.transformControls.getHelper() : ui.transformControls);
  // Capture listeners run before the transform controls see the same pointer event.
  for (const type of ["pointerdown", "pointermove"]) {
    ui.renderer.domElement.addEventListener(type, trackTransformModifiers, { capture: true });
  }
  window.addEventListener("keydown", trackTransformModifiers);
  window.addEventListener("keyup", trackTransformModifiers);
  window.addEventListener("blur", () => trackTransformModifiers({ type: "blur", altKey: false, shiftKey: false }));

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
    const mesh = ui.selectedBlockMesh;
    const block = selectedBlock();
    ui.transformDrag = {
      blockId: mesh.userData.blockId,
      mode: ui.transformControls.mode,
      axis: ui.transformControls.axis ?? "",
      // The handle grabbed on each axis decides which face a one-sided resize moves.
      side: ui.transformControls.pointStart.toArray().map((value) => (value < 0 ? -1 : 1)),
      baseSize: [...mesh.userData.baseSize],
      baseCenter: mesh.position.toArray(),
      baseLmin: [...block.lmin],
      baseLmax: [...block.lmax],
      baseFaces: AXES.map((axis) => {
        const spec = ensureAxis(block, axis);
        return spec.kind === "explicit" ? [...spec.faces] : null;
      }),
      travel: [0, 0, 0],
      ratio: [1, 1, 1],
      offset: [0, 0, 0],
      ghost: null,
    };
    if (ui.transformDrag.mode === "translate") {
      // Alt-dragging leaves a copy behind; the ghost shows where while Alt is held.
      ui.transformDrag.ghost = mesh.clone();
      ui.transformDrag.ghost.visible = false;
      ui.blockGroup.add(ui.transformDrag.ghost);
    }
    ui.controls.enabled = false;
    ui.boundaryFaceGroup.visible = false;
    ui.gridLineGroup.visible = false;
    ui.partitionGroup.visible = false;
    return;
  }

  ui.controls.enabled = true;
  if (!ui.isTransformDragging || !ui.transformDrag) return;
  const drag = ui.transformDrag;
  ui.isTransformDragging = false;
  ui.transformDrag = null;
  drawSnapGuides(null, []);
  const copy = drag.mode === "translate" && ui.transformModifiers.alt ? duplicateDraggedBlock(drag) : null;
  commitHistory();
  window.setTimeout(() => {
    renderAll();
    schedulePreview(drag.blockId, 0);
    if (!copy) return;
    schedulePreview(copy.id, 0);
    const source = state.project.blocks.find((block) => block.id === drag.blockId);
    if (source) setStatus(`Duplicated ${blockLabel(source)} as ${blockLabel(copy)}`);
  }, 0);
}

export function handleTransformObjectChange() {
  const mesh = ui.transformControls.object;
  const drag = ui.transformDrag;
  if (!ui.isTransformDragging || !drag || !mesh) return;
  const controls = ui.transformControls;
  // How far the pointer moved along each axis, in world units.
  drag.travel = controls.pointEnd.toArray().map((value, index) => value - controls.pointStart.getComponent(index));
  drag.ratio = mesh.scale.toArray();
  drag.offset = mesh.position.toArray().map((value, index) => value - drag.baseCenter[index]);
  applyTransformDrag();
}

// Alt and Shift change a drag in progress as soon as they are pressed or released.
function trackTransformModifiers(event) {
  const modifiers = { alt: event.altKey, shift: event.shiftKey };
  const changed = modifiers.alt !== ui.transformModifiers.alt || modifiers.shift !== ui.transformModifiers.shift;
  ui.transformModifiers = modifiers;
  if (!ui.isTransformDragging) return;
  // Keeps a lone Alt from opening the browser menu bar on Windows and Linux.
  if (event.key === "Alt") event.preventDefault?.();
  if (changed && !event.type.startsWith("pointer")) applyTransformDrag();
}

// Applies the drag to the block: the gizmo reports the pointer motion, and the
// modifiers decide which faces follow it.
function applyTransformDrag() {
  const drag = ui.transformDrag;
  const mesh = ui.transformControls.object;
  const block = state.project.blocks.find((item) => item.id === drag?.blockId);
  if (!drag || !mesh || !block) return;
  const { alt, shift } = ui.transformModifiers;
  const bounds = drag.mode === "scale" ? resizedBounds(drag, alt) : movedBounds(drag);
  const snaps = shift ? snapBounds(drag, bounds, alt) : [];
  block.lmin = bounds.map(([lower]) => lower);
  block.lmax = bounds.map(([, upper]) => upper);
  for (let index = 0; index < 3; index += 1) {
    rescaleExplicitFaces(
      block,
      index,
      drag.baseLmin[index],
      drag.baseLmax[index],
      block.lmin[index],
      block.lmax[index],
      drag.baseFaces[index]
    );
    // Draw the block where it now is, not where the gizmo alone would put it.
    mesh.position.setComponent(index, 0.5 * (block.lmin[index] + block.lmax[index]));
    mesh.scale.setComponent(index, (block.lmax[index] - block.lmin[index]) / drag.baseSize[index]);
  }
  if (drag.ghost) drag.ghost.visible = alt;
  drawSnapGuides(block, snaps);
  updateGeometryInputValues(block);
  if (!drag.dirty) {
    // Clearing check, repair, and decomposition state once per drag is enough.
    markProjectDirty(block.id);
    drag.dirty = true;
  } else {
    previewCache.delete(block.id);
  }
  renderBlockList();
}

// Resizing moves the grabbed face with the pointer, or with Alt both faces of
// that axis symmetrically. The centre handle always scales about the centre.
function resizedBounds(drag, symmetric) {
  return AXES.map((axis, index) => {
    const lower = drag.baseLmin[index];
    const upper = drag.baseLmax[index];
    const minimum = MIN_DRAG_FRACTION * (upper - lower);
    if (drag.axis === "XYZ") {
      const center = 0.5 * (lower + upper);
      const half = 0.5 * Math.max(minimum, (upper - lower) * Math.abs(drag.ratio[index]));
      return [center - half, center + half];
    }
    if (!dragsAxis(drag, index)) return [lower, upper];
    const growth = drag.travel[index] * drag.side[index];
    if (symmetric) return limitBounds(lower - growth, upper + growth, minimum, 0);
    return drag.side[index] > 0
      ? limitBounds(lower, upper + growth, minimum, 1)
      : limitBounds(lower - growth, upper, minimum, -1);
  });
}

function movedBounds(drag) {
  return AXES.map((axis, index) => [drag.baseLmin[index] + drag.offset[index], drag.baseLmax[index] + drag.offset[index]]);
}

// Keeps a block at least `minimum` wide, holding the face that is not moving.
function limitBounds(lower, upper, minimum, movingSide) {
  if (upper - lower >= minimum) return [lower, upper];
  if (movingSide > 0) return [lower, lower + minimum];
  if (movingSide < 0) return [upper - minimum, upper];
  const center = 0.5 * (lower + upper);
  return [center - 0.5 * minimum, center + 0.5 * minimum];
}

function dragsAxis(drag, index) {
  return drag.axis.includes(AXES[index].toUpperCase());
}

// Moves the dragged faces onto the nearest face of another visible block (or of
// the original when Alt-dragging a copy) within SNAP_PIXELS.
function snapBounds(drag, bounds, alt) {
  if (drag.mode === "scale" && drag.axis === "XYZ") return [];
  const center = new THREE.Vector3(...bounds.map(([lower, upper]) => 0.5 * (lower + upper)));
  const tolerance = SNAP_PIXELS * worldUnitsPerPixel(center);
  const snaps = [];
  for (let index = 0; index < 3; index += 1) {
    if (!dragsAxis(drag, index)) continue;
    const [lower, upper] = bounds[index];
    const moving = drag.mode === "translate" || alt ? [-1, 1] : [drag.side[index]];
    let best = null;
    for (const target of snapTargets(drag, index, alt)) {
      for (const side of moving) {
        const shift = target.value - (side < 0 ? lower : upper);
        if (Math.abs(shift) <= tolerance && (!best || Math.abs(shift) < Math.abs(best.shift))) best = { ...target, side, shift };
      }
    }
    if (!best) continue;
    let snapped;
    if (drag.mode === "translate") snapped = [lower + best.shift, upper + best.shift];
    else if (alt) snapped = best.side < 0 ? [lower + best.shift, upper - best.shift] : [lower - best.shift, upper + best.shift];
    else snapped = best.side < 0 ? [lower + best.shift, upper] : [lower, upper + best.shift];
    if (snapped[1] - snapped[0] < MIN_DRAG_FRACTION * (drag.baseLmax[index] - drag.baseLmin[index])) continue;
    bounds[index] = snapped;
    snaps.push({ index, value: best.value, target: best.target });
  }
  return snaps;
}

function snapTargets(drag, index, alt) {
  const targets = state.project.blocks.filter((block) => block.id !== drag.blockId && !ui.hiddenBlockIds.has(block.id));
  if (drag.mode === "translate" && alt) targets.push({ lmin: drag.baseLmin, lmax: drag.baseLmax });
  return targets.flatMap((target) => [
    { value: target.lmin[index], target },
    { value: target.lmax[index], target },
  ]);
}

function worldUnitsPerPixel(point) {
  const height = Math.max(1, ui.renderer.domElement.clientHeight);
  if (ui.camera.isOrthographicCamera) return (ui.camera.top - ui.camera.bottom) / ui.camera.zoom / height;
  const distance = ui.camera.position.distanceTo(point);
  return (2 * distance * Math.tan(THREE.MathUtils.degToRad(ui.camera.fov) / 2)) / height;
}

// Outlines each snapped plane across both the dragged block and the block snapped to.
function drawSnapGuides(block, snaps) {
  clearDisposableGroup(ui.snapGroup);
  for (const { index, value, target } of snaps) {
    const span = {
      lmin: block.lmin.map((lower, axis) => Math.min(lower, target.lmin[axis])),
      lmax: block.lmax.map((upper, axis) => Math.max(upper, target.lmax[axis])),
    };
    const line = new THREE.Line(
      new THREE.BufferGeometry().setFromPoints(planeOutlinePoints(span, index, value)),
      // Transparent so that it draws after, and over, the translucent blocks.
      new THREE.LineBasicMaterial({ color: SNAP_COLOR, transparent: true, opacity: 1, depthTest: false })
    );
    line.renderOrder = 9;
    ui.snapGroup.add(line);
  }
  els.scene.dataset.snapFaces = snaps.map(({ index, value }) => `${AXES[index]}=${value}`).join(",");
}

// Ends an Alt-drag of the Move tool: the copy keeps the new position and the
// original goes back to where the drag started.
function duplicateDraggedBlock(drag) {
  const block = state.project.blocks.find((item) => item.id === drag.blockId);
  const moved = block && block.lmin.some((value, index) => value !== drag.baseLmin[index]);
  if (!moved) return null;
  const copy = copyBlock(block, state.project.blocks, state.project.nscal);
  block.lmin = [...drag.baseLmin];
  block.lmax = [...drag.baseLmax];
  drag.baseFaces.forEach((faces, index) => {
    if (faces) ensureAxis(block, AXES[index]).faces = [...faces];
  });
  state.project.blocks.push(copy);
  state.selectedId = copy.id;
  state.selectedIds = new Set([copy.id]);
  markProjectDirty(copy.id);
  return copy;
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
