// Shared editor state, DOM handles, and selectors over the current selection.

import { loadLibrary } from "./library.js";

export const els = {};
export const previewCache = new Map();

// Editor and project state shared across modules.
export const state = {
  project: null,
  selectedId: 1,
  selectedIds: new Set([1]),
  currentAxis: "x",
  lastCheck: null,
  activeDiagnostic: null,
  pendingRepair: null,
  decompositionResult: null,
  previewTimer: null,
  previewSequence: 0,
  draggedBlockId: null,
  historyTimer: null,
  pendingRecovery: null,
  userLibrary: loadLibrary(),
  currentLibraryTab: "axisPresets",
  selectedLibraryItem: "builtin:builtin-uniform-32",
  repairSources: {},
};

// three.js scene objects and viewport display settings.
export const ui = {
  scene: null,
  camera: null,
  perspectiveCamera: null,
  orthographicCamera: null,
  orthographicViewSize: null,
  renderer: null,
  controls: null,
  raycaster: null,
  pointer: null,
  blockGroup: null,
  boundaryFaceGroup: null,
  gridHelper: null,
  axesGroup: null,
  partitionGroup: null,
  gridLineGroup: null,
  spacingJumpGroup: null,
  diagnosticGroup: null,
  transformControls: null,
  selectedBlockMesh: null,
  transformDrag: null,
  isTransformDragging: false,
  scenePointerDown: null,
  groundGridVisible: true,
  boundaryColorsVisible: true,
  partitionPlanesVisible: true,
  gridLinesMode: "selected",
  spacingJumpsVisible: true,
  blockTool: "select",
  clippingEnabled: false,
  clippingAxisIndex: 0,
  clippingFraction: 0.5,
  clippingDirection: 1,
  hiddenBlockIds: new Set(),
  issuesCollapsed: false,
};

export const COLORS = [0x43c989, 0x76a9ff, 0xd28a3e, 0xc789e8, 0xe46868, 0x82c6c2, 0xd6bf5b, 0x9ba7ff];

export function selectedBlock() {
  return state.project?.blocks?.find((block) => block.id === state.selectedId);
}

export function selectedBlockIds() {
  const ids = state.project.blocks.filter((block) => state.selectedIds.has(block.id)).map((block) => block.id);
  return ids.length ? ids : state.selectedId == null ? [] : [state.selectedId];
}

export function firstSelectedBlockId() {
  return state.project?.blocks?.find((block) => state.selectedIds.has(block.id))?.id ?? null;
}

export function visibleBlocks() {
  return state.project.blocks.filter((block) => !ui.hiddenBlockIds.has(block.id));
}

export function selectedVisibleBlocks() {
  const ids = new Set(selectedBlockIds());
  return state.project.blocks.filter((block) => ids.has(block.id) && !ui.hiddenBlockIds.has(block.id));
}

export function blockLabel(block) {
  return block.name ? `${block.id}: ${block.name}` : `Block ${block.id}`;
}
