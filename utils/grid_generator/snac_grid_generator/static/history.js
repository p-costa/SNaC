// Undo/redo snapshots, autosave, and session recovery.

import { postJson } from "./api.js";
import { loadProject, renderAll, setStatus } from "./core.js";
import {
  SnapshotHistory,
  clearStoredProject,
  listStoredProjects,
  loadStoredProject,
  saveStoredProject,
} from "./project-state.js";
import { els, previewCache, state, ui } from "./state.js";

// Each browser tab autosaves under its own key, so two open tabs never replace
// each other's recovery copy. The tab id lives in sessionStorage, which keeps it
// across reloads of the same tab.
const AUTOSAVE_PREFIX = "snac-grid-generator-autosave-v2:";
const LEGACY_AUTOSAVE_KEY = "snac-grid-generator-autosave-v1";
const MAX_AUTOSAVES = 5;
const AUTOSAVE_DELAY_MS = 400;
export const AUTOSAVE_KEY = AUTOSAVE_PREFIX + tabId();

let autosaveTimer = null;
let autosaveWarned = false;

function tabId() {
  const fallback = () => `${Date.now().toString(36)}${Math.random().toString(36).slice(2)}`;
  try {
    let id = sessionStorage.getItem("snac-grid-generator-tab");
    if (!id) {
      id = globalThis.crypto?.randomUUID?.() ?? fallback();
      sessionStorage.setItem("snac-grid-generator-tab", id);
    }
    return id;
  } catch (_error) {
    return fallback();
  }
}

export const projectHistory = new SnapshotHistory(100);

export function projectSnapshot() {
  return {
    project: JSON.stringify(state.project),
    decompositionResult: state.decompositionResult ? JSON.stringify(state.decompositionResult) : null,
    selectedId: state.selectedId,
    selectedIds: [...state.selectedIds],
    currentAxis: state.currentAxis,
  };
}

export function resetHistory(save = true) {
  window.clearTimeout(state.historyTimer);
  state.historyTimer = null;
  projectHistory.reset(projectSnapshot());
  updateHistoryButtons();
  if (save) saveAutosave();
}

export function scheduleHistoryCommit(delay = 250) {
  if (ui.isTransformDragging) return;
  saveAutosave();
  window.clearTimeout(state.historyTimer);
  state.historyTimer = window.setTimeout(commitHistory, delay);
}

export function commitHistory() {
  window.clearTimeout(state.historyTimer);
  state.historyTimer = null;
  const snapshot = projectSnapshot();
  saveAutosave();
  projectHistory.commit(snapshot);
  updateHistoryButtons();
}

export function undoProject() {
  commitHistory();
  const snapshot = projectHistory.undo();
  if (snapshot) restoreHistory(snapshot);
}

export function redoProject() {
  const snapshot = projectHistory.redo();
  if (snapshot) restoreHistory(snapshot);
}

export function restoreHistory(snapshot) {
  window.clearTimeout(state.historyTimer);
  state.historyTimer = null;
  state.project = JSON.parse(snapshot.project);
  state.selectedId = state.project.blocks.some((block) => block.id === snapshot.selectedId) ? snapshot.selectedId : state.project.blocks[0]?.id ?? null;
  state.selectedIds = new Set((snapshot.selectedIds ?? []).filter((id) => state.project.blocks.some((block) => block.id === id)));
  state.currentAxis = snapshot.currentAxis;
  state.lastCheck = null;
  state.pendingRepair = null;
  state.decompositionResult = snapshot.decompositionResult ? JSON.parse(snapshot.decompositionResult) : null;
  previewCache.clear();
  state.previewSequence += 1;
  renderAll();
  saveAutosave();
  setStatus("");
}

// Autosave writes are batched: typing does not serialize the project on every
// keystroke. Pending writes are flushed when the page is hidden or closed.
export function saveAutosave() {
  if (state.pendingRecovery) return; // keep the recovery copy until the user decides
  window.clearTimeout(autosaveTimer);
  autosaveTimer = window.setTimeout(writeAutosave, AUTOSAVE_DELAY_MS);
}

export function flushAutosave() {
  if (autosaveTimer == null) return;
  window.clearTimeout(autosaveTimer);
  writeAutosave();
}

function writeAutosave() {
  autosaveTimer = null;
  if (!state.project || state.pendingRecovery) return;
  if (saveStoredProject(AUTOSAVE_KEY, state.project)) {
    autosaveWarned = false;
    for (const stale of listStoredProjects(AUTOSAVE_PREFIX).slice(MAX_AUTOSAVES)) clearStoredProject(stale.key);
  } else if (!autosaveWarned) {
    autosaveWarned = true;
    setStatus("Autosave failed: browser storage is full or unavailable; download the project JSON to keep it");
  }
}

export function clearAutosave() {
  window.clearTimeout(autosaveTimer);
  autosaveTimer = null;
  clearStoredProject(AUTOSAVE_KEY);
}

// While a tab is open it holds a lock named after its autosave key, so other
// tabs can tell open work from work left by a closed or crashed tab.
export function claimAutosave() {
  navigator.locks?.request(AUTOSAVE_KEY, () => new Promise(() => {})).catch(() => {});
}

async function openTabAutosaveKeys() {
  try {
    const { held } = await navigator.locks.query();
    return new Set(held.map((lock) => lock.name));
  } catch (_error) {
    return new Set();
  }
}

// This tab's autosave after a reload; otherwise the newest autosave left by a
// tab that is no longer open, including the single key older versions used.
export async function loadAutosave() {
  const own = loadStoredProject(AUTOSAVE_KEY);
  if (own) return own;
  const open = await openTabAutosaveKeys();
  return listStoredProjects(AUTOSAVE_PREFIX, [LEGACY_AUTOSAVE_KEY]).find((item) => !open.has(item.key)) ?? null;
}

export function showRecoveryDialog(recovery) {
  state.pendingRecovery = recovery;
  const savedAt = new Date(recovery.savedAt);
  const time = Number.isNaN(savedAt.getTime()) ? "an earlier session" : savedAt.toLocaleString();
  els.recoverySummary.textContent = `${recovery.project.name || "Untitled project"}, saved ${time}`;
  els.recoveryDialog.showModal();
}

export async function restoreRecovery() {
  if (!state.pendingRecovery) return;
  try {
    const recovery = state.pendingRecovery;
    const payload = await postJson("/api/migrate", { project: recovery.project });
    // The recovered project now autosaves under this tab's key.
    if (recovery.key !== AUTOSAVE_KEY) clearStoredProject(recovery.key);
    state.pendingRecovery = null;
    els.recoveryDialog.close();
    loadProject(payload.project);
    resetHistory();
    setStatus("Recovered autosaved project");
  } catch (error) {
    setStatus(error.message);
  }
}

export function discardRecovery() {
  if (state.pendingRecovery?.key) clearStoredProject(state.pendingRecovery.key);
  state.pendingRecovery = null;
  clearAutosave();
  els.recoveryDialog.close();
  resetHistory(false);
  setStatus("Started with the sample project");
}

export function updateHistoryButtons() {
  if (!els.undoProject || !els.redoProject) return;
  els.undoProject.disabled = !projectHistory.canUndo;
  els.redoProject.disabled = !projectHistory.canRedo;
}
