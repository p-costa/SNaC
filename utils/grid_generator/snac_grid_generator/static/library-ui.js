// The grid library dialog for axis presets and case templates.

import { updateStructure } from "./actions.js";
import { postJson } from "./api.js";
import { axisKindLabel } from "./axis-text.js";
import { AXES, canonicalizeAxis, ensureAxis } from "./block-model.js";
import { loadProject, markProjectDirty, renderAll, setStatus } from "./core.js";
import { downloadText, escapeHtml, formatInteger } from "./format.js";
import { commitHistory } from "./history.js";
import {
  LIBRARY_SCHEMA_VERSION,
  applyAxisPreset,
  captureAxisPreset,
  captureProjectTemplate,
  mergeLibraries,
  normalizeLibrary,
  removeLibraryItem,
  renameLibraryItem,
  replaceLibraryItem,
  saveLibrary,
} from "./library.js";
import { schedulePreview } from "./previews.js";
import { builtInAxisPresets, builtInProjectTemplates } from "./snac-templates.js";
import { blockLabel, els, selectedBlock, state } from "./state.js";

export function showLibraryDialog() {
  renderLibraryDialog();
  els.libraryDialog.showModal();
}

export function renderLibraryDialog() {
  for (const button of els.libraryTabs.querySelectorAll("button")) {
    button.classList.toggle("active", button.dataset.libraryTab === state.currentLibraryTab);
  }
  const entries = libraryEntries();
  if (!entries.some((entry) => entry.key === state.selectedLibraryItem)) state.selectedLibraryItem = entries[0]?.key ?? null;
  els.libraryItems.innerHTML = entries
    .map((entry) => `<option value="${escapeHtml(entry.key)}">${escapeHtml(entry.item.name)}${entry.builtIn ? " (built-in)" : ""}</option>`)
    .join("");
  els.libraryItems.value = state.selectedLibraryItem ?? "";
  renderLibrarySelection();
  if (window.lucide) window.lucide.createIcons();
}

export function renderLibrarySelection() {
  const entry = selectedLibraryEntry();
  const axisMode = state.currentLibraryTab === "axisPresets";
  const block = selectedBlock();
  const axisName = state.currentAxis.toUpperCase();
  // Save stores the current axis grid or project, so name it after that.
  const suggestedName = axisMode
    ? block ? `${block.name || `Block ${block.id}`} ${axisName} grid` : ""
    : state.project.name || "Case template";
  const source = axisMode
    ? block ? `the ${axisName} grid of ${blockLabel(block)}` : null
    : "the current project";
  if (!entry) {
    els.libraryName.value = suggestedName;
    els.librarySummary.textContent = "Empty library";
  } else {
    els.libraryName.value = entry.builtIn ? suggestedName : entry.item.name;
    els.librarySummary.innerHTML = libraryItemSummary(entry);
  }
  els.applyLibraryItem.disabled = !entry || (axisMode && !block);
  els.saveLibraryItem.disabled = axisMode && !block;
  els.saveLibraryItem.title = !source
    ? "Select a block to save one of its axis grids"
    : entry && !entry.builtIn
      ? `Update ${entry.item.name} with ${source}`
      : `Save ${source} as a new library item`;
  els.renameLibraryItem.disabled = !entry || entry.builtIn;
  els.deleteLibraryItem.disabled = !entry || entry.builtIn;
}

export function libraryEntries() {
  const builtIns = state.currentLibraryTab === "axisPresets" ? builtInAxisPresets() : builtInProjectTemplates();
  return [
    ...builtIns.map((item) => ({ key: `builtin:${item.id}`, item, builtIn: true })),
    ...state.userLibrary[state.currentLibraryTab].map((item) => ({ key: `user:${item.id}`, item, builtIn: false })),
  ];
}

export function selectedLibraryEntry() {
  return libraryEntries().find((entry) => entry.key === state.selectedLibraryItem) ?? null;
}

export function libraryItemSummary(entry) {
  const source = entry.builtIn ? "Built-in" : "User";
  if (state.currentLibraryTab === "axisPresets") {
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

export async function applySelectedLibraryItem() {
  const entry = selectedLibraryEntry();
  if (!entry) return;
  if (state.currentLibraryTab === "axisPresets") {
    const block = selectedBlock();
    if (!block) return;
    try {
      const axisIndex = AXES.indexOf(state.currentAxis);
      const applied = applyAxisPreset(entry.item, block.lmin[axisIndex], block.lmax[axisIndex]);
      block.axes[state.currentAxis] = canonicalizeAxis(applied.axis);
      block.ng[axisIndex] = applied.ng;
      markProjectDirty(block.id);
      renderAll();
      schedulePreview(block.id, 0);
      setStatus(`Applied ${entry.item.name} to ${state.currentAxis.toUpperCase()}`);
    } catch (error) {
      setStatus(error.message);
    }
    return;
  }
  setStatus(`Applying ${entry.item.name}...`);
  try {
    commitHistory();
    const payload = await postJson("/api/migrate", { project: entry.item.project });
    loadProject(payload.project);
    commitHistory();
    els.libraryDialog.close();
    const structured = state.project.inferConnectivity ? await updateStructure({ silent: true }) : true;
    if (structured) setStatus(`Applied ${entry.item.name}`);
  } catch (error) {
    setStatus(error.message);
  }
}

export function saveCurrentLibraryItem() {
  const selected = selectedLibraryEntry();
  const existingId = selected && !selected.builtIn ? selected.item.id : undefined;
  try {
    assertLibraryNameAvailable(els.libraryName.value, existingId);
    let item;
    if (state.currentLibraryTab === "axisPresets") {
      const block = selectedBlock();
      if (!block) throw new Error("Select a block before saving an axis preset");
      const axisIndex = AXES.indexOf(state.currentAxis);
      item = captureAxisPreset(
        els.libraryName.value,
        ensureAxis(block, state.currentAxis),
        block.ng[axisIndex],
        block.lmin[axisIndex],
        block.lmax[axisIndex],
        existingId
      );
    } else {
      item = captureProjectTemplate(els.libraryName.value, state.project, existingId);
    }
    state.userLibrary = saveLibrary(replaceLibraryItem(state.userLibrary, state.currentLibraryTab, item));
    state.selectedLibraryItem = `user:${item.id}`;
    renderLibraryDialog();
    setStatus(`${existingId ? "Updated" : "Saved"} ${item.name}`);
  } catch (error) {
    setStatus(error.message);
  }
}

export function renameSelectedLibraryItem() {
  const entry = selectedLibraryEntry();
  if (!entry || entry.builtIn) return;
  try {
    assertLibraryNameAvailable(els.libraryName.value, entry.item.id);
    state.userLibrary = saveLibrary(renameLibraryItem(state.userLibrary, state.currentLibraryTab, entry.item.id, els.libraryName.value));
    renderLibraryDialog();
    setStatus(`Renamed library item to ${els.libraryName.value}`);
  } catch (error) {
    setStatus(error.message);
  }
}

export function deleteSelectedLibraryItem() {
  const entry = selectedLibraryEntry();
  if (!entry || entry.builtIn) return;
  const name = entry.item.name;
  try {
    state.userLibrary = saveLibrary(removeLibraryItem(state.userLibrary, state.currentLibraryTab, entry.item.id));
    state.selectedLibraryItem = null;
    renderLibraryDialog();
    setStatus(`Deleted ${name}`);
  } catch (error) {
    setStatus(error.message);
  }
}

export async function importLibraryJson(event) {
  const file = event.target.files[0];
  if (!file) return;
  event.target.value = "";
  try {
    const incoming = normalizeLibrary(JSON.parse(await file.text()));
    const merged = mergeLibraries(state.userLibrary, incoming);
    state.userLibrary = saveLibrary(merged.library);
    state.selectedLibraryItem = null;
    renderLibraryDialog();
    const count = merged.counts.axisPresets + merged.counts.projectTemplates;
    setStatus(`Imported ${count} library item${count === 1 ? "" : "s"}`);
  } catch (error) {
    setStatus(error.message);
  }
}

export function exportLibraryJson() {
  const content = `${JSON.stringify({ ...state.userLibrary, schemaVersion: LIBRARY_SCHEMA_VERSION }, null, 2)}\n`;
  downloadText(content, "snac-grid-library.json", "application/json");
  setStatus("Downloaded grid library");
}

export function assertLibraryNameAvailable(name, exceptId) {
  const target = String(name ?? "").trim().toLocaleLowerCase();
  const duplicate = libraryEntries().some((entry) => entry.item.id !== exceptId && entry.item.name.toLocaleLowerCase() === target);
  if (duplicate) throw new Error("Library item names must be unique");
}
