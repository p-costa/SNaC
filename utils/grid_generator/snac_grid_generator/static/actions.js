// Server-backed project actions and block editing commands.

import { postJson } from "./api.js";
import { describeAxis } from "./axis-text.js";
import { AXES, defaultBlock, nextBlockId, rescaleExplicitFaces, resetFriendBoundaries } from "./block-model.js";
import { keepModelInView } from "./camera.js";
import {
  loadProject,
  markProjectDirty,
  projectRequestSnapshot,
  projectResponseIsCurrent,
  renderAll,
  setStatus,
} from "./core.js";
import { downloadText, escapeHtml, formatNumber, positiveNumber } from "./format.js";
import { clearAutosave, commitHistory, resetHistory, scheduleHistoryCommit } from "./history.js";
import { renderCheckResults, renderDecompositionSummary, renderQualitySummary } from "./results.js";
import { emptyProject } from "./snac-templates.js";
import {
  blockLabel,
  els,
  firstSelectedBlockId,
  previewCache,
  selectedBlock,
  selectedBlockIds,
  state,
} from "./state.js";
import { renderScene } from "./viewport.js";

export async function applyAxisToAlignedBlocks() {
  const block = selectedBlock();
  if (!block) return;
  setStatus(`Applying ${state.currentAxis.toUpperCase()} grid...`);
  const requestProject = projectRequestSnapshot();
  try {
    const payload = await postJson("/api/apply-axis", {
      project: state.project,
      sourceBlockId: block.id,
      axis: state.currentAxis,
    });
    if (!projectResponseIsCurrent(requestProject, "Grid propagation")) return;
    state.project = payload.project;
    previewCache.clear();
    state.previewSequence += 1;
    state.lastCheck = null;
    state.decompositionResult = null;
    renderAll();
    scheduleHistoryCommit();
    const count = payload.changedBlockIds?.length ?? 0;
    setStatus(count ? `Applied ${state.currentAxis.toUpperCase()} grid to ${count} aligned block${count === 1 ? "" : "s"}` : "No other aligned blocks");
  } catch (error) {
    setStatus(error.message);
  }
}

// Opens the repair preview with the default source grids.
export function openGridRepair() {
  state.repairSources = {};
  return previewGridRepair();
}

export async function previewGridRepair() {
  if (!state.project.blocks.length) {
    setStatus("Add a block before repairing grids");
    return;
  }
  setStatus("Checking grid congruence...");
  const requestProject = projectRequestSnapshot();
  try {
    const payload = await postJson(
      "/api/repair",
      { project: state.project, sourceBlockId: state.selectedId, sources: state.repairSources },
      { allowInvalid: true }
    );
    if (!projectResponseIsCurrent(requestProject, "Grid repair")) return;
    if (!payload.ok) {
      state.pendingRepair = null;
      if (els.repairDialog.open) els.repairDialog.close();
      state.lastCheck = payload;
      renderCheckResults();
      setStatus(`${payload.errors.length} repair issue${payload.errors.length === 1 ? "" : "s"}`);
      return;
    }
    if (!payload.changes?.length) {
      state.pendingRepair = null;
      setStatus("Grids and interface spacing are already congruent");
      return;
    }
    state.pendingRepair = payload;
    renderRepairSummary(payload);
    if (!els.repairDialog.open) els.repairDialog.showModal();
    setStatus(`${payload.changes.length} proposed grid repair${payload.changes.length === 1 ? "" : "s"}`);
  } catch (error) {
    setStatus(error.message);
  }
}

function blockNames(blockIds) {
  return blockIds
    .map((blockId) => {
      const block = state.project.blocks.find((item) => item.id === blockId);
      return block ? blockLabel(block) : `Block ${blockId}`;
    })
    .join(", ");
}

function renderRepairSummary(payload) {
  // Each conflict lists the distinct grids; the checked one is kept.
  const choices = (payload.choices ?? [])
    .map((choice, index) => {
      const hasLocked = choice.candidates.some((candidate) => candidate.locked);
      const options = choice.candidates
        .map((candidate) => {
          const source = candidate.blockIds[0];
          const checked = candidate.blockIds.includes(choice.sourceBlockId);
          const disabled = hasLocked && !candidate.locked;
          return `<label class="repair-option${disabled ? " disabled" : ""}">
            <input type="radio" name="repair-choice-${index}" data-axis="${escapeHtml(choice.axis)}" value="${source}"${checked ? " checked" : ""}${disabled ? " disabled" : ""}>
            <span><strong>${escapeHtml(blockNames(candidate.blockIds))}</strong>${candidate.locked ? " · locked" : ""}<br>${escapeHtml(describeAxis(candidate.axis, candidate.cells))}</span>
          </label>`;
        })
        .join("");
      return `<fieldset class="repair-choice"><legend>Keep which ${choice.axis.toUpperCase()} grid?</legend>${options}</fieldset>`;
    })
    .join("");
  const changes = payload.changes
    .map((change) => {
      const block = escapeHtml(blockNames([change.blockId]));
      if (change.kind === "spacing") {
        return `<div class="repair-change"><strong>${block} · ${change.face.toUpperCase()}</strong><br>first cell ${formatNumber(change.oldSpacing)} → ${formatNumber(change.newSpacing)}, matched to ${escapeHtml(blockNames([change.sourceBlockId]))}</div>`;
      }
      return `<div class="repair-change"><strong>${block} · ${change.axis.toUpperCase()}</strong><br>${escapeHtml(describeAxis(change.oldAxis, change.oldN))} → ${escapeHtml(describeAxis(change.newAxis, change.newN))}</div>`;
    })
    .join("");
  const warnings = (payload.warnings ?? []).map((warning) => `<div class="repair-warning">${escapeHtml(warning)}</div>`).join("");
  els.repairSummary.innerHTML = `${choices}<div class="repair-changes">${changes}</div>${warnings}`;
  for (const input of els.repairSummary.querySelectorAll('input[type="radio"]')) {
    input.addEventListener("change", () => {
      state.repairSources = { ...state.repairSources, [input.dataset.axis]: [Number(input.value)] };
      previewGridRepair();
    });
  }
}

export function closeRepairDialog() {
  state.pendingRepair = null;
  state.repairSources = {};
  if (els.repairDialog.open) els.repairDialog.close();
}

export async function applyGridRepair() {
  if (!state.pendingRepair?.project) return;
  const count = state.pendingRepair.changes?.length ?? 0;
  state.project = state.pendingRepair.project;
  state.pendingRepair = null;
  state.decompositionResult = null;
  if (els.repairDialog.open) els.repairDialog.close();
  previewCache.clear();
  state.previewSequence += 1;
  state.lastCheck = null;
  renderAll();
  scheduleHistoryCommit();
  const options = { successMessage: `Applied ${count} grid repair${count === 1 ? "" : "s"}` };
  if (state.project.inferConnectivity) await updateStructure(options);
  else {
    const valid = await checkProject({ silent: true });
    setStatus(valid ? options.successMessage : "Grid repair applied; checks still need attention");
  }
}

export async function balanceDecomposition() {
  setStatus("Balancing MPI decomposition...");
  const requestProject = projectRequestSnapshot();
  try {
    const payload = await postJson("/api/decompose", { project: state.project }, { allowInvalid: true });
    if (!projectResponseIsCurrent(requestProject, "Decomposition")) return;
    state.decompositionResult = payload;
    if (!payload.ok) {
      renderDecompositionSummary();
      scheduleHistoryCommit();
      setStatus(payload.errors?.[0] ?? "No exact decomposition found");
      return;
    }
    state.project = payload.project;
    previewCache.clear();
    state.previewSequence += 1;
    state.lastCheck = null;
    renderAll();
    scheduleHistoryCommit();
    const warningText = payload.warnings?.length ? `, ${payload.warnings.length} warning${payload.warnings.length === 1 ? "" : "s"}` : "";
    setStatus(`Balanced ${payload.totalRanks} MPI ranks${warningText}`);
    if (state.project.inferConnectivity) await updateStructure({ silent: true });
    else await checkProject({ silent: true });
  } catch (error) {
    setStatus(error.message);
  }
}

export async function applyBoundaryPreset() {
  const blockIds = selectedBlockIds();
  if (!blockIds.length) return;
  setStatus("Applying boundary preset...");
  const requestProject = projectRequestSnapshot();
  try {
    const payload = await postJson("/api/bc-preset", {
      project: state.project,
      blockIds,
      preset: els.bcPreset.value,
      face: els.bcPresetFace.value,
      velocity: [els.bcPresetU.value, els.bcPresetV.value, els.bcPresetW.value].map(Number),
      pressure: Number(els.bcPresetP.value) || 0,
    });
    if (!projectResponseIsCurrent(requestProject, "Boundary preset")) return;
    state.project = payload.project;
    markProjectDirty();
    renderAll();
    setStatus(`Applied preset to ${blockIds.length} block${blockIds.length === 1 ? "" : "s"}`);
  } catch (error) {
    setStatus(error.message);
  }
}

export function addFreeBlock() {
  const source = selectedBlock() ?? state.project.blocks[0];
  const copy = source ? structuredClone(source) : defaultBlock(nextBlockId(state.project.blocks));
  copy.id = nextBlockId(state.project.blocks);
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
  resetFriendBoundaries(copy, state.project.nscal);
  state.project.blocks.push(copy);
  state.selectedId = copy.id;
  state.selectedIds = new Set([copy.id]);
  markProjectDirty(copy.id);
  renderAll();
  keepModelInView();
}

export function addAdjacentBlock() {
  const source = selectedBlock();
  if (!source) return;
  const face = els.adjacentFace.value;
  const axis = AXES.indexOf(face[0]);
  const sign = face[1] === "+" ? 1 : -1;
  const thickness = positiveNumber(els.adjacentSize.value, source.lmax[axis] - source.lmin[axis]);
  const copy = structuredClone(source);
  copy.id = nextBlockId(state.project.blocks);
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
  resetFriendBoundaries(copy, state.project.nscal);
  state.project.blocks.push(copy);
  state.selectedId = copy.id;
  state.selectedIds = new Set([copy.id]);
  markProjectDirty(copy.id);
  renderAll();
  keepModelInView();
}

export async function applyGeometryOperation(operation) {
  if (!selectedBlock()) return;
  let blockIds = selectedBlockIds();
  const payload = {
    project: state.project,
    operation,
    blockIds,
    sourceBlockId: state.selectedId,
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
    blockIds = blockIds.filter((id) => id !== state.selectedId);
    if (blockIds.length !== 1) {
      setStatus("Select one target block and make the source block primary");
      return;
    }
    payload.blockIds = blockIds;
    payload.face = els.snapFace.value;
  }
  setStatus(`${operation[0].toUpperCase()}${operation.slice(1)} blocks...`);
  const requestProject = projectRequestSnapshot();
  try {
    const result = await postJson("/api/geometry", payload);
    if (!projectResponseIsCurrent(requestProject, "Geometry operation")) return;
    state.project = result.project;
    const changed = result.changedBlockIds ?? [];
    if (["duplicate", "mirror"].includes(operation) && changed.length) {
      state.selectedIds = new Set(changed);
      state.selectedId = changed[0];
    }
    previewCache.clear();
    state.previewSequence += 1;
    markProjectDirty();
    renderAll();
    keepModelInView();
    setStatus(`${operation[0].toUpperCase()}${operation.slice(1)} updated ${changed.length} block${changed.length === 1 ? "" : "s"}`);
  } catch (error) {
    setStatus(error.message);
  }
}

export function removeSelectedBlock() {
  if (!selectedBlock()) return;
  const removedId = state.selectedId;
  state.project.blocks = state.project.blocks.filter((block) => block.id !== state.selectedId);
  previewCache.delete(removedId);
  state.selectedIds.delete(removedId);
  state.selectedId = firstSelectedBlockId() ?? state.project.blocks[0]?.id ?? null;
  if (state.selectedId != null) state.selectedIds.add(state.selectedId);
  markProjectDirty();
  renderAll();
}

// Starts a new empty project with default settings; Undo brings back the old one.
export function resetProject() {
  commitHistory();
  loadProject(emptyProject(), { fit: false });
  commitHistory();
  setStatus("Started a new empty project; Undo brings back the previous one");
}

export async function checkProject(options = {}) {
  if (!options.silent) setStatus("Checking grid...");
  const requestProject = projectRequestSnapshot();
  try {
    const payload = await postJson("/api/check", { project: state.project }, { allowInvalid: true });
    if (!projectResponseIsCurrent(requestProject, "Grid check")) return false;
    state.lastCheck = payload;
    renderCheckResults();
    renderQualitySummary();
    renderScene();
    if (!payload.ok) {
      setStatus(`${payload.errors.length} check issue${payload.errors.length === 1 ? "" : "s"}`);
      return false;
    }
    const warningText = payload.warnings?.length ? `, ${payload.warnings.length} warnings` : "";
    const staleText = payload.connectivityStale ? "; Update applies the inferred friend boundaries" : "";
    if (!options.silent) setStatus(`Grid checks passed${warningText}${staleText}`);
    return true;
  } catch (error) {
    state.lastCheck = { ok: false, errors: [error.message], warnings: [] };
    renderCheckResults();
    renderQualitySummary();
    renderScene();
    setStatus(error.message);
    return false;
  }
}

export async function updateStructure(options = {}) {
  if (!options.silent) setStatus("Updating structure...");
  const requestProject = projectRequestSnapshot();
  try {
    const payload = await postJson("/api/update", { project: state.project, sourceBlockId: state.selectedId }, { allowInvalid: true });
    if (!projectResponseIsCurrent(requestProject, "Structure update")) return false;
    if (payload.project) {
      const oldSelectedId = state.selectedId;
      state.project = payload.project;
      state.selectedId = state.project.blocks.find((block) => block.id === oldSelectedId)?.id ?? state.project.blocks[0]?.id ?? null;
      previewCache.clear();
      state.previewSequence += 1;
    }
    state.lastCheck = payload;
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

export async function exportCase() {
  const valid = state.project.inferConnectivity ? await updateStructure({ silent: true }) : await checkProject({ silent: true });
  if (!valid) {
    setStatus("Fix grid checks before writing");
    return;
  }
  setStatus("Writing files...");
  const requestProject = projectRequestSnapshot();
  try {
    const payload = await postJson("/api/export", { outputDir: els.outputDir.value, project: state.project });
    if (!projectResponseIsCurrent(requestProject, "Export")) {
      setStatus(`Wrote an earlier project state to ${payload.outputDir}; export again`);
      return;
    }
    commitHistory();
    clearAutosave();
    const warningText = payload.warnings?.length ? `, ${payload.warnings.length} warnings` : "";
    setStatus(`Wrote ${payload.files.length} files to ${payload.outputDir}${warningText}`);
  } catch (error) {
    setStatus(error.message);
  }
}

export async function importCase() {
  const caseDir = els.outputDir.value.trim();
  if (!caseDir) {
    setStatus("Set the case folder in Output");
    return;
  }
  setStatus("Importing SNaC case...");
  try {
    const payload = await postJson("/api/import", { caseDir });
    loadProject(payload.project);
    clearAutosave();
    resetHistory(false);
    const warningText = payload.warnings?.length
      ? `, ${payload.warnings.length} warning${payload.warnings.length === 1 ? "" : "s"}`
      : "";
    setStatus(`Imported ${state.project.blocks.length} blocks from ${payload.source}${warningText}`);
    await checkProject({ silent: true });
  } catch (error) {
    setStatus(error.message);
  }
}

export function downloadProjectJson() {
  const blob = new Blob([JSON.stringify(state.project, null, 2) + "\n"], { type: "application/json" });
  const link = document.createElement("a");
  link.href = URL.createObjectURL(blob);
  link.download = `${state.project.name || "snac-grid"}.json`;
  link.click();
  URL.revokeObjectURL(link.href);
  commitHistory();
  clearAutosave();
}

export async function downloadCaseReport(format) {
  setStatus(`Building ${format === "json" ? "JSON" : "Markdown"} report...`);
  const requestProject = projectRequestSnapshot();
  try {
    const payload = await postJson("/api/report", { project: state.project });
    if (!projectResponseIsCurrent(requestProject, "Report")) return;
    const extension = format === "json" ? "json" : "md";
    const type = format === "json" ? "application/json" : "text/markdown";
    const content = format === "json"
      ? `${JSON.stringify(payload.report, null, 2)}\n`
      : payload.markdown;
    const base = (state.project.name || "snac-grid").replace(/[^A-Za-z0-9._-]+/g, "-");
    downloadText(content, `${base}-quality.${extension}`, type);
    setStatus(`Downloaded ${extension.toUpperCase()} quality report`);
  } catch (error) {
    setStatus(error.message);
  }
}

export async function openProjectJson(event) {
  const file = event.target.files[0];
  if (!file) return;
  event.target.value = "";
  try {
    const loaded = JSON.parse(await file.text());
    const payload = await postJson("/api/migrate", { project: loaded });
    loadProject(payload.project);
    clearAutosave();
    resetHistory(false);
    setStatus(`Loaded ${file.name}`);
  } catch (error) {
    setStatus(error.message);
  }
}
