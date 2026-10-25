// The block list: selection, drag-and-drop ordering, and renumbering.

import { renumberBlocksByOrder } from "./block-model.js";
import { markProjectDirty, renderAll, setStatus } from "./core.js";
import { escapeHtml } from "./format.js";
import { COLORS, blockLabel, els, firstSelectedBlockId, previewCache, state, ui } from "./state.js";

export function renderBlockList() {
  els.blockList.innerHTML = "";
  if (!state.project.blocks.length) {
    els.blockList.innerHTML = `<div class="empty-state">No blocks</div>`;
    return;
  }
  for (const block of state.project.blocks) {
    const item = document.createElement("div");
    const selected = state.selectedIds.has(block.id);
    item.className = [
      "block-item",
      selected ? "selected" : "",
      block.id === state.selectedId ? "primary" : "",
      ui.hiddenBlockIds.has(block.id) ? "hidden-block" : "",
    ].filter(Boolean).join(" ");
    item.dataset.blockId = String(block.id);
    item.innerHTML = `
      <button type="button" class="block-drag-handle" draggable="true" aria-label="Reorder ${escapeHtml(blockLabel(block))}" title="Drag to reorder"><i data-lucide="grip-vertical"></i></button>
      <input type="checkbox" aria-label="Select ${escapeHtml(blockLabel(block))}" ${selected ? "checked" : ""} />
      <span class="swatch" style="background:#${COLORS[(block.id - 1) % COLORS.length].toString(16).padStart(6, "0")}"></span>
      <button type="button" class="block-main"><strong>${escapeHtml(blockLabel(block))}</strong><span>${block.lmin.join(", ")} -> ${block.lmax.join(", ")}</span></button>
    `;
    item.querySelector("input").addEventListener("change", (event) => {
      setBlockSelected(block.id, event.target.checked, false);
    });
    item.querySelector(".block-main").addEventListener("click", (event) => {
      selectBlock(block.id, event.shiftKey);
    });
    const dragHandle = item.querySelector(".block-drag-handle");
    dragHandle.addEventListener("dragstart", (event) => startBlockDrag(event, block.id));
    dragHandle.addEventListener("dragend", clearBlockDrag);
    dragHandle.addEventListener("keydown", (event) => handleBlockOrderKey(event, block.id));
    item.addEventListener("dragover", (event) => showBlockDropPosition(event, block.id));
    item.addEventListener("dragleave", (event) => clearBlockDropPosition(event.currentTarget));
    item.addEventListener("drop", (event) => dropBlock(event, block.id));
    els.blockList.appendChild(item);
  }
}

export function selectBlock(blockId, toggle = false) {
  if (!state.project.blocks.some((block) => block.id === blockId)) return;
  if (toggle) {
    setBlockSelected(blockId, !state.selectedIds.has(blockId), true);
    return;
  }
  state.selectedId = blockId;
  state.selectedIds = new Set([blockId]);
  state.activeDiagnostic = null;
  renderAll();
}

export function setBlockSelected(blockId, selected, makePrimary) {
  if (!state.project.blocks.some((block) => block.id === blockId)) return;
  if (selected) {
    state.selectedIds.add(blockId);
    if (makePrimary || state.selectedId == null) state.selectedId = blockId;
  } else {
    state.selectedIds.delete(blockId);
    if (state.selectedId === blockId) state.selectedId = firstSelectedBlockId();
  }
  state.activeDiagnostic = null;
  renderAll();
}

export function startBlockDrag(event, blockId) {
  state.draggedBlockId = blockId;
  event.dataTransfer.effectAllowed = "move";
  event.dataTransfer.setData("text/plain", String(blockId));
  els.blockList.querySelector(`[data-block-id="${blockId}"]`)?.classList.add("dragging");
}

export function showBlockDropPosition(event, targetId) {
  if (state.draggedBlockId == null || state.draggedBlockId === targetId) return;
  event.preventDefault();
  event.dataTransfer.dropEffect = "move";
  clearBlockDropPositions();
  const item = event.currentTarget;
  const after = event.clientY >= item.getBoundingClientRect().top + item.offsetHeight / 2;
  item.classList.add(after ? "drop-after" : "drop-before");
}

export function clearBlockDropPosition(item) {
  item.classList.remove("drop-before", "drop-after");
}

export function clearBlockDropPositions() {
  for (const item of els.blockList.querySelectorAll(".block-item")) clearBlockDropPosition(item);
}

export function dropBlock(event, targetId) {
  if (state.draggedBlockId == null || state.draggedBlockId === targetId) return;
  event.preventDefault();
  const after = event.currentTarget.classList.contains("drop-after");
  const blockId = state.draggedBlockId;
  clearBlockDrag();
  reorderBlock(blockId, targetId, after);
}

export function clearBlockDrag() {
  state.draggedBlockId = null;
  clearBlockDropPositions();
  for (const item of els.blockList.querySelectorAll(".dragging")) item.classList.remove("dragging");
}

export function handleBlockOrderKey(event, blockId) {
  if (!["ArrowUp", "ArrowDown"].includes(event.key)) return;
  event.preventDefault();
  const index = state.project.blocks.findIndex((block) => block.id === blockId);
  const targetIndex = index + (event.key === "ArrowUp" ? -1 : 1);
  if (index < 0 || targetIndex < 0 || targetIndex >= state.project.blocks.length) return;
  const [block] = state.project.blocks.splice(index, 1);
  state.project.blocks.splice(targetIndex, 0, block);
  finishBlockReorder(block, targetIndex);
}

export function reorderBlock(blockId, targetId, after) {
  const sourceIndex = state.project.blocks.findIndex((block) => block.id === blockId);
  if (sourceIndex < 0) return;
  const [block] = state.project.blocks.splice(sourceIndex, 1);
  const targetIndex = state.project.blocks.findIndex((item) => item.id === targetId);
  if (targetIndex < 0) {
    state.project.blocks.splice(sourceIndex, 0, block);
    return;
  }
  const insertionIndex = targetIndex + (after ? 1 : 0);
  state.project.blocks.splice(insertionIndex, 0, block);
  finishBlockReorder(block, insertionIndex);
}

export function finishBlockReorder(block, index) {
  markProjectDirty();
  renderAll();
  setStatus(`Moved ${blockLabel(block)} to position ${index + 1}`);
}

export function renumberBlocksFromOrder() {
  const changed = state.project.blocks.some((block, index) =>
    block.id !== index + 1 || block.name !== `block-${index + 1}`
  );
  if (!changed) {
    setStatus("Block IDs and names already follow the list order");
    return;
  }
  const idMap = renumberBlocksByOrder(state.project.blocks);
  state.selectedId = idMap.get(state.selectedId) ?? null;
  state.selectedIds = new Set([...state.selectedIds].map((id) => idMap.get(id)).filter((id) => id != null));
  ui.hiddenBlockIds = new Set([...ui.hiddenBlockIds].map((id) => idMap.get(id)).filter((id) => id != null));
  state.activeDiagnostic = null;
  previewCache.clear();
  state.previewSequence += 1;
  markProjectDirty();
  renderAll();
  setStatus("Renumbered block IDs and names from the current list order");
}
