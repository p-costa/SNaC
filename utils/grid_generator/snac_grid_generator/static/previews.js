// Debounced server previews of block grids.

import { postJson } from "./api.js";
import { AXES } from "./block-model.js";
import { markProjectDirty, setStatus } from "./core.js";
import { scheduleHistoryCommit } from "./history.js";
import { renderSpacingPreview, updateSegmentMetrics } from "./inspector.js";
import { els, previewCache, selectedVisibleBlocks, state, ui, visibleBlocks } from "./state.js";
import { renderScene } from "./viewport.js";

export function schedulePreview(blockId = state.selectedId, delay = 120) {
  const blockIds = [
    blockId,
    ...gridPreviewBlocks().filter((block) => !previewCache.has(block.id)).map((block) => block.id),
  ];
  schedulePreviews(blockIds, delay);
}

export function scheduleGridPreviews(delay = 120) {
  schedulePreviews(
    gridPreviewBlocks().filter((block) => !previewCache.has(block.id)).map((block) => block.id),
    delay,
  );
}

export function gridPreviewBlocks() {
  return ui.gridLinesMode === "all" ? visibleBlocks() : selectedVisibleBlocks();
}

// Cancels the requests of the previous preview batch when a new one starts.
let previewRequests = null;

export function schedulePreviews(blockIds, delay) {
  const uniqueIds = [...new Set(blockIds)].filter((blockId) =>
    blockId != null && state.project.blocks.some((block) => block.id === blockId)
  );
  if (!uniqueIds.length) return;
  window.clearTimeout(state.previewTimer);
  previewRequests?.abort();
  previewRequests = new AbortController();
  const { signal } = previewRequests;
  const sequence = ++state.previewSequence;
  state.previewTimer = window.setTimeout(() => requestPreviews(uniqueIds, sequence, signal), delay);
}

// Previews are requested in batches. The first batch holds only the block being
// edited, so its preview arrives quickly; later batches fill in the rest.
const PREVIEW_BATCH_SIZE = 24;

export async function requestPreviews(blockIds, sequence, signal) {
  const batches = [blockIds.slice(0, 1)];
  for (let start = 1; start < blockIds.length; start += PREVIEW_BATCH_SIZE) {
    batches.push(blockIds.slice(start, start + PREVIEW_BATCH_SIZE));
  }
  for (const batch of batches) {
    if (sequence !== state.previewSequence || signal?.aborted) return;
    await requestPreviewBatch(batch, sequence, signal);
  }
}

async function requestPreviewBatch(blockIds, sequence, signal) {
  let previews;
  try {
    const payload = await postJson("/api/preview", { project: state.project, blockIds }, { signal });
    previews = payload.previews ?? [];
  } catch (error) {
    if (sequence !== state.previewSequence || signal?.aborted) return;
    previews = blockIds.map((blockId) => ({ blockId, error: error.message }));
  }
  if (sequence !== state.previewSequence) return;
  let cellCountChanged = false;
  let errorMessage = null;
  for (const preview of previews) {
    if (preview.error) {
      previewCache.delete(preview.blockId);
      errorMessage ??= preview.error;
      continue;
    }
    const block = state.project.blocks.find((item) => item.id === preview.blockId);
    if (!block) continue;
    cellCountChanged ||= block.ng.some((value, index) => value !== preview.ng[index]);
    block.ng = preview.ng;
    previewCache.set(block.id, preview.axes);
    if (block.id === state.selectedId) showSelectedPreview(block, preview);
  }
  if (cellCountChanged) scheduleHistoryCommit();
  renderSpacingPreview();
  if (!ui.isTransformDragging) renderScene(); // the drag end re-renders everything
  if (errorMessage) setStatus(errorMessage);
}

// Server-derived cell counts and segment metrics of the selected block.
function showSelectedPreview(block, preview) {
  for (const input of els.resolutionGrid.querySelectorAll('input[data-field="ng"]')) {
    input.value = block.ng[Number(input.dataset.index)];
  }
  const monotoneCells = els.axisEditor.querySelector('input[data-control="n"]');
  if (monotoneCells) monotoneCells.value = block.ng[AXES.indexOf(state.currentAxis)];
  updateSegmentMetrics(preview.axes[state.currentAxis]);
}

export function gridDefinitionChanged(block) {
  if (!block) return;
  markProjectDirty(block.id);
  schedulePreview(block.id);
  renderSpacingPreview();
  renderScene();
}
