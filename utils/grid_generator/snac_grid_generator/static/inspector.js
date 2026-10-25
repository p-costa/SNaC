// The block inspector: geometry, axis grading editors, boundary rows, and spacing preview.

import { applyAxisToAlignedBlocks, updateStructure } from "./actions.js";
import { MONOTONE_CONTROL_LABELS } from "./axis-text.js";
import { renderBlockList } from "./block-list.js";
import {
  AXES,
  FACE_ORDER,
  ensureAxis,
  ensureBoundaryArrays,
  isCellCountDerived,
  propagateMpiPartition,
  rescaleExplicitFaces,
} from "./block-model.js";
import { keepModelInView } from "./camera.js";
import { markProjectDirty } from "./core.js";
import {
  bindNumberField,
  escapeHtml,
  formatInputNumber,
  formatNumber,
  formatPercent,
  positiveNumber,
} from "./format.js";
import { gridDefinitionChanged, schedulePreview } from "./previews.js";
import {
  SNAC_GEOMETRIC_UNIFORM_GR as GEOMETRIC_UNIFORM_GR,
  SNAC_GRID_FUNCTIONS as GRID_FUNCTIONS,
  SNAC_INITIAL_VELOCITY_FIELDS as INIT_FIELDS,
} from "./snac.js";
import { blockLabel, els, previewCache, selectedBlock, state } from "./state.js";
import { renderScene } from "./viewport.js";

export const PROFILE_OPTIONS = [
  ["geometric", "Geometric"],
  ["tanh", "tanh"],
  ["erf", "erf"],
];

export const MONOTONE_CONTROL_PAIRS = [
  ["n", "ratio"],
  ["n", "cellRatio"],
  ["n", "widthStart"],
  ["n", "widthEnd"],
  ["ratio", "cellRatio"],
  ["ratio", "widthStart"],
  ["ratio", "widthEnd"],
  ["cellRatio", "widthStart"],
  ["cellRatio", "widthEnd"],
  ["widthStart", "widthEnd"],
];

export function renderInspector() {
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

export function renderPeriodicInputs() {
  els.periodicGrid.innerHTML = `
    <span>Periodic</span>
    ${AXES.map(
      (axis, index) => `
        <label class="check periodic-check">
          <input type="checkbox" data-axis-index="${index}" ${state.project.periodicAxes[index] ? "checked" : ""} />
          <span>${axis.toUpperCase()}</span>
        </label>
      `
    ).join("")}
  `;
  for (const input of els.periodicGrid.querySelectorAll("input")) {
    input.onchange = () => {
      state.project.periodicAxes[Number(input.dataset.axisIndex)] = input.checked;
      updateStructure();
    };
  }
}

export function renderGeometryInputs(block) {
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
    const field = input.dataset.field;
    const index = Number(input.dataset.index);
    const counts = field === "ng" || field === "dims";
    bindNumberField(input, {
      integer: counts,
      valid: (value) => !counts || value >= 1,
      read: () => block[field][index],
      write: (value) => {
        const oldMin = block.lmin[index];
        const oldMax = block.lmax[index];
        block[field][index] = value;
        if (field === "lmin" || field === "lmax") {
          rescaleExplicitFaces(block, index, oldMin, oldMax, block.lmin[index], block.lmax[index]);
        }
        if (field === "dims") propagateMpiPartition(state.project.blocks, block.id, index, value);
        markProjectDirty(field === "ng" || field === "lmin" || field === "lmax" ? block.id : null);
        if (field === "ng" || field === "lmin" || field === "lmax") schedulePreview(block.id);
        renderScene();
        renderBlockList();
      },
    });
    if (field === "dims") {
      input.addEventListener("change", () => updateStructure({ silent: true }));
    }
    if (field === "lmin" || field === "lmax") input.addEventListener("change", keepModelInView);
  }
}

export function renderAxisTabs() {
  for (const button of els.axisTabs.querySelectorAll("button")) {
    button.classList.toggle("active", button.dataset.axis === state.currentAxis);
    button.onclick = () => {
      state.currentAxis = button.dataset.axis;
      renderInspector();
      renderSpacingPreview();
      renderScene();
    };
  }
}

export function renderAxisEditor(block) {
  const axis = ensureAxis(block, state.currentAxis);
  const axisIndex = AXES.indexOf(state.currentAxis);
  const locked = block.axisLocks[axisIndex];
  const selectedKind = axis.kind ?? "native";
  els.axisEditor.innerHTML = `
    <label class="field">
      <span>Mode</span>
      <select id="grid-kind">
        <option value="native">Native SNaC</option>
        <option value="monotone">Monotone grading</option>
        <option value="symmetric">Symmetric grading</option>
        <option value="multi">Multi-region grading</option>
        ${selectedKind === "explicit" ? '<option value="explicit">Explicit coordinates</option>' : ""}
      </select>
    </label>
    <div id="grid-kind-body"></div>
    <div class="axis-actions">
      <button id="axis-lock" class="button icon${locked ? " active" : ""}" type="button" title="${locked ? "Unlock" : "Lock"} ${state.currentAxis.toUpperCase()} grid"><i data-lucide="${locked ? "lock" : "unlock"}"></i></button>
      <button id="apply-axis" class="button" type="button"><i data-lucide="copy-check"></i><span>Apply to aligned blocks</span></button>
    </div>
  `;
  const kindSelect = document.getElementById("grid-kind");
  kindSelect.value = selectedKind;
  kindSelect.onchange = () => {
    axis.kind = kindSelect.value;
    if (axis.kind === "monotone") {
      axis.controls = axis.controls?.length === 2 ? axis.controls : ["n", "ratio"];
      setMonotoneControlsFromPreview(block, axis, axis.controls);
    }
    if (axis.kind === "symmetric" && !["both", "middle"].includes(axis.side)) axis.side = "both";
    renderAxisEditor(block);
    renderGeometryInputs(block);
    gridDefinitionChanged(block);
  };
  renderAxisKindBody(block, axis);
  document.getElementById("axis-lock").onclick = () => {
    block.axisLocks[axisIndex] = !block.axisLocks[axisIndex];
    markProjectDirty();
    renderAxisEditor(block);
    if (window.lucide) window.lucide.createIcons();
  };
  document.getElementById("apply-axis").onclick = applyAxisToAlignedBlocks;
}

export function profileSelect(axis) {
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

export function bindProfileSelect(block, axis) {
  const select = document.getElementById("grid-profile");
  if (!select) return;
  select.onchange = () => {
    axis.profile = select.value;
    if (axis.kind === "monotone") {
      if (axis.profile !== "geometric" && axis.controls?.includes("cellRatio")) {
        setMonotoneControlsFromPreview(block, axis, ["n", "ratio"]);
      }
      renderAxisKindBody(block, axis);
    }
    gridDefinitionChanged(block);
  };
}

export function renderAxisKindBody(block, axis) {
  const body = document.getElementById("grid-kind-body");
  const preview = previewCache.get(block.id)?.[state.currentAxis];
  if (axis.kind === "explicit") {
    const faces = axis.faces ?? [];
    const spacings = faces.slice(1).map((face, index) => face - faces[index]);
    const minimum = spacings.reduce((value, spacing) => Math.min(value, spacing), Number.POSITIVE_INFINITY);
    const maximum = spacings.reduce((value, spacing) => Math.max(value, spacing), Number.NEGATIVE_INFINITY);
    body.innerHTML = `
      <div class="explicit-grid-summary">
        <span>cells <strong>${Math.max(0, faces.length - 1)}</strong></span>
        <span>minimum <strong>${formatNumber(minimum)}</strong></span>
        <span>maximum <strong>${formatNumber(maximum)}</strong></span>
      </div>
    `;
  } else if (axis.kind === "multi") {
    axis.segments = axis.segments?.length ? axis.segments : [{ length: 1, cells: 1, ratio: 1, continuous: false }];
    body.innerHTML = `
      ${profileSelect(axis)}
      <div class="segment-table">
        ${axis.segments
          .map(
            (segment, index) => `
          <div class="segment-row" data-segment="${index}">
            <label class="field"><span>Length weight</span><input data-key="length" type="number" step="any" min="1e-12" value="${escapeHtml(segment.length)}"></label>
            <label class="field"><span>Cell weight</span><input data-key="cells" type="number" step="any" min="1e-12" value="${escapeHtml(segment.cells)}"></label>
            <label class="field"><span>End/start</span><input data-key="ratio" type="number" step="any" min="1e-12" value="${escapeHtml(segment.ratio)}" ${index > 0 && segment.continuous ? "disabled" : ""}></label>
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
        const key = input.dataset.key;
        bindNumberField(input, {
          valid: (value) => value > 0,
          read: () => axis.segments[index][key],
          write: (value) => {
            axis.segments[index][key] = value;
            gridDefinitionChanged(block);
          },
        });
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
  } else if (axis.kind === "monotone") {
    normalizeMonotoneAxis(axis);
    const pairs = MONOTONE_CONTROL_PAIRS.filter((pair) => axis.profile === "geometric" || !pair.includes("cellRatio"));
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
      const key = input.dataset.control;
      const axisIndex = AXES.indexOf(state.currentAxis);
      bindNumberField(input, {
        integer: key === "n",
        valid: (value) => (key === "n" ? value >= 1 : value > 0),
        read: () => (key === "n" ? block.ng[axisIndex] : axis[key]),
        write: (value) => updateMonotoneControl(block, axis, key, value),
      });
    }
  } else if (axis.kind === "symmetric") {
    body.innerHTML = `
      ${profileSelect(axis)}
      <div class="form-grid two">
        <label class="field"><span>Minimum spacing</span><input id="min-dx" type="number" step="any" min="1e-12" value="${escapeHtml(axis.min ?? 0.01)}"></label>
        <label class="field"><span>Maximum spacing</span><input id="max-dx" type="number" step="any" min="1e-12" value="${escapeHtml(axis.max ?? 0.02)}"></label>
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
      const key = id === "min-dx" ? "min" : "max";
      bindNumberField(document.getElementById(id), {
        valid: (value) => value > 0,
        read: () => axis[key],
        write: (value) => {
          axis[key] = value;
          gridDefinitionChanged(block);
        },
      });
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
      <label class="field"><span>gr</span><input id="snac-gr" type="number" step="any" min="0" value="${escapeHtml(axis.gr ?? 0)}"></label>
    `;
    document.getElementById("snac-gt").value = axis.gt ?? 0;
    document.getElementById("snac-gt").onchange = (event) => {
      const previousUniformGr = GEOMETRIC_UNIFORM_GR[axis.gt];
      axis.gt = Number(event.target.value);
      // The geometric functions divide by gr; start them from their uniform
      // value, and drop that value again when returning to a tanh function.
      const uniformGr = GEOMETRIC_UNIFORM_GR[axis.gt];
      if (uniformGr && !(axis.gr > 0)) axis.gr = uniformGr;
      if (!uniformGr && previousUniformGr && axis.gr === previousUniformGr) axis.gr = 0;
      document.getElementById("snac-gr").value = String(axis.gr ?? 0);
      gridDefinitionChanged(block);
    };
    bindNumberField(document.getElementById("snac-gr"), {
      valid: (value) => value >= 0,
      read: () => axis.gr ?? 0,
      write: (value) => {
        axis.gr = value;
        gridDefinitionChanged(block);
      },
    });
  }
  if (window.lucide) window.lucide.createIcons();
}

export function normalizeMonotoneAxis(axis) {
  axis.controls = Array.isArray(axis.controls) && axis.controls.length === 2 ? axis.controls : ["n", "ratio"];
  if (axis.profile !== "geometric" && axis.controls.includes("cellRatio")) axis.controls = ["n", "ratio"];
  axis.ratio = positiveNumber(axis.ratio, 1);
  axis.cellRatio = positiveNumber(axis.cellRatio, 1);
  axis.widthStart = positiveNumber(axis.widthStart, 0.01);
  axis.widthEnd = positiveNumber(axis.widthEnd, 0.02);
  axis.side = "end";
}

export function monotoneControlInput(block, axis, key) {
  const axisIndex = AXES.indexOf(state.currentAxis);
  const value = key === "n" ? block.ng[axisIndex] : axis[key];
  const step = key === "n" ? "1" : "any";
  const min = key === "n" ? "1" : "1e-12";
  return `<label class="field"><span>${MONOTONE_CONTROL_LABELS[key]}</span><input data-control="${key}" type="number" step="${step}" min="${min}" value="${escapeHtml(value)}"></label>`;
}

export function updateMonotoneControl(block, axis, key, value) {
  if (key === "n") {
    const axisIndex = AXES.indexOf(state.currentAxis);
    block.ng[axisIndex] = value;
    const overview = els.resolutionGrid.querySelector(`input[data-field="ng"][data-index="${axisIndex}"]`);
    if (overview) overview.value = block.ng[axisIndex];
  } else {
    axis[key] = value;
  }
  gridDefinitionChanged(block);
}

export function setMonotoneControlsFromPreview(block, axis, controlKeys) {
  const preview = previewCache.get(block.id)?.[state.currentAxis];
  axis.controls = controlKeys;
  const axisIndex = AXES.indexOf(state.currentAxis);
  const values = {
    n: preview?.n ?? block.ng[axisIndex],
    ratio: preview?.expansion ?? axis.ratio ?? 1,
    cellRatio: preview?.cellRatio ?? axis.cellRatio ?? 1,
    widthStart: preview?.lower ?? axis.widthStart ?? 0.01,
    widthEnd: preview?.upper ?? axis.widthEnd ?? 0.02,
  };
  for (const key of controlKeys) {
    if (key === "n") block.ng[axisIndex] = Math.max(1, Math.round(values.n));
    else axis[key] = positiveNumber(values[key], axis[key] ?? 1);
  }
  axis.side = "end";
}

export function segmentMetricMarkup(metric, fallbackIndex = "") {
  const index = metric?.index ?? fallbackIndex;
  if (!metric) return `<div class="segment-metrics" data-segment-metric="${index}"></div>`;
  return `<div class="segment-metrics" data-segment-metric="${index}">
    <span>N <strong>${metric.cells}</strong></span>
    <span>lower <strong>${formatNumber(metric.widthStart)}</strong></span>
    <span>upper <strong>${formatNumber(metric.widthEnd)}</strong></span>
    <span>jump <strong>${formatNumber(metric.jump)}</strong></span>
  </div>`;
}

export function updateSegmentMetrics(preview) {
  for (const metric of preview?.segments ?? []) {
    const container = els.axisEditor.querySelector(`[data-segment-metric="${metric.index}"]`);
    if (container) container.outerHTML = segmentMetricMarkup(metric);
    if (metric.continuous) {
      const ratio = els.axisEditor.querySelector(`.segment-row[data-segment="${metric.index}"] input[data-key="ratio"]`);
      if (ratio) ratio.value = formatNumber(metric.ratio);
    }
  }
}

export function renderBoundaryTable(block) {
  ensureBoundaryArrays(block, state.project.nscal);
  els.boundaryTable.innerHTML = FACE_ORDER.map((face, index) => {
    return `
      <div class="boundary-row">
        <span class="axis-label">${face.toUpperCase()}</span>
        <div class="boundary-components">
          ${["u", "v", "w"]
            .map((name, component) => boundaryComponent(name, index, block.cbcvel[component][index], block.bcvel[component][index], component))
            .join("")}
          ${boundaryComponent("p", index, block.cbcpre[index], block.bcpre[index])}
          ${Array.from({ length: state.project.nscal }, (_, iscal) =>
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
    const face = Number(input.dataset.face);
    const component = Number(input.dataset.component);
    const row = () => {
      if (input.dataset.kind === "pre") return block.bcpre;
      return input.dataset.kind === "scal" ? block.bcscal[component] : block.bcvel[component];
    };
    bindNumberField(input, {
      read: () => row()[face],
      write: (value) => {
        row()[face] = value;
        markProjectDirty();
      },
    });
  }
}

export function numberInput(field, index, value, min = null, disabled = false) {
  const minAttr = min == null ? "" : `min="${min}"`;
  const disabledAttr = disabled ? 'disabled title="Derived from min/max spacing"' : "";
  return `<label class="field"><input data-field="${field}" data-index="${index}" type="number" step="any" ${minAttr} ${disabledAttr} value="${escapeHtml(value)}"></label>`;
}

export function boundaryComponent(name, face, code, value, component = null, forcedKind = null) {
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

export function bcSelect(kind, face, selected, componentAttr = "") {
  return `<select data-kind="${kind}" data-face="${face}" ${componentAttr}>
    ${["D", "N", "F"].map((value) => `<option value="${value}" ${value === selected ? "selected" : ""}>${value}</option>`).join("")}
  </select>`;
}

export function bcValueInput(kind, face, value, componentAttr = "") {
  return `
    <input data-kind="${kind}" data-face="${face}" ${componentAttr} type="number" step="any" value="${escapeHtml(value ?? 0)}">
  `;
}

export function updateGeometryInputValues(block) {
  for (const field of ["lmin", "lmax"]) {
    for (const input of els.geometryGrid.querySelectorAll(`input[data-field="${field}"]`)) {
      input.value = formatInputNumber(block[field][Number(input.dataset.index)]);
    }
  }
}

export function renderSpacingPreview() {
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
  const preview = previewCache.get(block.id)?.[state.currentAxis];
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
