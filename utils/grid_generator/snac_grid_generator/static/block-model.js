export const AXES = ["x", "y", "z"];
// Mirrors DEFAULT_MIN_LOCAL_CELLS in model.py.
export const DEFAULT_MIN_LOCAL_CELLS = 4;
export const FACE_ORDER = ["x-", "x+", "y-", "y+", "z-", "z+"];

export function blockSize(block) {
  return block.lmax.map((value, index) => Math.max(1e-9, value - block.lmin[index]));
}

export function blockCenter(block) {
  return block.lmin.map((value, index) => 0.5 * (value + block.lmax[index]));
}

// Legacy axis names, as in axis_format.py: older project files, autosaves, and
// library items may still use them.
const LEGACY_KINDS = { snac: "native", simple_ratio: "monotone", max_min: "symmetric" };
const LEGACY_KEYS = { cell_ratio: "cellRatio", width_start: "widthStart", width_end: "widthEnd" };

export function canonicalizeAxis(axis) {
  if (!axis) return axis;
  axis.kind = LEGACY_KINDS[axis.kind] ?? axis.kind;
  for (const [legacy, key] of Object.entries(LEGACY_KEYS)) {
    if (!(legacy in axis)) continue;
    if (!(key in axis)) axis[key] = axis[legacy];
    delete axis[legacy];
  }
  if (Array.isArray(axis.controls)) axis.controls = axis.controls.map((name) => LEGACY_KEYS[name] ?? name);
  return axis;
}

export function ensureAxis(block, axis) {
  block.axes = block.axes || {};
  block.axes[axis] = canonicalizeAxis(block.axes[axis] || defaultAxis());
  block.axes[axis].segments = block.axes[axis].segments || [defaultSegment()];
  for (const segment of block.axes[axis].segments) segment.continuous = Boolean(segment.continuous);
  block.axes[axis].profile = block.axes[axis].profile || "geometric";
  block.axes[axis].controls = block.axes[axis].controls || ["n", "ratio"];
  return block.axes[axis];
}

export function isCellCountDerived(axis) {
  const kind = LEGACY_KINDS[axis.kind] ?? axis.kind;
  return kind === "explicit" || kind === "symmetric" || (kind === "monotone" && !axis.controls?.includes("n"));
}

export function rescaleExplicitFaces(block, axisIndex, oldMin, oldMax, newMin, newMax, sourceFaces = null) {
  const axis = ensureAxis(block, AXES[axisIndex]);
  if (axis.kind !== "explicit" || !axis.faces?.length || oldMax <= oldMin || newMax <= newMin) return;
  axis.faces = (sourceFaces ?? axis.faces).map(
    (face) => newMin + ((face - oldMin) / (oldMax - oldMin)) * (newMax - newMin)
  );
  axis.faces[0] = newMin;
  axis.faces[axis.faces.length - 1] = newMax;
}

export function defaultAxis() {
  return {
    kind: "native",
    gt: 0,
    gr: 0,
    ratio: 1,
    cellRatio: 1,
    widthStart: 0.01,
    widthEnd: 0.02,
    controls: ["n", "ratio"],
    profile: "geometric",
    min: 0.01,
    max: 0.02,
    side: "end",
    segments: [defaultSegment()],
  };
}

export function defaultBlock(id) {
  return {
    id,
    name: "",
    dims: [1, 1, 1],
    ng: [32, 32, 2],
    lmin: [0, 0, 0],
    lmax: [1, 1, 0.05],
    axes: { x: defaultAxis(), y: defaultAxis(), z: defaultAxis() },
    axisLocks: [false, false, false],
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

export function ensureBoundaryArrays(block, nscal) {
  const count = normalizedScalarCount(nscal);
  block.cbcvel = Array.from({ length: 3 }, (_, component) => completeArray(block.cbcvel?.[component], "D"));
  block.cbcpre = completeArray(block.cbcpre, "N");
  block.bcvel = Array.from({ length: 3 }, (_, component) => completeArray(block.bcvel?.[component], 0));
  block.bcpre = completeArray(block.bcpre, 0);
  block.cbcscal = Array.from({ length: count }, (_, scalar) => completeArray(block.cbcscal?.[scalar], "N"));
  block.bcscal = Array.from({ length: count }, (_, scalar) => completeArray(block.bcscal?.[scalar], 0));
  block.inflow = completeArray(block.inflow, 0);
}

export function normalizedPeriodicAxes(project) {
  const values = Array.isArray(project.periodicAxes) ? project.periodicAxes : project.periodic_axes ?? [];
  return AXES.map((_, index) => Boolean(values[index]));
}

export function normalizedFftAxis(project) {
  const value = String(project.fftAxis ?? project.fft_axis ?? "").trim().toLowerCase();
  // SNaC's FFT_AXIS build values 1 to 3 name the same axes.
  const axis = { 1: "x", 2: "y", 3: "z" }[value] ?? value;
  return AXES.includes(axis) ? axis : null;
}

export function normalizedDecomposition(project) {
  const value = project.decomposition ?? {};
  const targetRanks = Math.max(0, Math.round(Number(value.targetRanks ?? value.target_ranks) || 0));
  const mode = ["auto", "1d", "2d", "3d"].includes(value.mode) ? value.mode : "auto";
  const sourceAxes = Array.isArray(value.axes) ? value.axes : [true, true, true];
  const axes = AXES.map((_, index) => (index < sourceAxes.length ? Boolean(sourceAxes[index]) : true));
  // SNaC never splits the FFT axis, so it is never a partition axis.
  const fftIndex = AXES.indexOf(normalizedFftAxis(project));
  if (fftIndex >= 0) axes[fftIndex] = false;
  if (!axes.some(Boolean)) axes.forEach((_, index) => { axes[index] = index !== fftIndex; });
  const minLocalCells = Math.max(1, Math.round(Number(value.minLocalCells ?? value.min_local_cells) || DEFAULT_MIN_LOCAL_CELLS));
  const maxLocalAspect = Math.max(0, Number(value.maxLocalAspect ?? value.max_local_aspect) || 0);
  return { targetRanks, mode, axes, minLocalCells, maxLocalAspect };
}

export function normalizedAxisLocks(values) {
  const source = Array.isArray(values) ? values : [];
  return AXES.map((_, index) => Boolean(source[index]));
}

export function normalizedScalarCount(value) {
  return Math.max(0, Math.round(Number(value) || 0));
}

export function resetFriendBoundaries(block, nscal) {
  block.cbcvel = Array.from({ length: 3 }, () => ["D", "D", "D", "D", "D", "D"]);
  block.cbcpre = ["N", "N", "N", "N", "N", "N"];
  block.bcvel = Array.from({ length: 3 }, () => [0, 0, 0, 0, 0, 0]);
  block.bcpre = [0, 0, 0, 0, 0, 0];
  block.cbcscal = Array.from({ length: normalizedScalarCount(nscal) }, () => ["N", "N", "N", "N", "N", "N"]);
  block.bcscal = Array.from({ length: normalizedScalarCount(nscal) }, () => [0, 0, 0, 0, 0, 0]);
}

export function renumberBlocksByOrder(blocks) {
  const idMap = new Map(blocks.map((block, index) => [block.id, index + 1]));
  for (const block of blocks) {
    block.id = idMap.get(block.id);
    block.name = `block-${block.id}`;
  }
  for (const block of blocks) {
    remapFriendRows(block.cbcvel, block.bcvel, idMap);
    remapFriendRow(block.cbcpre, block.bcpre, idMap);
    remapFriendRows(block.cbcscal, block.bcscal, idMap);
  }
  return idMap;
}

export function propagateMpiPartition(blocks, sourceId, axisIndex, value) {
  const graph = new Map(blocks.map((block) => [block.id, new Set()]));
  for (const connection of fullFaceConnections(blocks)) {
    if (connection.axisIndex === axisIndex) continue;
    graph.get(connection.a.id).add(connection.b.id);
    graph.get(connection.b.id).add(connection.a.id);
  }
  const queue = [sourceId];
  const seen = new Set(queue);
  while (queue.length) {
    const blockId = queue.shift();
    const block = blocks.find((item) => item.id === blockId);
    if (block) block.dims[axisIndex] = value;
    for (const neighborId of graph.get(blockId) ?? []) {
      if (!seen.has(neighborId)) {
        seen.add(neighborId);
        queue.push(neighborId);
      }
    }
  }
}

export function fullFaceConnections(blocks) {
  const connections = [];
  for (let i = 0; i < blocks.length; i += 1) {
    for (let j = i + 1; j < blocks.length; j += 1) {
      const a = blocks[i];
      const b = blocks[j];
      for (let axisIndex = 0; axisIndex < 3; axisIndex += 1) {
        const scale = Math.max(axisExtent(a, axisIndex), axisExtent(b, axisIndex));
        if (touches(a.lmax[axisIndex], b.lmin[axisIndex], scale) && sameCrossSection(a, b, axisIndex)) {
          connections.push({ a, b, axisIndex });
        }
        if (touches(b.lmax[axisIndex], a.lmin[axisIndex], scale) && sameCrossSection(a, b, axisIndex)) {
          connections.push({ a: b, b: a, axisIndex });
        }
      }
    }
  }
  return connections;
}

export function nextBlockId(blocks) {
  return Math.max(0, ...blocks.map((block) => block.id)) + 1;
}

function defaultSegment() {
  return { length: 1, cells: 1, ratio: 1, continuous: false };
}

function completeArray(values, fallback) {
  return Array.from({ length: 6 }, (_, index) => values?.[index] ?? fallback);
}

function remapFriendRows(codeRows, valueRows, idMap) {
  for (let index = 0; index < (codeRows?.length ?? 0); index += 1) {
    remapFriendRow(codeRows[index], valueRows?.[index], idMap);
  }
}

function remapFriendRow(codes, values, idMap) {
  if (!codes || !values) return;
  for (let face = 0; face < codes.length; face += 1) {
    if (codes[face] !== "F") continue;
    const friendId = idMap.get(Math.round(Number(values[face])));
    if (friendId != null) values[face] = friendId;
  }
}

function sameCrossSection(a, b, axisIndex) {
  for (let index = 0; index < 3; index += 1) {
    if (index === axisIndex) continue;
    const scale = Math.max(axisExtent(a, index), axisExtent(b, index));
    if (!touches(a.lmin[index], b.lmin[index], scale) || !touches(a.lmax[index], b.lmax[index], scale)) {
      return false;
    }
  }
  return true;
}

function axisExtent(block, axisIndex) {
  return Math.abs(block.lmax[axisIndex] - block.lmin[axisIndex]);
}

function touches(a, b, scale) {
  const roundoff = 8 * Number.EPSILON * Math.max(Math.abs(a), Math.abs(b), Math.abs(scale));
  return Math.abs(a - b) <= Math.max(1e-10 * Math.abs(scale), roundoff);
}
