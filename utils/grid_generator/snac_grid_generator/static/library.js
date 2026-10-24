export const LIBRARY_SCHEMA_VERSION = 1;
export const LIBRARY_STORAGE_KEY = "snac-grid-generator-library-v1";

const ITEM_KINDS = ["axisPresets", "projectTemplates"];

export function emptyLibrary() {
  return {
    schemaVersion: LIBRARY_SCHEMA_VERSION,
    axisPresets: [],
    projectTemplates: [],
  };
}

export function loadLibrary(storage = localStorage) {
  try {
    const text = storage.getItem(LIBRARY_STORAGE_KEY);
    return text ? normalizeLibrary(JSON.parse(text)) : emptyLibrary();
  } catch (_error) {
    return emptyLibrary();
  }
}

export function saveLibrary(library, storage = localStorage) {
  const normalized = normalizeLibrary(library);
  storage.setItem(LIBRARY_STORAGE_KEY, JSON.stringify(normalized));
  return normalized;
}

export function normalizeLibrary(value) {
  if (!value || typeof value !== "object" || Array.isArray(value)) throw new Error("Library JSON must be an object");
  const version = Number(value.schemaVersion ?? 1);
  if (version !== LIBRARY_SCHEMA_VERSION) throw new Error(`Unsupported library schema version ${version}`);
  const result = emptyLibrary();
  for (const kind of ITEM_KINDS) {
    if (!Array.isArray(value[kind] ?? [])) throw new Error(`${kind} must be an array`);
    result[kind] = (value[kind] ?? []).map((item) => normalizeItem(kind, item));
    ensureUniqueIds(result[kind], kind);
  }
  return result;
}

export function mergeLibraries(current, incoming) {
  const merged = normalizeLibrary(current);
  const addition = normalizeLibrary(incoming);
  const counts = { axisPresets: 0, projectTemplates: 0 };
  for (const kind of ITEM_KINDS) {
    const byId = new Map(merged[kind].map((item) => [item.id, item]));
    for (const item of addition[kind]) {
      byId.set(item.id, item);
      counts[kind] += 1;
    }
    merged[kind] = [...byId.values()];
  }
  return { library: merged, counts };
}

export function captureAxisPreset(name, axis, ng, lower, upper, id = createId("axis")) {
  const cleanName = normalizeName(name);
  const cells = positiveInteger(ng, "Axis preset cell count");
  const cleanAxis = cloneJson(axis, "Axis preset");
  const coordinateSpace = cleanAxis.kind === "explicit" ? "relative" : "native";
  if (coordinateSpace === "relative") {
    const extent = finiteExtent(lower, upper);
    if (!Array.isArray(cleanAxis.faces) || cleanAxis.faces.length !== cells + 1) {
      throw new Error("Explicit axis preset coordinates must contain one more face than cells");
    }
    cleanAxis.faces = cleanAxis.faces.map((value) => finiteNumber(value, "Explicit face"));
    const tolerance = Math.max(1, Math.abs(extent.lower), Math.abs(extent.upper)) * 1e-10;
    if (
      Math.abs(cleanAxis.faces[0] - extent.lower) > tolerance ||
      Math.abs(cleanAxis.faces[cleanAxis.faces.length - 1] - extent.upper) > tolerance
    ) {
      throw new Error("Explicit axis preset coordinates must span the source block");
    }
    cleanAxis.faces = cleanAxis.faces.map((value) => (value - extent.lower) / extent.length);
    validateIncreasingFaces(cleanAxis.faces);
    cleanAxis.faces[0] = 0;
    cleanAxis.faces[cleanAxis.faces.length - 1] = 1;
  }
  return normalizeItem("axisPresets", { id, name: cleanName, ng: cells, coordinateSpace, axis: cleanAxis });
}

export function applyAxisPreset(preset, lower, upper) {
  const item = normalizeItem("axisPresets", preset);
  const axis = cloneJson(item.axis, "Axis preset");
  if (item.coordinateSpace === "relative") {
    const extent = finiteExtent(lower, upper);
    axis.faces = axis.faces.map((value) => extent.lower + value * extent.length);
    axis.faces[0] = extent.lower;
    axis.faces[axis.faces.length - 1] = extent.upper;
  }
  return { axis, ng: item.ng };
}

export function captureProjectTemplate(name, project, id = createId("case")) {
  return normalizeItem("projectTemplates", {
    id,
    name: normalizeName(name),
    project: cloneJson(project, "Case template"),
  });
}

export function renameLibraryItem(library, kind, id, name) {
  const normalized = normalizeLibrary(library);
  const item = normalized[libraryKind(kind)].find((candidate) => candidate.id === id);
  if (!item) throw new Error("Library item no longer exists");
  item.name = normalizeName(name);
  return normalized;
}

export function removeLibraryItem(library, kind, id) {
  const normalized = normalizeLibrary(library);
  const key = libraryKind(kind);
  normalized[key] = normalized[key].filter((item) => item.id !== id);
  return normalized;
}

export function replaceLibraryItem(library, kind, item) {
  const normalized = normalizeLibrary(library);
  const key = libraryKind(kind);
  const cleanItem = normalizeItem(key, item);
  const index = normalized[key].findIndex((candidate) => candidate.id === cleanItem.id);
  if (index < 0) normalized[key].push(cleanItem);
  else normalized[key][index] = cleanItem;
  return normalized;
}

function normalizeItem(kind, value) {
  if (!value || typeof value !== "object" || Array.isArray(value)) throw new Error("Library items must be objects");
  const item = { id: normalizeId(value.id), name: normalizeName(value.name) };
  if (kind === "axisPresets") {
    item.ng = positiveInteger(value.ng, "Axis preset cell count");
    item.axis = cloneJson(value.axis, "Axis preset");
    if (!item.axis || typeof item.axis !== "object" || Array.isArray(item.axis)) throw new Error("Axis preset data must be an object");
    item.coordinateSpace = value.coordinateSpace === "relative" ? "relative" : "native";
    if (item.axis.kind === "explicit") {
      if (!Array.isArray(item.axis.faces) || item.axis.faces.length !== item.ng + 1) {
        throw new Error("Explicit axis preset coordinates must contain one more face than cells");
      }
      item.axis.faces = item.axis.faces.map((face) => finiteNumber(face, "Explicit face"));
      validateIncreasingFaces(item.axis.faces);
      if (item.coordinateSpace === "relative" && (
        Math.abs(item.axis.faces[0]) > 1e-10 || Math.abs(item.axis.faces[item.axis.faces.length - 1] - 1) > 1e-10
      )) {
        throw new Error("Relative axis preset coordinates must span zero to one");
      }
    }
    return item;
  }
  if (kind === "projectTemplates") {
    item.project = cloneJson(value.project, "Case template");
    if (!item.project || typeof item.project !== "object" || !Array.isArray(item.project.blocks)) {
      throw new Error("Case template must contain a project with a blocks array");
    }
    return item;
  }
  throw new Error(`Unknown library item kind ${kind}`);
}

function libraryKind(kind) {
  if (!ITEM_KINDS.includes(kind)) throw new Error(`Unknown library item kind ${kind}`);
  return kind;
}

function createId(prefix) {
  const suffix = globalThis.crypto?.randomUUID?.() ?? `${Date.now()}-${Math.random().toString(16).slice(2)}`;
  return `${prefix}-${suffix}`;
}

function normalizeId(value) {
  const id = String(value ?? "").trim();
  if (!id || id.length > 160) throw new Error("Library item id is invalid");
  return id;
}

function normalizeName(value) {
  const name = String(value ?? "").trim();
  if (!name) throw new Error("Library item name is required");
  if (name.length > 80) throw new Error("Library item name must be 80 characters or fewer");
  return name;
}

function positiveInteger(value, label) {
  const number = Math.round(Number(value));
  if (!Number.isFinite(number) || number < 1) throw new Error(`${label} must be positive`);
  return number;
}

function finiteNumber(value, label) {
  const number = Number(value);
  if (!Number.isFinite(number)) throw new Error(`${label} must be finite`);
  return number;
}

function finiteExtent(lower, upper) {
  const lo = finiteNumber(lower, "Axis lower bound");
  const hi = finiteNumber(upper, "Axis upper bound");
  if (hi <= lo) throw new Error("Axis upper bound must exceed its lower bound");
  return { lower: lo, upper: hi, length: hi - lo };
}

function validateIncreasingFaces(faces) {
  for (let index = 1; index < faces.length; index += 1) {
    if (faces[index] <= faces[index - 1]) throw new Error("Axis preset coordinates must be strictly increasing");
  }
}

function cloneJson(value, label) {
  try {
    return JSON.parse(JSON.stringify(value));
  } catch (_error) {
    throw new Error(`${label} must contain JSON-compatible data`);
  }
}

function ensureUniqueIds(items, kind) {
  const ids = new Set();
  for (const item of items) {
    if (ids.has(item.id)) throw new Error(`${kind} contains duplicate id ${item.id}`);
    ids.add(item.id);
  }
}
