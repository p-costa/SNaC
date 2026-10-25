// Plain-language names and descriptions of axis gradings, shared by the
// inspector, the grid library, and the repair dialog.

import { canonicalizeAxis } from "./block-model.js";
import { formatCompact } from "./format.js";

export const MONOTONE_CONTROL_LABELS = {
  n: "Cells",
  ratio: "Upper/lower ratio",
  cellRatio: "Cell-to-cell ratio",
  widthStart: "Lower spacing",
  widthEnd: "Upper spacing",
};

const KIND_LABELS = {
  native: "Native SNaC",
  monotone: "Monotone grading",
  symmetric: "Symmetric grading",
  multi: "Multi-region grading",
  explicit: "Explicit coordinates",
};

export function axisKindLabel(kind) {
  return KIND_LABELS[canonicalizeAxis({ kind }).kind] ?? String(kind ?? "Axis grid");
}

// For example "monotone tanh, lower spacing 0.005, 24 cells".
export function describeAxis(axis, cells) {
  const spec = canonicalizeAxis({ ...(axis ?? {}) });
  const parts = [];
  if (spec.kind === "native") {
    const uniform = (spec.gt ?? 0) < 4 && !(spec.gr > 0);
    parts.push(uniform ? "uniform" : `native function ${spec.gt}, gr ${formatCompact(spec.gr)}`);
  } else if (spec.kind === "monotone") {
    parts.push(`monotone ${spec.profile ?? "geometric"}`);
    for (const key of spec.controls ?? ["n", "ratio"]) {
      if (key !== "n") parts.push(`${MONOTONE_CONTROL_LABELS[key].toLowerCase()} ${formatCompact(spec[key])}`);
    }
  } else if (spec.kind === "symmetric") {
    const where = spec.side === "middle" ? "middle" : "both ends";
    parts.push(`symmetric ${spec.profile ?? "geometric"}, spacing ${formatCompact(spec.min)} to ${formatCompact(spec.max)} at ${where}`);
  } else if (spec.kind === "multi") {
    const count = spec.segments?.length ?? 1;
    parts.push(`multi-region, ${count} segment${count === 1 ? "" : "s"}`);
  } else if (spec.kind === "explicit") {
    parts.push("explicit coordinates");
  } else {
    parts.push(axisKindLabel(spec.kind).toLowerCase());
  }
  if (Number.isFinite(cells)) parts.push(`${cells} cells`);
  return parts.join(", ");
}
