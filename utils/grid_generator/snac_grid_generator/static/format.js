// Pure formatting, parsing, escaping, and download helpers.

export function positiveNumber(value, fallback) {
  const number = Number(value);
  return Number.isFinite(number) && number > 0 ? number : fallback;
}

export function formatNumber(value) {
  if (!Number.isFinite(value)) return "n/a";
  if (Math.abs(value) >= 1000 || Math.abs(value) < 0.001) return value.toExponential(3);
  return value.toPrecision(5);
}

export function formatInteger(value) {
  return Number.isFinite(value) ? Math.round(value).toLocaleString() : "n/a";
}

export function formatPercent(value) {
  if (!Number.isFinite(value)) return "n/a";
  return `${(100 * value).toPrecision(value < 0.001 ? 2 : 3)}%`;
}

export function formatBytes(value) {
  let size = Number(value) || 0;
  const units = ["B", "KiB", "MiB", "GiB"];
  let index = 0;
  while (size >= 1024 && index < units.length - 1) {
    size /= 1024;
    index += 1;
  }
  return `${formatNumber(size)} ${units[index]}`;
}

export function formatInputNumber(value) {
  return Number.isFinite(value) ? value.toString() : "";
}

// Bind a number field so that empty, partial, or out-of-range input never reaches
// the model. While typing, the last valid value is kept; when the field loses
// focus with an invalid entry, it shows that value again.
export function bindNumberField(input, { read, write, valid = () => true, integer = false }) {
  const parse = () => {
    const parsed = parsedNumber(input.value);
    const value = integer ? Math.round(parsed) : parsed;
    return Number.isFinite(value) && valid(value) ? value : null;
  };
  input.addEventListener("input", () => {
    const value = parse();
    if (value !== null) write(value);
  });
  input.addEventListener("change", () => {
    if (parse() === null) input.value = String(read());
  });
}

// Up to five significant digits without trailing zeros, for prose: 4, 0.005.
export function formatCompact(value) {
  return Number.isFinite(value) ? String(Number(value.toPrecision(5))) : "n/a";
}

export function parsedNumber(text) {
  const trimmed = String(text).trim();
  return trimmed === "" ? Number.NaN : Number(trimmed);
}

export function escapeHtml(value) {
  return String(value).replace(/[&<>"']/g, (char) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#039;" })[char]);
}

export function toCamel(id) {
  return id.replace(/-([a-z])/g, (_, char) => char.toUpperCase());
}

export function downloadText(content, filename, type) {
  const blob = new Blob([content], { type });
  const link = document.createElement("a");
  link.href = URL.createObjectURL(blob);
  link.download = filename;
  link.click();
  URL.revokeObjectURL(link.href);
}
