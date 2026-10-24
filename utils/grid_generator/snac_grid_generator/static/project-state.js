export class SnapshotHistory {
  constructor(limit = 100) {
    this.limit = limit;
    this.entries = [];
    this.index = -1;
  }

  reset(snapshot) {
    this.entries = [snapshot];
    this.index = 0;
  }

  commit(snapshot) {
    const current = this.entries[this.index];
    if (current?.project === snapshot.project && current?.decompositionResult === snapshot.decompositionResult) {
      return false;
    }
    this.entries = this.entries.slice(0, this.index + 1);
    this.entries.push(snapshot);
    if (this.entries.length > this.limit) this.entries.shift();
    this.index = this.entries.length - 1;
    return true;
  }

  undo() {
    if (!this.canUndo) return null;
    this.index -= 1;
    return this.entries[this.index];
  }

  redo() {
    if (!this.canRedo) return null;
    this.index += 1;
    return this.entries[this.index];
  }

  get canUndo() {
    return this.index > 0;
  }

  get canRedo() {
    return this.index >= 0 && this.index < this.entries.length - 1;
  }
}

// Returns false when browser storage is unavailable or full.
export function saveStoredProject(key, project) {
  try {
    localStorage.setItem(key, JSON.stringify({ savedAt: new Date().toISOString(), project }));
    return true;
  } catch (_error) {
    return false;
  }
}

export function clearStoredProject(key) {
  try {
    localStorage.removeItem(key);
  } catch (_error) {
    // Browser storage can be unavailable.
  }
}

export function loadStoredProject(key) {
  try {
    const value = JSON.parse(localStorage.getItem(key));
    return value?.project ? { ...value, key } : null;
  } catch (_error) {
    return null;
  }
}

// Stored projects whose keys start with *prefix* or equal one of *extraKeys*,
// newest first.
export function listStoredProjects(prefix, extraKeys = []) {
  let keys = [];
  try {
    keys = Object.keys(localStorage).filter((key) => key.startsWith(prefix) || extraKeys.includes(key));
  } catch (_error) {
    return [];
  }
  return keys
    .map((key) => loadStoredProject(key))
    .filter(Boolean)
    .sort((a, b) => String(b.savedAt).localeCompare(String(a.savedAt)));
}
