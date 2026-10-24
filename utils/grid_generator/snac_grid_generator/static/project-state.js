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

export function saveStoredProject(key, project) {
  try {
    localStorage.setItem(key, JSON.stringify({ savedAt: new Date().toISOString(), project }));
  } catch (_error) {
    // Browser storage can be unavailable or full.
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
    return value?.project ? value : null;
  } catch (_error) {
    return null;
  }
}
