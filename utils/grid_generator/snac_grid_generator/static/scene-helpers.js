import * as THREE from "three";

const axisLabelMaterials = new Map();

export function clearDisposableGroup(group, disposeGeometries = true) {
  const geometries = new Set();
  const materials = new Set();
  for (const child of [...group.children]) {
    child.traverse((object) => {
      if (object.userData.sharedMaterial) return;
      if (disposeGeometries && object.geometry) geometries.add(object.geometry);
      if (!object.material) return;
      for (const material of Array.isArray(object.material) ? object.material : [object.material]) materials.add(material);
    });
    group.remove(child);
  }
  for (const geometry of geometries) geometry.dispose();
  for (const material of materials) material.dispose();
}

export function partitionFaceIndices(n, partitions) {
  if (partitions <= 1) return [];
  const base = Math.floor(n / partitions);
  const remainder = n % partitions;
  return Array.from({ length: partitions - 1 }, (_, rank) => {
    const count = rank + 1;
    return count * base + Math.min(count, remainder);
  });
}

export function planeOutlinePoints(block, axisIndex, value) {
  const [xmin, ymin, zmin] = block.lmin;
  const [xmax, ymax, zmax] = block.lmax;
  if (axisIndex === 0) {
    return [
      new THREE.Vector3(value, ymin, zmin),
      new THREE.Vector3(value, ymax, zmin),
      new THREE.Vector3(value, ymax, zmax),
      new THREE.Vector3(value, ymin, zmax),
      new THREE.Vector3(value, ymin, zmin),
    ];
  }
  if (axisIndex === 1) {
    return [
      new THREE.Vector3(xmin, value, zmin),
      new THREE.Vector3(xmax, value, zmin),
      new THREE.Vector3(xmax, value, zmax),
      new THREE.Vector3(xmin, value, zmax),
      new THREE.Vector3(xmin, value, zmin),
    ];
  }
  return [
    new THREE.Vector3(xmin, ymin, value),
    new THREE.Vector3(xmax, ymin, value),
    new THREE.Vector3(xmax, ymax, value),
    new THREE.Vector3(xmin, ymax, value),
    new THREE.Vector3(xmin, ymin, value),
  ];
}

export function axisLabel(label, color) {
  const key = `${label}-${color}`;
  if (axisLabelMaterials.has(key)) {
    const sprite = new THREE.Sprite(axisLabelMaterials.get(key));
    sprite.userData.sharedMaterial = true;
    return sprite;
  }
  const canvas = document.createElement("canvas");
  canvas.width = 64;
  canvas.height = 64;
  const context = canvas.getContext("2d");
  context.clearRect(0, 0, canvas.width, canvas.height);
  context.fillStyle = `#${color.toString(16).padStart(6, "0")}`;
  context.font = "700 42px Inter, Arial, sans-serif";
  context.textAlign = "center";
  context.textBaseline = "middle";
  context.fillText(label, 32, 34);
  const texture = new THREE.CanvasTexture(canvas);
  texture.colorSpace = THREE.SRGBColorSpace;
  const material = new THREE.SpriteMaterial({ map: texture, transparent: true, depthTest: false });
  axisLabelMaterials.set(key, material);
  const sprite = new THREE.Sprite(material);
  sprite.userData.sharedMaterial = true;
  return sprite;
}

export function projectBox(blocks) {
  const box = new THREE.Box3();
  for (const block of blocks) {
    box.expandByPoint(new THREE.Vector3(...block.lmin));
    box.expandByPoint(new THREE.Vector3(...block.lmax));
  }
  if (box.isEmpty()) {
    box.expandByPoint(new THREE.Vector3(0, 0, 0));
    box.expandByPoint(new THREE.Vector3(1, 1, 1));
  }
  return box;
}
