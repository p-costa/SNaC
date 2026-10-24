// Camera framing, orthographic axis views, and renderer sizing.

import * as THREE from "three";
import { AXES } from "./block-model.js";
import { setStatus } from "./core.js";
import { projectBox } from "./scene-helpers.js";
import { els, selectedVisibleBlocks, ui, visibleBlocks } from "./state.js";

export function resizeRenderer() {
  const rect = els.scene.parentElement.getBoundingClientRect();
  const width = Math.max(1, Math.floor(rect.width));
  const height = Math.max(1, Math.floor(rect.height));
  ui.renderer.setSize(width, height, false);
  if (ui.camera.isOrthographicCamera) {
    updateOrthographicProjection(width / height);
  } else {
    ui.camera.aspect = width / height;
    ui.camera.updateProjectionMatrix();
  }
}

export function fitView() {
  activateCamera(ui.perspectiveCamera, "3d");
  ui.controls.enableRotate = true;
  framePerspective(
    visibleBlocks(),
    new THREE.Vector3(1.6, -2.2, 1.45),
  );
}

export function focusSelection() {
  const blocks = selectedVisibleBlocks();
  if (!frameBlocks(blocks)) {
    setStatus("Select a visible block to fit");
    return;
  }
  setStatus(`Fitted ${blocks.length} selected block${blocks.length === 1 ? "" : "s"}`);
}

// Frame *blocks* while keeping the current view direction and projection.
export function frameBlocks(blocks) {
  if (!blocks.length) return false;
  if (ui.camera.isOrthographicCamera) {
    if (axisViewDirection(els.scene.dataset.view)) {
      setAxisView(els.scene.dataset.view, blocks);
    } else {
      frameOrthographic(blocks);
    }
  } else {
    framePerspective(blocks, ui.camera.position.clone().sub(ui.controls.target));
  }
  return true;
}

// After geometry changes, re-frame the visible blocks when part of them left
// the view or when they shrank to a small corner of it.
export function keepModelInView() {
  const blocks = visibleBlocks();
  if (!blocks.length) return;
  const box = projectBox(blocks);
  ui.camera.updateMatrixWorld();
  const corners = [];
  for (const x of [box.min.x, box.max.x]) {
    for (const y of [box.min.y, box.max.y]) {
      for (const z of [box.min.z, box.max.z]) corners.push(new THREE.Vector3(x, y, z).project(ui.camera));
    }
  }
  const outside = corners.some((point) => Math.abs(point.x) > 1 || Math.abs(point.y) > 1 || Math.abs(point.z) > 1);
  const spanX = Math.max(...corners.map((point) => point.x)) - Math.min(...corners.map((point) => point.x));
  const spanY = Math.max(...corners.map((point) => point.y)) - Math.min(...corners.map((point) => point.y));
  const tiny = spanX < 0.25 && spanY < 0.25; // under an eighth of the view each way
  if (outside || tiny) frameBlocks(blocks);
}

export function framePerspective(blocks, direction) {
  const box = projectBox(blocks);
  const center = box.getCenter(new THREE.Vector3());
  const size = box.getSize(new THREE.Vector3());
  const radius = Math.max(size.x, size.y, size.z, 1e-9);
  const viewDirection = direction.lengthSq() > 0
    ? direction.clone().normalize()
    : new THREE.Vector3(1.6, -2.2, 1.45).normalize();
  ui.camera.position.copy(center).addScaledVector(viewDirection, radius * 3.1);
  ui.controls.target.copy(center);
  ui.controls.update();
  ui.gridHelper.position.copy(center);
  ui.gridHelper.scale.setScalar(Math.max(1, radius / 5));
  els.scene.dataset.framedBlocks = blocks.map((block) => block.id).join(",");
}

export function setAxisView(view, blocks = visibleBlocks()) {
  const direction = axisViewDirection(view);
  if (!direction) return;
  const axisIndex = AXES.indexOf(view.slice(1));

  const box = projectBox(blocks);
  const center = box.getCenter(new THREE.Vector3());
  const size = box.getSize(new THREE.Vector3());
  const up = axisIndex === 2
    ? new THREE.Vector3(0, 1, 0)
    : new THREE.Vector3(0, 0, 1);
  const projectedWidth = axisIndex === 0 ? size.y : size.x;
  const projectedHeight = axisIndex === 2 ? size.y : size.z;
  const radius = Math.max(size.x, size.y, size.z, 1);
  const distance = radius * 2.5;

  ui.orthographicViewSize = {
    width: Math.max(projectedWidth, 1e-9),
    height: Math.max(projectedHeight, 1e-9),
  };
  ui.orthographicCamera.position.copy(center).addScaledVector(direction, distance);
  ui.orthographicCamera.up.copy(up);
  ui.orthographicCamera.near = 0.01;
  ui.orthographicCamera.far = Math.max(10000, distance * 4);
  ui.orthographicCamera.zoom = 1;
  activateCamera(ui.orthographicCamera, view);
  ui.controls.enableRotate = true;
  ui.controls.target.copy(center);
  updateOrthographicProjection(ui.renderer.domElement.clientWidth / ui.renderer.domElement.clientHeight);
  ui.controls.update();
  ui.gridHelper.position.copy(center);
  ui.gridHelper.scale.setScalar(Math.max(1, radius / 5));
  els.scene.dataset.framedBlocks = blocks.map((block) => block.id).join(",");
}

export function frameOrthographic(blocks) {
  const box = projectBox(blocks);
  const center = box.getCenter(new THREE.Vector3());
  const size = box.getSize(new THREE.Vector3());
  const offset = ui.camera.position.clone().sub(ui.controls.target);
  const direction = offset.lengthSq() > 0
    ? offset.normalize()
    : new THREE.Vector3(1.6, -2.2, 1.45).normalize();
  const radius = Math.max(size.x, size.y, size.z, 1);

  ui.orthographicCamera.position.copy(center).addScaledVector(direction, radius * 2.5);
  ui.orthographicCamera.near = 0.01;
  ui.orthographicCamera.far = Math.max(10000, radius * 10);
  ui.orthographicCamera.zoom = 1;
  ui.controls.target.copy(center);
  ui.controls.update();

  const right = new THREE.Vector3(1, 0, 0).applyQuaternion(ui.camera.quaternion);
  const up = new THREE.Vector3(0, 1, 0).applyQuaternion(ui.camera.quaternion);
  ui.orthographicViewSize = {
    width: Math.max(projectedExtent(size, right), 1e-9),
    height: Math.max(projectedExtent(size, up), 1e-9),
  };
  updateOrthographicProjection(ui.renderer.domElement.clientWidth / ui.renderer.domElement.clientHeight);
  ui.gridHelper.position.copy(center);
  ui.gridHelper.scale.setScalar(Math.max(1, radius / 5));
  els.scene.dataset.framedBlocks = blocks.map((block) => block.id).join(",");
}

export function projectedExtent(size, direction) {
  return size.x * Math.abs(direction.x)
    + size.y * Math.abs(direction.y)
    + size.z * Math.abs(direction.z);
}

export function axisViewDirection(view) {
  const axisIndex = AXES.indexOf(view?.slice(1));
  if (axisIndex < 0) return null;
  const sign = view.startsWith("-") ? -1 : 1;
  return new THREE.Vector3().setComponent(axisIndex, sign);
}

export function handleCameraChange() {
  if (!ui.camera?.isOrthographicCamera) return;
  const expectedDirection = axisViewDirection(els.scene.dataset.view);
  if (!expectedDirection) return;
  const direction = ui.camera.position.clone().sub(ui.controls.target);
  if (direction.lengthSq() === 0) return;
  if (direction.normalize().dot(expectedDirection) < 1.0 - 1e-10) {
    setViewState("free-orthographic");
  }
}

export function activateCamera(nextCamera, view) {
  ui.camera = nextCamera;
  ui.controls.object = ui.camera;
  ui.transformControls.camera = ui.camera;
  if (ui.camera.isPerspectiveCamera) {
    ui.camera.aspect = ui.renderer.domElement.clientWidth / ui.renderer.domElement.clientHeight;
    ui.camera.updateProjectionMatrix();
  }
  els.scene.dataset.projection = ui.camera.isOrthographicCamera ? "orthographic" : "perspective";
  setViewState(view);
}

export function setViewState(view) {
  els.scene.dataset.view = view;
  for (const button of document.querySelectorAll("[data-axis-view]")) {
    const active = button.dataset.axisView === view;
    button.classList.toggle("active", active);
    button.setAttribute("aria-pressed", String(active));
  }
}

export function updateOrthographicProjection(aspect) {
  if (!ui.orthographicViewSize) return;
  const safeAspect = Math.max(aspect, 1e-9);
  const halfHeight = 0.56 * Math.max(
    ui.orthographicViewSize.height,
    ui.orthographicViewSize.width / safeAspect,
  );
  const halfWidth = halfHeight * safeAspect;
  ui.orthographicCamera.left = -halfWidth;
  ui.orthographicCamera.right = halfWidth;
  ui.orthographicCamera.top = halfHeight;
  ui.orthographicCamera.bottom = -halfHeight;
  ui.orthographicCamera.updateProjectionMatrix();
}
