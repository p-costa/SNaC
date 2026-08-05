import { defaultBlock } from "./block-model.js";

export function builtInAxisPresets() {
  return [
    {
      id: "builtin-uniform-32",
      name: "Uniform 32",
      ng: 32,
      coordinateSpace: "native",
      axis: {
        kind: "simple_ratio",
        profile: "geometric",
        controls: ["n", "ratio"],
        ratio: 1,
        side: "end",
      },
    },
    {
      id: "builtin-wall-normal-64",
      name: "Wall-normal 64",
      ng: 64,
      coordinateSpace: "native",
      axis: {
        kind: "simple_ratio",
        profile: "tanh",
        controls: ["n", "ratio"],
        ratio: 8,
        side: "end",
      },
    },
    {
      id: "builtin-channel-three-region",
      name: "Channel three-region",
      ng: 96,
      coordinateSpace: "native",
      axis: {
        kind: "multi",
        profile: "geometric",
        segments: [
          { length: 20, cells: 30, ratio: 4, continuous: false },
          { length: 60, cells: 40, ratio: 1, continuous: false },
          { length: 20, cells: 30, ratio: 0.25, continuous: false },
        ],
      },
    },
  ];
}

export function builtInProjectTemplates() {
  return [
    template("builtin-single-block", "Single block", singleBlockProject()),
    template("builtin-two-block-channel", "Two-block channel", twoBlockChannelProject()),
    template("builtin-periodic-channel", "Periodic channel", periodicChannelProject()),
  ];
}

function template(id, name, project) {
  return { id, name, project };
}

function projectDefaults(name, blocks) {
  return {
    schemaVersion: 2,
    name,
    nscal: 0,
    periodicAxes: [false, false, false],
    decomposition: {
      targetRanks: 0,
      mode: "auto",
      axes: [true, true, true],
      minLocalCells: 4,
      maxLocalAspect: 0,
    },
    inferConnectivity: true,
    writeExternalGrid: true,
    externalGridSource: "grid",
    blocks,
  };
}

function singleBlockProject() {
  const block = defaultBlock(1);
  block.name = "domain";
  block.ng = [32, 32, 32];
  block.lmax = [1, 1, 1];
  return projectDefaults("single-block", [block]);
}

function twoBlockChannelProject() {
  const left = defaultBlock(1);
  left.name = "left";
  left.ng = [48, 64, 16];
  left.lmax = [1, 1, 0.25];
  left.axes.y = channelAxis();
  const right = JSON.parse(JSON.stringify(left));
  right.id = 2;
  right.name = "right";
  right.lmin[0] = 1;
  right.lmax[0] = 2;
  return projectDefaults("two-block-channel", [left, right]);
}

function periodicChannelProject() {
  const block = defaultBlock(1);
  block.name = "periodic-domain";
  block.ng = [64, 64, 16];
  block.lmax = [2, 1, 0.25];
  block.axes.y = channelAxis();
  const project = projectDefaults("periodic-channel", [block]);
  project.periodicAxes = [true, false, false];
  return project;
}

function channelAxis() {
  return {
    kind: "multi",
    profile: "geometric",
    segments: [
      { length: 20, cells: 30, ratio: 4, continuous: false },
      { length: 60, cells: 40, ratio: 1, continuous: false },
      { length: 20, cells: 30, ratio: 0.25, continuous: false },
    ],
  };
}
