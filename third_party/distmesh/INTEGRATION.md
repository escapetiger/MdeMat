# DistMesh in MdeMat

This directory contains an unmodified upstream v1.2 snapshot. The exact commit,
original file list and SHA-256 identities are in [provenance.json](provenance.json).
The official source is [popersson/DistMesh](https://github.com/popersson/DistMesh).
This version is pure MATLAB and uses the included MIT license.

`setupMdeMat` adds only the normal toolbox namespace root. The callable
`toolbox/examples/geometry/ex04_nonuniform_trimesh.m` resolves this snapshot in
both a checkout and an installed toolbox, temporarily adds only `src`, and
restores the caller path and RNG. It owns the figure used by the upstream mesh
progress display; `MakePlots=false` uses a hidden temporary figure and removes
it before returning. Upstream meshing and convergence rules are unchanged.

From a loaded MdeMat toolbox:

```matlab
folder = fullfile(mdemat.tools.environment.toolboxRoot(), 'examples', 'geometry');
addpath(folder);
result = ex04_nonuniform_trimesh(ElementSize=0.25, MakePlots=false, Seed=0);
mesh = result.mesh;
```

To use the complete original API directly in a checkout, explicitly add
`third_party/distmesh/src` to the MATLAB path. Upstream callbacks use points as
**N-by-d**, which differs from MdeMat field callbacks' usual d-by-N convention.
The upstream generator displays mesh progress and consumes the active RNG.

The installable toolbox includes every `src/*.m` function and license/provenance
records. Upstream tests, examples, startup, CI and documentation are retained in
the repository snapshot but are not packaged or added to the default path.
MdeMat tests cover its real circular-mesh integration; they do not certify every
upstream geometry or guarantee convergence for arbitrary distance/size inputs.
