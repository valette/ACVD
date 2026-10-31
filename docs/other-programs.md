# Other Programs & Utilities in ACVD Repository

The ACVD repository contains auxiliary executables beyond the main remeshing tools. These programs provide supplementary functionality for volume processing, label image segmentation, mesh conversion, analysis, and cleaning.

---

## 1. Volume Processing Programs

### `VolumeAnalysis`
Generates surface meshes from 3D label volume segmentations. Reads MetaImage `.mhd`/`.mha`, NIFTI `.nii`, and NRRD `.nrrd`/`.nhdr` files (via `vtkRobustImageReader`), and extracts surfaces for each label using marching cubes, then simplifies them using ACVD or quadric methods.

**Usage:**
```bash
./build/bin/VolumeAnalysis file.mhd [options]
# NRRD and NIFTI volumes are also supported:
./build/bin/VolumeAnalysis file.nrrd [options]
./build/bin/VolumeAnalysis file.nii [options]
```

**Options:**
- `-n <number>`: Desired number of simplified mesh vertices (default: `32000`).
- `-v <number>`: Maximum number of voxels before image resampling (default: `100000000`).
- `-r <number>`: Maximum subsampling factor (default: `4`).
- `-j <number>`: Number of threads (default: auto).
- `-s <0|1>`: Simplification method (`0` = ACVD, `1` = Quadric).
- `-sm <number>`: Number of smoothing iterations (default: `0`).
- `-g <number>`: Gradation parameter for ACVD remeshing.
- `-o <directory>`: Output directory.
- `-f <format>`: Output mesh format (`vtk`, `ply`, `stl`, `obj`).
- `-m <0|1>`: Force 2-manifold output (default: `0`).
- `-a <0|1>`: Anisotropic coarsening (`1` = enabled, `0` = disabled).
- `-c <0|1>`: Retain only the largest connected component (`1` = enabled).
- `-t <threshold>`: Label threshold value.
- `-x <xmlfile>`: XML color palette file for labels.

---

## 2. Mesh Analysis & Manipulation Programs

### `ManifoldSimplification`
Performs manifold mesh simplification using the `vtkManifoldSimplification` engine.

**Usage:**
```bash
./build/bin/ManifoldSimplification input_mesh.ply n_vertices
```

**Functionality:**
- Reduces mesh complexity while preserving 2-manifold topological guarantees.
- Uses `vtkSurfaceBase` memory cleanup (`CleanMemory()`).

---

## 3. Auxiliary Examples & Tools (`Common/Examples/`)

### Mesh Format Converters
- `mesh2ply.cxx`: Convert meshes to PLY format.
- `mesh2vtk.cxx`: Convert meshes to VTK format.
- `mesh2stl.cxx`: Convert meshes to STL format.
- `mesh2obj.cxx`: Convert meshes to Wavefront OBJ format.
- `mesh2off.cxx`: Convert meshes to OFF format.
- `stl2ply.cxx`: Direct conversion from STL to PLY.

### Mesh Analysis & Inspection
- `meshCurvature.cxx`: Compute discrete principal curvature properties.
- `meshComponents.cxx`: Compute and extract connected components.
- `meshDifference.cxx`: Compute geometric difference between two meshes.
- `CheckOrientation.cxx`: Verify triangle normal orientation.
- `sampleMesh.cxx`: Uniform point sampling on surface meshes.

### Mesh Processing
- `subdivideMesh.cxx`: Apply Loop/Butterfly subdivision to meshes.
- `CleanMesh.cxx`: Clean topology, merge duplicate vertices, purge degenerate cells.
- `clipMesh.cxx`: Clip meshes using planes or bounding boxes.
- `viewer.cxx` / `viewer2.cxx`: Interactive VTK mesh visualization tools.
- `icp.cxx`: Iterative Closest Point 3D surface registration.

### Volume Processing Utilities (`VolumeProcessing/`)
- `VolumeSlice.cxx`: Extract 2D slices from label volumes.
- `VolumeCrop.cxx`: Crop 3D volume bounding boxes.
- `VolumeMedian.cxx`: Apply 3D median filtering to label volumes.
- `VolumeSubsample.cxx`: Subsample volume grids.
- `VolumeCleanLabels.cxx`: Clean isolated voxels and noisy label boundaries.
- `VolumeOOCSlice.cxx`: Out-of-core volume slice extraction.

---

## 4. Build Instructions

All auxiliary programs are built automatically with the standard CMake build process:
```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

Executables are written to `build/bin/` (or `bin/`).

---

## 5. Related Documentation

- [Root AGENTS.md](../AGENTS.md): Master repository guide & architectural context.
- [ACVD CLI Reference](../Readme.md#4-command-line-parameters-reference): Core ACVD remeshing CLI tool options.
- [vtkSurface Documentation](vtkSurface.md): Core mesh data structure documentation.