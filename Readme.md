[![CI](https://github.com/valette/ACVD/actions/workflows/ci.yml/badge.svg)](https://github.com/valette/ACVD/actions/workflows/ci.yml)

# ACVD - Approximated Centroidal Voronoi Diagrams

<p align="center">
  <img src="https://www.creatis.insa-lyon.fr/~valette/public/project/acvd/featured.jpg" alt="ACVD Featured Image">
</p>

**ACVD** is a C++ library and executable suite for adaptive coarsening, surface remeshing, and simplification of 3D triangular meshes built on top of **VTK (Visualization Toolkit)**. It implements discrete Centroidal Voronoi Diagram (CVD) algorithms based on publications by Sebastien Valette et al.

---

## 1. Documentation Links

- **[Agentic Coding Guide](./AGENTS.md)**: Master reference guide containing repository architecture, C++/VTK memory safety, $O(1)$ mesh topology mechanics (`vtkSurfaceBase`), CVD remeshing algorithms, metrics, and CMake build/test instructions.
- **[vtkSurface API Documentation](./docs/vtkSurface.md)**: Detailed API reference for `vtkSurfaceBase` and `vtkSurface` classes.
- **[Other Programs Documentation](./docs/other-programs.md)**: Guide for `VolumeAnalysis`, `ManifoldSimplification`, and mesh conversion tools.

---

## 2. Core Executables Overview

The build process generates command-line tools under `build/bin/` (or `bin/`):
- **`ACVD`**: Sequential isotropic discrete surface remesher.
- **`ACVDQ`**: Quadric-enhanced surface remesher (preserves sharp feature edges, corners, and boundaries).
- **`ACVDP`**: Multi-threaded (OpenMP) parallel version of `ACVD`.
- **`ACVDQP`**: Multi-threaded (OpenMP) quadric-enhanced parallel version of `ACVDQ`.
- **`AnisotropicRemeshingQ`**: Anisotropic metric surface remesher with quadric error bounds.
- **`AnisotropicRemeshingQP`**: Multi-threaded version of `AnisotropicRemeshingQ`.

---

## 3. Dependencies & Compilation

### Prerequisites
- **VTK** 9.0+ (`libvtk9-dev` on Debian/Ubuntu)
- **CMake** 3.11+
- **C++ Compiler** (GCC 7+ or Clang 6+) with C++11 support
- **OpenMP** (Optional, for multi-threaded executables `ACVDP`, `ACVDQP`, `AnisotropicRemeshingQP`)

### Build Instructions
```bash
git clone https://github.com/valette/ACVD.git
cd ACVD

# Out-of-source build configuration
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

Executables will be placed in `build/bin/` (or `bin/`).

---

## 4. Command Line Parameters Reference

### `ACVD` & `ACVDP` (Basic Isotropic Remeshing)

```bash
./build/bin/ACVD <input_mesh> <n_vertices> <curvature_weight> [options]
```

**Required Arguments:**
- `<input_mesh>`: Path to input 3D model (`.ply`, `.obj`, `.vtk`, `.off`, `.wrl`, `.stl`).
- `<n_vertices>`: Desired number of vertices in output mesh (e.g. `3000`).
- `<curvature_weight>`: Gamma parameter controlling curvature gradation (`0.0` = uniform density, `1.0` - `2.0` = higher density in high-curvature areas).

**Options:**
- `-s <ratio>`: Subsampling threshold ratio (default: `10`). Setting `-s 100` significantly improves output mesh quality.
- `-l <ratio>`: Edge splitting ratio (`-l 3` splits long edges exceeding 3 times the average edge length).
- `-m <0|1>`: Enforce 2-manifold output topology (`1` = enabled, `0` = disabled).
- `-d <0|1|2>`: Interactive visualization level (`0` = off, `1` = window on completion, `2` = step-by-step iteration display). Press key `'e'` in the display window to advance iterations.
- `-o <directory>`: Custom output directory.
- `-of <filename>`: Custom output filename (default: `simplification.ply`).
- `-b <0|1>`: Boundary vertex fixing (`1` = fixed, `0` = free).
- `-q <0|1|2|3>`: Quadric optimization level (default: `3`).
- `-np <n_threads>`: Set number of threads (vtkSMPTools) for parallel executables (`ACVDP`, `ACVDQP`).

---

### `ACVDQ` & `ACVDQP` (Quadric-Enhanced Remeshing)

```bash
./build/bin/ACVDQ <input_mesh> <n_vertices> <curvature_weight> [options]
```

**Options:** Includes all `ACVD` options plus:
- `-q <1|2|3>`: Number of quadric matrix eigenvalues used for vertex positioning (default: `3`).
- `-fv <file>`: Text file containing IDs of fixed vertices.
- `-ft <file>`: Text file containing IDs of fixed triangles.
- `-np <n_threads>`: Set number of threads for `ACVDQP`.
- `-ipr <n_items>`: Set the minimum number of items per parallel region for `ACVDQP` (default 512). The result does not depend on the number of threads.

---

## 5. Output Mesh Quality vs. Complexity Tuning

ACVD works by clustering input mesh vertices. Output quality depends on two key parameters:

1. **Subsampling Ratio (`-s ratio`)**: Controls the ratio of intermediate samples to output vertices. Defaults to `10`. Increasing to `-s 100` improves sampling accuracy and output triangle quality.
2. **Edge Splitting Ratio (`-l ratio`)**: Pre-splits input edges longer than ratio X average edge length. Use `-l 3` to eliminate long, thin triangles in input models.
3. **Manifold Enforcement (`-m 1`)**: Guarantees that the output mesh is topologically 2-manifold (no self-intersections or non-manifold edges).

---

## 6. Usage Examples

### Remeshing the Stanford Bunny to 3000 Vertices (Uniform)
```bash
wget https://github.com/alecjacobson/common-3d-test-models/raw/master/data/stanford-bunny.obj
./build/bin/ACVD stanford-bunny.obj 3000 0.0
```

### Curvature-Adapted High-Quality Remeshing
```bash
./build/bin/ACVD stanford-bunny.obj 3000 1.5 -s 100 -l 3 -m 1
```

### Feature-Preserving Quadric Remeshing (`ACVDQ`)
```bash
wget https://github.com/alecjacobson/common-3d-test-models/raw/master/data/fandisk.obj
./build/bin/ACVDQ fandisk.obj 3000 0.0 -s 100 -l 3 -m 1
```

### Anisotropic Metric Remeshing (`AnisotropicRemeshingQ`)
```bash
wget https://github.com/alecjacobson/common-3d-test-models/raw/master/data/horse.obj
./build/bin/AnisotropicRemeshingQ horse.obj 1000 1.5
```

### Parallel Quadric Remeshing (`ACVDQP`)
```bash
wget http://graphics.stanford.edu/data/3Dscanrep/xyzrgb/xyzrgb_statuette.ply.gz
gunzip xyzrgb_statuette.ply.gz
./build/bin/ACVDQP xyzrgb_statuette.ply 100000 1.5 -np 8 -s 50
```

---

## 7. Python Port

A Python wrapper and PyVista port of ACVD is available at [pyacvd](https://github.com/pyvista/pyacvd).

---

## 8. Publications & References

1. **S. Valette, J.-M. Chassery and R. Prost**, *Generic remeshing of 3D triangular meshes with metric-dependent discrete Voronoi Diagrams*, IEEE Transactions on Visualization and Computer Graphics, Volume 14, no. 2, pages 369-381, 2008.
2. **S. Valette and J.-M. Chassery**, *Approximated Centroidal Voronoi Diagrams for Uniform Polygonal Mesh Coarsening*, Computer Graphics Forum (Eurographics 2004 proceedings), Vol. 23, No. 3, September 2004, pp. 381-389.
3. **M. Audette, D. Rivière, M. Ewend, A. Enquobahrie, and S. Valette**, *Approach-guided controlled resolution brain meshing for FE-based interactive neurosurgery simulation*, Workshop on Mesh Processing in Medical Image Analysis, MICCAI 2011, pp. 176-186, 2011.
