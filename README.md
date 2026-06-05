# TITAN Engine 🚀

**TITAN** is a high-performance, HPC-grade molecular trajectory analysis engine written in C++23. It is designed to process massive Molecular Dynamics (MD) trajectories (PDB, GRO, XTC, TRR) with a focus on memory efficiency, streaming accumulation, and SIMD-accelerated calculations.

Unlike generic MD analysis tools that load full trajectories into RAM, TITAN employs a highly optimized chunk-streaming architecture. It performs cheap, SIMD-friendly work inside the trajectory loop and orchestrates heavy mathematical solves (like Eigenvalue decomposition for TICA or clustering for MSM) during post-processing.

## 🌟 Key Features

TITAN provides an end-to-end pipeline covering both structural thermodynamics and advanced kinetic network models:

### 1. Thermodynamics & Solvation
* **Interaction Energy & Force Analysis:** Highly optimized electrostatic and van der Waals (vdW) calculations.
* **GBSA / ALPB:** Implicit-solvent free-energy terms and decomposition.
* **GIST:** 3D grid-based mapping of water thermodynamics, solvent occupancy, and entropy.
* **Hydrogen Bonds & IFP:** Interaction Fingerprint (IFP) for tracking H-bonds, salt bridges, hydrophobic, $\pi-\pi$, and cation-$\pi$ contacts.

### 2. Kinetics & Dimensionality Reduction (The Markov Pipeline)
* **TICA (Time-lagged Independent Component Analysis):** Extracts slow collective coordinates. Features auto-scaling backends (`exact_dense`, `exact_dual`, `randomized`) for massive systems.
* **MSM (Markov State Models):** Builds discrete transition probability matrices (TPM), computes stationary distributions, and performs PCCA macrostate lumping.
* **DHAM (Dynamic Histogram Analysis Method):** Reweights MSM state populations using replica temperatures and GROMACS `.xvg` energy traces.
* **UMAP:** Projects TICA coordinates into a non-linear 2D/3D Free Energy Surface (FES) landscape.

### 3. Networks & Allostery
* **DCCM (Dynamic Cross-Correlation Matrix):** Network analysis from aligned $C\alpha$/backbone trajectories.
* **PCA / FES:** Principal Component Analysis coupled with automated basin extraction.
* **Allostery & LSP-MD:** Estimates nonlinear coupling and directed information flow (Transfer Entropy / Generalized Correlation).

### 4. Specialized Modules
* **LLPS (Liquid-Liquid Phase Separation):** Detects dense/dilute phases, percolation, and spatial clustering (optimized for MARTINI coarse-grained systems).
* **LLPS Paper Fractal Morphology:** With `-llps-paper -llps-paper-fractal=auto`, TITAN exports dense-phase bead/residue fractal descriptors for MARTINI protein/RNA condensates. The citeable sidecar files are `llps_paper_fractal_timeseries.csv`, `llps_paper_fractal_scales.csv`, `llps_paper_fractal_summary.dat`, and `llps_paper_fractal_readme.md`. These descriptors report morphology only; interpret them together with phase partition, convergence, contact, RNA compaction/topology, and thermodynamic evidence rather than as standalone LLPS/binodal/spinodal proof.
* **CAS:** Computational Alanine Scanning for hotspot identification.

## 🛠️ Build Instructions

TITAN is built for modern architectures, heavily utilizing `std::jthread` for scoped parallelism and AVX2 for SIMD vectorization.

### Requirements
* **C++ Standard:** C++23
* **Compiler:** GCC 13+ (Linux) or Clang (Windows)

### Windows (MinGW/Clang)
A batch script is provided for standard builds:
```bash
./compile.bat
(Uses mingw32-make.exe. Use the release target for LTO and aggressive optimization).

Linux (GCC)
Bash
make release
🚀 Quick Start & Examples
TITAN supports both direct CLI flags and argument file (-a arg.in) parsing.

1. Basic Energy Decomposition (GBSA)
Calculate binding free energy components from a trajectory using a reference structure:

Bash
titan -s complex.pdb -traj md.xtc -g 1\&2 -ff AMBER99SB -gbsa -igb 5
2. The Complete Kinetics Pipeline (TICA + UMAP + MSM)
Process multiple replicas directly from a text list, extract slow dynamics, and build a Markov model:

Bash
# Provide a list of .xtc files in traj_list.txt
titan -s reference.pdb -i traj_list.txt -tica --tica-lag=10 --tica-components=5 -umap -msm -msm-clusters 150 -msm-lag 10
3. Condensate Analysis (LLPS)
Analyze a MARTINI coarse-grained system for phase separation:

Bash
titan -s system.pdb -traj md.xtc -ff MARTINI -llps -llps-dbscan-eps 12 -llps-dbscan-minpts 8 -llps-fuzzy

4. LLPS Paper Fractal Morphology
Generate citeable dense-phase fractal morphology outputs for a MARTINI protein/RNA condensate:

Bash
titan -s system.pdb -traj md.xtc -ff MARTINI -llps-paper -llps-paper-fractal=auto -llps-dbscan-eps 12 -llps-dbscan-minpts 8
🧠 Architecture Overview
TITAN's runtime philosophy revolves around:

Chunk Streamer (chunk_streamer.cpp): Reads XDR/PDB chunks and partitions coordinates.

Engine Loop (impl_engine.cpp): The hot path. Routes per-frame work to modules.

Post-Processing (engine_postproc.cpp): Executes heavy solves (e.g., K-means, Eigen solvers) after the streaming loop finishes, avoiding memory bloat.
