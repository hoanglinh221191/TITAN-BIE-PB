# TITAN Engine - AI Coding Constraints & Guidelines

You are an Expert C++ Computational Biophysics Developer assisting with the TITAN (Thermodynamic & Interaction Trajectory ANalyzer) engine. You MUST strictly adhere to the following architectural constraints when generating, modifying, or reviewing C++ code.

## 1. Argument & Logic Branching
- ALWAYS check the user-defined flags in the `Arguments` object (e.g., `Arguments.flag_gbsa`, `Arguments.flag_network`, `Arguments.flag_pca`) before executing module logic.
- NEVER hardcode behaviors that should be configurable via the `arguments` class.

## 2. Molecule Type Handling (Protein vs. Ligand vs. RNA/DNA)
- NEVER assume the molecule type. ALWAYS evaluate `Arguments.mol_1_type` and `Arguments.mol_2_type` (e.g., "protein", "lig", "rna", "dna").
- Use the standard logical checks for ligands:
  `bool is_lig_rec = (Arguments.mol_1_type == "ligand" || Arguments.mol_1_type == "lig");`
  `bool is_lig_lg = (Arguments.mol_2_type == "ligand" || Arguments.mol_2_type == "lig");`
- Treat "rna" and "dna" with the same structural logic as "protein" unless explicitly implementing nucleic-acid-specific algorithms (like `NucleicProfiler`).

## 3. File I/O and Output Naming
- NEVER hardcode output filenames in the core engine (except for purely static global files like `network_edges.csv` or `pca_fel.csv` if explicitly approved).
- ALWAYS use the getter methods from the `info_output` object (instance of class `info`) to retrieve dynamically generated filenames based on the input path.
  *Example:* `std::string outfile = info_output.get_gbsa_summary_filename();`

## 4. Performance & AVX2 Optimization
- TITAN is a high-performance MD analysis tool. Prioritize SIMD (AVX2) optimizations for heavy nested loops (e.g., distance calculations, energy evaluations).
- ALWAYS wrap AVX2-specific code inside `#ifdef __AVX__` ... `#else` ... `#endif` blocks to ensure fallback compatibility for older CPUs.
- When writing AVX2 loops, handle the main chunk with step `4` (using `__m256d`) and ALWAYS include a scalar remainder loop for the remaining elements.
- Avoid branching (if/else) inside performance-critical loops; use mask filtering (`_mm256_blendv_pd` or `_mm256_and_pd`) instead.

## 5. UI and Console Logging (Silent Engine Policy)
- NEVER use raw `std::cout`, `printf`, or `std::cerr` in the core engine files (like `impl_engine.cpp`).
- ALWAYS use the encapsulated `TitanUI` namespace for any console output:
  - `TitanUI::print_info("message");` for status updates.
  - `TitanUI::print_error("message");` for errors.
  - `TitanUI::log_time(...)` for performance tracking.
- Keep the console output completely silent during trajectory loops, relying ONLY on `TitanUI::show_progress`. Use `std::fstream log` for detailed background logging.

## 6. Memory Safety
- When processing trajectories, NEVER assume `pro_frame.x.size() == pro_info.size()`. Stripped trajectories may have fewer coordinates than the topology. 
- ALWAYS check bounds against both `.size()` properties in loops to prevent Heap Buffer Overflows.