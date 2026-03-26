# CRITICAL DIRECTIVES FOR AGENT
You are an expert C++23 systems programmer specializing in high-performance scientific computing, SIMD intrinsics, and lock-free concurrency.

Architecture/codebase map:
- Read `AGENT_ONBOARDING.md` first for TITAN workflow, file ownership, TICA/MSM/DHAM data flow, and HPC conventions.

## 1. FILE OPERATIONS & SEARCHING (NO PYTHON SCRIPTS)
- **NEVER** write ad-hoc Python scripts to do basic Search and Replace. 
- Always use fast, native Linux CLI tools (`sed`, `perl -pi -e`, `awk`, `find`, `rg`) directly in the terminal for refactoring the codebase.
- Example for replacement: `find . -type f \( -name "*.cpp" -o -name "*.h" \) -exec perl -pi -e 's/OLD/NEW/g' {} +`
- **RESPECT .gitignore:** Never attempt to read, scan, or modify files inside `build/`, `bin/`, or files ending in `.xtc`, `.pqr`, `.o`, `.trr`.

## 2. TERMINAL DISCIPLINE
- **SUPPRESS OUTPUT:** If a command produces more than 50 lines, pipe it or redirect it (e.g., `> output.log`). DO NOT flood the context window with raw build logs or massive file lists.
- **COMPILE FAST:** Always compile using multi-threading: `make -j$(nproc)` or `ninja`.

## 3. C++ CODING STANDARDS (SoA & AVX2)
- The project is using a strict **Structure of Arrays (SoA)** layout (`vector<double> x, y, z`).
- Write perfectly vectorizable loops using flat arrays and intrinsic functions (AVX2/SIMD).
-đảm bảo các câu lệnh in ra màn hình, log file theo chuẩn titanUI
1. The GCC 13 C++23 ALLOWLIST (Features you MUST use):

String inspection: Use .contains(...) for std::string and std::string_view instead of .find().

Modern Ranges: Use <ranges> (specifically std::views::zip and std::views::enumerate) to optimize and modernize loops.

Multidimensional Subscript: Use the C++23 operator[](size_t r, size_t c) when designing custom flat-matrix wrappers.

Error Handling: Use <expected> (std::expected<T, E>) for clean, exception-free execution flows.

2. The GCC 13 C++23 BLOCKLIST (Features you MUST NEVER use):

NO <mdspan>: It is not supported natively in GCC 13. Write a lightweight FlatMatrix struct utilizing a 1D std::vector and operator[](x, y) instead.

NO <print>: Do not use std::print or std::println. Stick to <iostream> or std::format.

NO C++20 Modules: Do not use import std;. Stick to traditional #include.

3. HPC & Memory Safety Constraints:

AVX2 Alignment: When modifying memory allocation, you must preserve the #ifdef _WIN32 logic for _aligned_malloc. For the #else (Linux) fallback, you MUST use std::aligned_alloc(Alignment, size), and CRITICALLY, you must bitwise round up the size to a multiple of the alignment: size_t size = (n * sizeof(T) + Alignment - 1) & ~(Alignment - 1);.

Flatten Matrices: Never use nested std::vector<std::vector<T>>.
##4.How to compile and run
# BUILD & COMPILATION RULES
1. **Always use the Batch script for compiling:** Whenever you need to compile, build, or test the project, YOU MUST execute `compile.bat` in the terminal.
2. **Strictly NO raw g++ commands:** DO NOT ever use or suggest direct `g++` commands (e.g., `g++ *.cpp`) in the terminal. This project is extremely heavy and raw compilation causes 16GB+ RAM overflow via `cc1plus`.
3. **Build Architecture:** The project uses a custom `Makefile` via MSYS2 (`mingw32-make.exe -j`) for incremental builds to save RAM and compilation time. 
4. **Modifying Build Steps:** If you need to add a new library, flag, or change compilation behavior, modify the `Makefile` directly. DO NOT modify `compile.bat` unless absolutely necessary.
5. **Release Build:** If the user explicitly asks to build an optimized release version, instruct them to run `C:\msys64\ucrt64\bin\mingw32-make.exe release` in the terminal.
##5.examples:
thư mục example/1 là ví dụ protein-ligand và example/3 là protein-ligand có nước
example/2 là ví dụ protein-protein
example/4 là ví dụ protein-rna
sau khi sửa header files (.h), nên chạy make clean trước khi make để tránh stale object files. Makefile hiện tại chỉ track thay đổi .cpp, không có dependency tracking tự động cho headers.
