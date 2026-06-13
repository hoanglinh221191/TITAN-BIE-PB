# FMM/M2L Scope G Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Compare the FMM/M2L K1-K4 prototype with TITAN's production direct JCP operator on production-built real PQR surfaces.

**Architecture:** A dedicated offline C++ probe builds the production TABI surface, samples panels deterministically, computes production direct JCP output, and sends the same geometry and unknowns through a new prototype parity API plus existing Scope E/F gates. A PowerShell harness runs bounded real fixtures and emits the Scope G report.

**Tech Stack:** C++23, `TitanPBE::Tabi::build_surface`, production JCP direct fallback, FMM/M2L prototype, `TitanExecutor`, PowerShell, CSV/report outputs.

---

### Task 1: Direct-Parity API

**Files:**
- Modify: `src/tabi/tabi_fmm_m2l_prototype.h`
- Modify: `src/tabi/tabi_fmm_m2l_prototype.cpp`
- Test: `test/titan_fmm_m2l_scope_g_probe.cpp`

- [x] Write a failing C++ probe that requests `compare_jcp_direct_reference(...)` for production-built panels.
- [x] Build the probe and confirm failure for the missing type/function.
- [x] Add `JcpDirectParityMetrics` and compare combined K1+K2 / K3+K4 analytic references against supplied production vectors.
- [x] Rebuild and require finite parity with aggregate relative RMS `<= 1e-10`.

### Task 2: Real-Surface Probe

**Files:**
- Create: `test/titan_fmm_m2l_scope_g_probe.cpp`
- Create: `test/pbe_fmm_m2l_scope_g_probe.ps1`

- [x] Parse PQR coordinates, charge, and radius from the final five numeric columns.
- [x] Build a triangle surface with production `build_surface` and deterministic panel sampling.
- [x] Run production direct JCP fallback with zero LHS, then Scope E/F on the converted panels.
- [x] Export one CSV row with parity, K1, derivative, timing, work-count, and gate columns.
- [x] Run a small generated-surface smoke and require a passing row.

### Task 3: Real Fixture Gate And Documentation

**Files:**
- Create: `TITAN_FMM_M2L_SCOPE_G_REPORT.txt`
- Modify: `TITAN_FMM_M2L_ROADMAP.md`
- Modify: `TITAN_PBE_TABI_ROADMAP.md`
- Modify: `.gitignore`

- [x] Run the bounded `research-medium` probe and record the result.
- [x] Run `research-large` only if the medium row passes within the configured timeout.
- [x] Mark Scope G implemented only if production parity and M2L accuracy gates pass.
- [x] Run native build, Scope A smoke, Scope G smoke, and `git diff --check`.
