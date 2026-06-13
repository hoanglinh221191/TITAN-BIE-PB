# FMM/M2L Scope I Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the Scope H finite-difference derivative bridge with a hidden offline analytical derivative operator and gate it against production direct JCP.

**Architecture:** Scope I stays inside `src/tabi/tabi_fmm_m2l_prototype.{h,cpp}` and test/probe files. It adds a new offline analytical result path that computes K1/K2/K3/K4 in fewer tree passes by evaluating analytical derivative kernels over the same near/far partition. It remains research-only and never wires into GMRES or production CLI.

**Tech Stack:** C++23, TITAN TABI prototype tree/worklist, `TitanExecutor`, PowerShell probe harnesses, CSV/report outputs.

---

### Task 1: Scope I Boundary

**Files:**
- Modify: `TITAN_FMM_M2L_ROADMAP.md`
- Modify: `TITAN_PBE_TABI_ROADMAP.md`
- Modify: `.gitignore`

- [ ] Add Scope I roadmap section after Scope H with these boundaries:
  - hidden offline analytical derivative operator only;
  - no GMRES wiring;
  - no production backend or CLI changes;
  - pass gate requires accuracy and speedup over production direct JCP.
- [ ] Whitelist:
  - `TITAN_FMM_M2L_SCOPE_I_REPORT.txt`
  - `test/pbe_fmm_m2l_scope_i_probe.ps1`
  - `test/titan_fmm_m2l_scope_i_probe.cpp`

### Task 2: RED Probe

**Files:**
- Create: `test/titan_fmm_m2l_scope_i_probe.cpp`
- Create: `test/pbe_fmm_m2l_scope_i_probe.ps1`

- [ ] Add a C++ probe that calls `evaluate_jcp_offline_operator_analytic(...)`.
- [ ] Verify it fails to compile because `JcpAnalyticOperatorResult` and `evaluate_jcp_offline_operator_analytic(...)` do not exist.

### Task 3: Analytical Operator

**Files:**
- Modify: `src/tabi/tabi_fmm_m2l_prototype.h`
- Modify: `src/tabi/tabi_fmm_m2l_prototype.cpp`

- [ ] Add `JcpAnalyticOperatorResult` with `phi`, `phin`, near/far counts, setup/apply timing, finite status, and fallback counters.
- [ ] Build a single panel-source tree and target tree for `PrototypePanel` geometry.
- [ ] For near leaf pairs, compute direct analytical K1/K2/K3/K4 formulas exactly per panel pair.
- [ ] For accepted far leaf pairs, compute cell-moment analytical approximations for K1/K2/K3/K4 using source cell scalar moments for `phi` and `phin`.
- [ ] Reuse `TitanExecutor` over target leaves; use worker-local vectors and deterministic reduce if multiple tasks can touch a target.

### Task 4: Gate And Reports

**Files:**
- Modify: `test/pbe_fmm_m2l_scope_i_probe.ps1`
- Create: `TITAN_FMM_M2L_SCOPE_I_REPORT.txt`

- [ ] Run generated-surface smoke.
- [ ] Run `research-medium`; only run `research-large` if medium passes.
- [ ] Gate:
  - parity aggregate relative RMS `<= 1e-10`;
  - analytical aggregate relative RMS `<= 1e-1`;
  - analytical max of `phi/phin` relative RMS `<= 2e-1`;
  - far-pair work nonzero;
  - analytical apply seconds faster than production direct, preferably `speedup > 1.2`.

### Task 5: Verification

**Commands:**

```powershell
mingw32-make PROFILE=native TARGET=titan.exe all -j1
powershell -NoProfile -ExecutionPolicy Bypass -File test\fmm_m2l_scope_a_smoke.ps1
powershell -NoProfile -ExecutionPolicy Bypass -File test\pbe_fmm_m2l_scope_h_probe.ps1 -Quick -Repeats 1
powershell -NoProfile -ExecutionPolicy Bypass -File test\pbe_fmm_m2l_scope_i_probe.ps1 -Repeats 5
git diff --check
```

Expected final decision:
- If Scope I passes speed and accuracy, plan a hidden single-matvec production adapter next.
- If Scope I is accurate but still `no_speedup`, stop FMM/M2L integration or redesign the far apply/cache architecture before any GMRES work.
