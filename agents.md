# Agents.md — Normal→Displacement Converter (CPU, Windows 11/Linux)

**Audience:** Code-generating/implementing agents (Codex-like)

**Goal:** Build a **CPU‑only** tool that converts tangent‑space normal maps to **scalar displacement** via **Poisson integration**, with robust handling of UV islands, UDIMs, and multiple UV sets.

**Target OS/Toolchain:** Windows 11, Visual Studio 2022, CMake (vcpkg for deps). Linux (Ubuntu 24.04, GCC/Clang, CMake/Ninja).

**Scope:** CLI first; **no GPU**; clean internal API enabling a future GUI.

**Primary Output:** EXR displacement map(s). **Default:** single‑channel 32‑bit float EXR (channel name `height`). Optional RGBA with R=height, G/B/A=0.

**Key features:** Mesh+material import; normal map detection/assignment (per material); UDIM patterns; sidecar generation (chart masks, mirrored flags, optional metric); per‑island masked Poisson solve; deterministic multi‑threading option.

---

## 1) High‑Level Overview

**Problem.** Normal maps encode a *direction* **n = (x,y,z)** in tangent space; we want a *height field* **h(u,v)** whose gradients match the normal‑implied slopes. We:

1. **Decode normals → slopes** in UV space:
   ∂h/∂u=−x/max(z,ε), ∂h/∂v=−y/max(z,ε).
2. **(Optional) metric scaling** from texel units → millimeters.
3. **Compute divergence** of the slope field within each UV island.
4. **Solve Poisson** Δh=∇⋅g on each island with **Neumann (zero‑flux)** boundaries (outer and **inner boundaries/holes**) and a single **nullspace constraint** (anchor or mean‑zero).
5. **Export** displacement EXR with metadata.

**Mesh awareness.** We load the mesh to derive UV charts/islands, island mirror flags, UDIM tile mapping, material assignments, and (optionally) a per‑texel metric. Computation happens in image space, guided by these derived assets.

**Out of scope (v1):** GPU compute; vector displacement; DCC plugins; Vulkan/DX backends.

---

## 2) Tech Stack & Libraries

* **Language:** C++20
* **Build:** CMake (Visual Studio generator on Windows, Ninja on Linux); **vcpkg** pinned (triplet: `x64-windows` or `x64-linux`)
* **Core libs:**

  * **Eigen 3.4** — sparse/dense linear algebra (Poisson solve)
  * **oneTBB** — parallelism (per‑island / per‑UDIM)
  * **OpenImageIO (OIIO)** — image I/O (PNG/TIF/EXR; UDIM patterns)
  * **Assimp** — mesh/material import (OBJ/FBX/GLTF/…)
  * **CLI11** — CLI parsing
  * **spdlog** — logging
  * **nlohmann/json** — sidecars & config
  * *(Optional)* **amgcl** — algebraic multigrid preconditioner (faster Poisson)

**Agent note:** Control this via a CMake option `N2D_USE_AMGCL` (default **OFF**). Guard usage with `#if N2D_USE_AMGCL` and add `option(N2D_USE_AMGCL "Enable amgcl" OFF)` in the top-level CMakeLists.

All dependency versions are pinned in `vcpkg.json` and locked via `vcpkg-configuration.json`.

**vcpkg baseline pinning:** Use a real commit SHA for the builtin baseline. Either set `"builtin-baseline": "<SHA>"` in `vcpkg.json` or set the default-registry baseline in `vcpkg-configuration.json`—if both are present, they **must** be identical. The SHA must exist in the local `VCPKG_ROOT` clone used for the build. Do **not** use placeholders in commits/PRs.

---

## 3) Repository Layout

```
normal2disp/
  CMakeLists.txt
  vcpkg.json
  vcpkg-configuration.json
  include/n2d/
    core/*.hpp        # math, types, errors
    io/*.hpp          # mesh & image I/O facades
    uv/*.hpp          # charts, UDIM, masks, metrics
    solve/*.hpp       # laplacian assembly, poisson
    bake/*.hpp        # normals→slopes→div→solve pipeline
    util/*.hpp        # logging, threads, hashing
  src/
    core/*.cpp
    io/*.cpp
    uv/*.cpp
    solve/*.cpp
    bake/*.cpp
    util/*.cpp
  apps/
    n2d_cli.cpp       # main() for CLI tool
  tests/
    CMakeLists.txt
    unit/*.cpp        # GoogleTest tests
    data/             # meshes, textures, UDIM samples
  docs/
    Agents.md
    DESIGN.md         # API surfaces for future GUI
  out/                # build outputs (git-ignored)

```

**Note (local test data):** See **§8A Local Test Data & PR Rules** for detailed rules on using and generating test data.

---

## 8) Development Plan (Phases & Tests)

**P0 — Bootstrap**
Deliver: buildable skeleton, `n2d --help`, logging.
Tests: build in x64‑Release; smoke tests.

**P1 — Mesh & Material Inspection**
Deliver: materials, UV sets, UDIMs, normal guesses, Y‑convention guess, overlap diagnostics.
Tests: OBJ/FBX/GLTF ground truths; UDIM detection.

**P2 — Texture I/O & UDIM Handling**
Deliver: load 2/3/4‑ch normals; UDIM pattern expansion; sRGB warnings.
Tests: 8‑bit→float correctness; tile mapping; Y‑flip unit test.

**P3 — UV Chart Masks & Mirroring**
Deliver: island rasterization, mirrored flags, sidecars. **Islands with holes must be represented; treat inner boundaries (holes) as Neumann (zero‑flux) during divergence/solve.**
Tests: island counts; mirrored shell detection; coverage vs UV area (±1%). Add an edge‑case test for **donut islands (islands with holes)** to ensure correct handling of non‑simply‑connected domains.

**P4 — Normals→Slopes & Divergence**
Deliver: normalization policy, slope computation, staggered divergence with zero‑flux boundaries.
Tests: flat normal → div≈0; synthetic gradient → analytic divergence match.

**P5 — Poisson Solver**
Deliver: masked Laplacian, anchor/mean‑zero, IC0 preconditioner. **Inner-boundary (hole) treatment is Neumann zero‑flux; ensure masked operator topology handles non‑simply‑connected islands.**
Tests: procedural height → normals → recon RMS < 1e‑3 (512²); border flatness (no bevel); solver convergence. Include an additional donut island case here if not covered in P3.

**P6 — UDIM & Multithreading**
Deliver: per‑tile + per‑island parallelism; deterministic reductions.
Tests: two tiles with different island counts; `--threads 1` vs default identical hashes.

**P7 — CLI Polish & Sidecars**
Deliver: `n2d inspect --json`; bake flags incl. `--normalization`, `--max-slope`, `--replicate-udims`, `--debug-dumps`, `--deterministic`.
Tests: end‑to‑end sample asset; EXR metadata correct.

**P8 — Metric Map (Optional)**
Deliver: rasterized `(mm_per_texel_u,v)`; flags/IO present; pipeline can skip.
Tests: two regions with 2× density → amplitude equalized with metric.

**Additional Tests**

Use `testdata/` for integration/smoke coverage; generate minimal fixtures at runtime for golden/edge-case tests requiring analytic ground truth or specific conditions not represented by the real assets.

* **Channel variants:** 2‑ch (BC5) vs 3‑ch yield matching results within tolerance.
* **Compression sim:** quantize/perturb channels; compare **raw**, **xyz‑norm**, **xy‑recon** pipelines.
* **Max‑slope clamp:** near‑horizon normals remain stable; bounded energy.
* **Determinism:** `--deterministic` bit‑for‑bit across runs/cores; regular mode nearly identical.
* **Units:** metadata correctness (`metric`/`units`).

---

## 8A) Local Test Data & PR Rules

* Local assets live under `testdata/` and are expected to be present. Prefer them for integration and I/O smoke tests (e.g., `inspect` on a real mesh, Y-channel guess behavior).
* Some tests intentionally generate fixtures even when `testdata/` exists (analytic golden cases, BC5 2-channel normals, UDIM tile synthesis, mirrored/overlap UVs, donut islands, near-horizon normals). Keep generated fixtures tiny (≤16×16) and sparse.
* Discovery and selection: if environment variable `N2D_TESTDATA_DIR` is set, use that; otherwise use repository `testdata/`. Tests may still choose runtime-generated fixtures when the scenario requires ground truth or specific edge conditions.
* PR policy: do not add or modify binary assets under `testdata/`. Do not commit any generated files. Ensure `.gitignore` ignores `build/` and patterns like `*.exr`, `*.png`, `*.tif`, `*.bin`.
* Generated fixtures must be written to a per-test temporary directory (e.g., under the CTest working dir) and deleted at test end. Never write inside the repo tree.
* When available, use the head normal map in `testdata/` to exercise +Y/−Y decoding and normalization policies.

---

## 13) Phase‑wise "Done" Checks (Quick)

* **P0:** VS builds; `--help` works.
* **P1:** Materials & UV sets; UDIMs; overlaps.
* **P2:** Normal maps load (2/3/4‑ch); UDIM expansion; sRGB warnings.
* **P3:** Chart masks/flags exported; coverage sane; donut island handled.
* **P4:** Divergence centered; flat → near‑zero.
* **P5:** Poisson passes golden RMS; border remains flat; donut island solved correctly.
* **P6:** UDIM end‑to‑end; multi‑thread determinism.
* **P7:** CLI polish; EXR metadata; JSON schema validated.
* **P8:** Metric map optional path proven.

---

## 14) CI & Container (Docker) Setup

This project builds cleanly in small cloud/Docker runners. CI should bootstrap or prebake vcpkg and point CMake at the vcpkg toolchain.

**Triplets:** Use `x64-windows` on Windows runners, `x64-linux` on Linux runners (for smoke builds).

**Option A — Bootstrap at runtime (portable):**

1. Install base tools in the container: git, curl, unzip/zip, CMake, Ninja (recommended), gcc/g++ (or MSVC on Windows), pkg-config.
2. Clone vcpkg in a fixed path (e.g., `/opt/vcpkg` or `C:/src/vcpkg`).
3. Ensure the clone contains the pinned baseline SHA from `vcpkg.json`/`vcpkg-configuration.json`:

   * `git -C <VCPKG_ROOT> fetch origin <BASELINE_SHA> --depth=1` (if needed)
   * `git -C <VCPKG_ROOT> checkout <BASELINE_SHA>`
4. Bootstrap vcpkg: `<VCPKG_ROOT>/bootstrap-vcpkg.sh` (Linux) or `.bat` (Windows).
5. Set env for the build step:

   * `VCPKG_ROOT=<VCPKG_ROOT>`
   * `VCPKG_FEATURE_FLAGS=manifests,binarycaching`
   * Optional cache: `VCPKG_DEFAULT_BINARY_CACHE=/work/.vcpkg-cache` (mount as a volume between runs)
6. Configure with the toolchain:

   * Linux (Ninja): `cmake -S . -B build -G Ninja -DCMAKE_TOOLCHAIN_FILE=$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake -DVCPKG_TARGET_TRIPLET=x64-linux`
   * Windows (VS): `cmake -S . -B build -G "Visual Studio 17 2022" -A x64 -DCMAKE_TOOLCHAIN_FILE=%VCPKG_ROOT%/scripts/buildsystems/vcpkg.cmake -DVCPKG_TARGET_TRIPLET=x64-windows`
7. Build & test: `cmake --build build --config Release` then `ctest --test-dir build -C Release`.

**Option B — Prebake vcpkg into the image (fast):**

* Create a base image with vcpkg cloned, bootstrapped, and checked out at the baseline SHA. Install compilers, CMake, Ninja.
* Mount a persistent volume for `VCPKG_DEFAULT_BINARY_CACHE` to avoid rebuilding ports each run.
* CI jobs then run only the CMake configure/build/test commands from Option A step 6–7.

**CMakePresets (recommended):** Provide a `CMakePresets.json` with presets that bake in the toolchain path, triplet, and generator (e.g., `ci-ninja`, `vs-windows`). CI can call `cmake --preset ci-ninja` and `cmake --build --preset ci-ninja`.

**Build hygiene (CI‑enforced):**

* Treat warnings as errors: **MSVC** `/W4 /WX`; **GCC/Clang** `-Wall -Wextra -Wpedantic -Werror`.
* Provide a repository-wide **.clang-format** and/or **.editorconfig**; add a CI step to check formatting (e.g., `clang-format --dry-run --Werror`).
* For third-party headers, consider suppressing noisy diagnostics locally (e.g., push/pop pragmas) rather than globally relaxing warning levels.

**Common CI pitfalls:**

* “invalid builtin-baseline”: the baseline SHA in the manifest/config does not exist in the container’s vcpkg clone. Fetch/checkout that SHA before bootstrap.
* Running `vcpkg install` manually: unnecessary in manifest mode—CMake + toolchain will drive installation.
* Missing binary cache: enable `VCPKG_FEATURE_FLAGS=binarycaching` and set a persistent `VCPKG_DEFAULT_BINARY_CACHE` to speed up builds.

---

## 15) Glossary

* **Neumann (zero‑flux) boundary:** derivative normal to boundary is zero; implemented by reusing interior edge flux in divergence at the boundary.
* **Staggered‑flux divergence:** difference of face‑centered fluxes; more stable pairing with 5‑point Laplacian than naive 4‑neighbor differences.
* **BC5/ATI2N:** 2‑channel compressed normal format; requires reconstructing Z from X,Y.
* **IC0:** Incomplete Cholesky preconditioner; accelerates CG for SPD systems.
