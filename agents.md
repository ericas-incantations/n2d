# Agents.md — Normal→Displacement Converter (CPU, Windows 11)

**Audience:** Code-generating/implementing agents (Codex-like)

**Goal:** Build a **CPU‑only** tool that converts tangent‑space normal maps to **scalar displacement** via **Poisson integration**, with robust handling of UV islands, UDIMs, and multiple UV sets.

**Target OS/Toolchain:** Windows 11, Visual Studio 2022, CMake (vcpkg for deps).

**Scope:** CLI first; **no GPU**; clean internal API enabling a future GUI.

**Primary Output:** EXR displacement map(s). **Default:** single‑channel 32‑bit float EXR (channel name `height`). Optional RGBA with R=height, G/B/A=0.

**Key features:** Mesh+material import; normal map detection/assignment (per material); UDIM patterns; sidecar generation (chart masks, mirrored flags, optional metric); per‑island masked Poisson solve; deterministic multi‑threading option.

---

## 1) High‑Level Overview

**Problem.** Normal maps encode a *direction* **n = (x,y,z)** in tangent space; we want a *height field* **h(u,v)** whose gradients match the normal‑implied slopes. We:

1. **Decode normals → slopes** in UV space:
   $\partial h/\partial u = -x/\max(z,\varepsilon)$, $\partial h/\partial v = -y/\max(z,\varepsilon)$.
2. **(Optional) metric scaling** from texel units → millimeters.
3. **Compute divergence** of the slope field within each UV island.
4. **Solve Poisson** $\Delta h = \nabla\cdot g$ on each island with **Neumann (zero‑flux)** boundaries and a single **nullspace constraint** (anchor or mean‑zero).
5. **Export** displacement EXR with metadata.

**Mesh awareness.** We load the mesh to derive UV charts/islands, island mirror flags, UDIM tile mapping, material assignments, and (optionally) a per‑texel metric. Computation happens in image space, guided by these derived assets.

**Out of scope (v1):** GPU compute; vector displacement; DCC plugins; Vulkan/DX backends.

---

## 2) Tech Stack & Libraries

* **Language:** C++20
* **Build:** CMake (Visual Studio generator); **vcpkg** pinned (triplet: `x64-windows`)
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

**Note (local test data):** A repository-level `testdata/` directory is present with binary assets (e.g., a character FBX and a head normal map). Use these assets for integration/smoke tests where appropriate. For unit/golden or edge-case scenarios that the real assets do not cover, generate tiny, ephemeral fixtures at runtime (e.g., minimal OBJ/MTL and small normal images via OIIO) under a per-test temp folder and clean up afterward. PRs must not add or modify files under `testdata/`. Never commit any generated assets.

---

## 4) Coding Conventions & Error Handling

* **Style:** Near‑Google C++ (brace on next line OK). Provide `.clang-format`.
* **Naming:** `snake_case` for funcs/vars; `CamelCase` for types; `kConstantCase` for constants; avoid globals.
* **Headers:** Keep light; forward‑declare; PIMPL where heavy.
* **Determinism:** Same inputs/flags → same outputs. Provide a **`--deterministic`** mode (see §11).
* **Errors:** Use `tl::expected<T, N2DError>` for fallible ops. Define `include/n2d/core/error.hpp`:

```cpp
enum class ErrorCode {
  IoError,
  MeshParseError,
  InvalidAsset,
  InvalidArgs,
  IncompatibleTextures,
  AmbiguousInput,
  SolverFailed,
  UserCancelled
};

struct N2DError {
  ErrorCode code;
  std::string message; // one-line actionable description
};
```

Log actionable messages; never crash on malformed assets.

---

## 5) CLI Design

**Binary:** `n2d`

### Examples

```bash
# Inspect materials, UV sets, UDIMs
n2d inspect --mesh model.fbx

# Single-tile bake with explicit normal map
n2d bake --mesh model.fbx --material "Body" \
         --normal body_N.png --uv-set "UVMap" \
         --out out/Body_disp.exr

# UDIM bake (auto-detect 1001..)
n2d bake --mesh model.fbx --material "Cloth" \
         --normal "tex/Cloth_N_<UDIM>.exr" --uv-set "UVSet2" \
         --out "out/Cloth_disp_<UDIM>.exr" --amplitude-mm 2.0
```

### Subcommands & Flags (minimum)

**inspect**

* `--mesh <path>` (required)
* Output: materials (index+name), UV sets per mesh/material, UDIM tiles per material, guessed normal paths, **y‑channel convention guess (+Y/−Y)**, overlap diagnostics.
* `--inspect-json <file>`: write JSON report.

**bake**

* `--mesh <path>` (required)
* `--material <name|index>` (required unless `--normal` supplied per subset)
* `--normal <file|pattern>` (supports `<UDIM>`, `%04d`, OIIO tiled specs)
* `--uv-set <name|index>` (default first/0)
* `--out <file|pattern>` (EXR; supports `<UDIM>`)
* `--amplitude-mm <float>` (post scale in mm)
* `--y-is-up` / `--y-is-down` (default auto)
* `--normalization {xyz,xy,none}` (default: **auto** → `xy` for 2‑ch, `xyz` for 3/4‑ch)
* `--max-slope <float>` (default: 10)
  **Prescriptive:** enforce minimum denominator by clamping `z` before slope calculation: `z = max(z, 1/sqrt(1+S^2))` to prevent numerical instability.
* `--threads <N>` (default: HW concurrency)
* `--export-sidecars` (write mask EXRs + JSON chart table)
* `--cache <dir>` (reuse sidecars)
* `--replicate-udims <list|range>` (e.g. `1001-1008,1010`) replicate a single map across tiles
* `--debug-dumps <dir>` (dump `(du,dv)`, divergence, masks)
* `--deterministic` (stable scheduling & FP; see §11)
* `--interactive` (prompt only when set; else fail fast on ambiguity)

**Behavioral rules:**

* Treat normals as **linear** (no sRGB). Warn if file is tagged sRGB.
* Accept 2, 3, or 4 channels. For 4‑channel, ignore A.
* If mesh UVs span multiple UDIMs but a single `--normal` is provided, fail unless `--replicate-udims` is specified.

---

## 6) Data & Intermediates

### 6.1 UV Chart Sidecars

* **Chart ID mask:** `uint16` EXR per tile: 0 = outside; 1..N = island id.
* **Boundary mask (optional):** `uint8` EXR; 1 = pixel has neighbor with different ID.
* **Chart table (JSON):**

```json
{
  "tile": 1001,
  "uv_set": "UVMap",
  "charts": [
    {"id": 1, "flip_u": false, "flip_v": false, "pixel_count": 12345},
    {"id": 2, "flip_u": true,  "flip_v": false, "pixel_count": 6789}
  ],
  "y_channel": "+Y"  // or "-Y"
}
```

* **Metric map (optional v1 scaffolding):** 2‑channel float EXR
  R = `mm_per_texel_u`, G = `mm_per_texel_v` (B/A = 0)

### 6.2 Output

* **Default:** single‑channel 32‑bit float EXR (`height`).
  **Optional:** RGBA EXR with `R=height, G/B/A=0` via `--rgba`.
* **Metadata:**

  * `n2d:space = tangent`
  * `n2d:metric = true|false`
  * `n2d:units = "mm"|"texel"` (do **not** claim mm unless metric was applied)
  * `n2d:amplitude = <float>`
  * `Software = normal2disp`
  * `SourceMesh = <path>`

---

## 7) Algorithm (Concise)

Per UV island (per UDIM tile):

1. **Decode & orient normal**

   * Map channels from `[0,1] → [-1,1]`: `x=2R-1`, `y=2G-1`, `z=2B-1` (if present).
   * Apply **Y‑flip** (+Y/−Y) to `y`.
   * Apply island **mirroring flips** to `(x,y)` (see §9.3).
   * **Normalization policy** (after flips):

     * **2‑channel:** `z = +sqrt(max(0, 1 - x*x - y*y))`; then `normalize([x,y,z])`.
     * **3/4‑channel:** `normalize([x,y,z])` using all three channels (default).
       `--normalization xy` forces reconstruction instead.

2. **Compute slopes**
   `du = -x / max(z, z_min)`, `dv = -y / max(z, z_min)`
   where `z_min` enforces **`--max-slope`** (S): `z_min = 1/sqrt(1+S*S)`.

3. **Divergence (zero‑flux Neumann)**

   * Use **staggered‑flux** form for stability:
     `div[i,j] = (p[i+1/2,j] - p[i-1/2,j]) + (q[i,j+1/2] - q[i,j-1/2])`,
     with edge fluxes as neighbor averages: `p[i+1/2,j]=0.5*(du[i+1,j]+du[i,j])`, etc.
   * At island boundaries, set the exterior edge flux **equal to the interior** (zero normal flux) rather than dropping the neighbor.

4. **Assemble masked Laplacian (5‑point)**

   * For interior pixel `p`: `A[p,p] = degree_in_mask(p)`; for neighbor `q` in same island: `A[p,q] = -1`.
   * RHS `b[p] = div[p]`.

5. **Fix nullspace**

   * **Anchor** one pixel: pick largest‑`z` (flattest) pixel; **tie‑break** with smallest linear index `(y*W+x)`.
     Overwrite row/col: set `A[p0,*]=0`, `A[p0,p0]=1`, `b[p0]=0`.
   * *(Alternative later)* add a **soft mean‑zero** constraint row.

6. **Solve**

   * Eigen **ConjugateGradient** with **Incomplete Cholesky (IC0)**, or **amgcl** if enabled.
   * Tolerance/iterations configurable.
   * (Optional) **coarse‑to‑fine** pyramid to initialize and cut iterations.

7. **Post**

   * If anchor not used: subtract mean(h).
   * Apply `--amplitude-mm`.
   * Write EXR; fill sidecars/metadata.

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
Deliver: island rasterization, mirrored flags, sidecars.
Tests: island counts; mirrored shell detection; coverage vs UV area (±1%).

**P4 — Normals→Slopes & Divergence**
Deliver: normalization policy, slope computation, staggered divergence with zero‑flux boundaries.
Tests: flat normal → div≈0; synthetic gradient → analytic divergence match.

**P5 — Poisson Solver**
Deliver: masked Laplacian, anchor/mean‑zero, IC0 preconditioner.
Tests: procedural height → normals → recon RMS < 1e‑3 (512²); border flatness (no bevel); solver convergence.

**P6 — UDIM & Multithreading**
Deliver: per‑tile + per‑island parallelism; deterministic reductions.
Tests: two tiles with different island counts; `--threads 1` vs default identical hashes.

**P7 — CLI Polish & Sidecars**
Deliver: `inspect --inspect-json`; bake flags incl. `--normalization`, `--max-slope`, `--replicate-udims`, `--debug-dumps`, `--deterministic`.
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

## 8A) Local Test Data & PR Rules

* Local assets live under `testdata/` and are expected to be present. Prefer them for integration and I/O smoke tests (e.g., `inspect` on a real mesh, Y-channel guess behavior).
* Some tests intentionally generate fixtures even when `testdata/` exists (analytic golden cases, BC5 2-channel normals, UDIM tile synthesis, mirrored/overlap UVs, donut islands, near-horizon normals). Keep generated fixtures tiny (≤16×16) and sparse.
* Discovery and selection: if environment variable `N2D_TESTDATA_DIR` is set, use that; otherwise use repository `testdata/`. Tests may still choose runtime-generated fixtures when the scenario requires ground truth or specific edge conditions.
* PR policy: do not add or modify binary assets under `testdata/`. Do not commit any generated files. Ensure `.gitignore` ignores `build/` and patterns like `*.exr`, `*.png`, `*.tif`, `*.bin`.
* Generated fixtures must be written to a per-test temporary directory (e.g., under the CTest working dir) and deleted at test end. Never write inside the repo tree.
* When available, use the head normal map in `testdata/` to exercise +Y/−Y decoding and normalization policies.

---

## 9) Detailed Implementation Notes

### 9.1 Normal Decoding & Normalization

* Accept 2/3/4 channels.
* Map to `[-1,1]` first.
* Apply **Y convention** and **island flips** to `(x,y)` **before** normalization.
* **Policy (after flips):**

  * **2‑channel**: compute `z = +sqrt(max(0, 1 - x*x - y*y))`; then `normalize([x,y,z])`.

**Note:** In the 2-channel path the vector is unit-length by construction; the extra `normalize([x,y,z])` is an intentional "belt-and-suspenders" to counter small floating-point drift and keep parity with the 3-channel path.

* **3/4‑channel**: compute `z = 2*B-1`, then `normalize([x,y,z])` (default).
  `--normalization xy` forces reconstruction; `--normalization none` uses decoded values as‑is (debug).
* Guard with **`--max-slope`**: either clamp `du,dv` or enforce `z_min`.

### 9.2 Poisson Assembly

* Masked 5‑point stencil; divergence as in §7.3; Neumann via edge flux reuse.
* **Anchor** by overwriting row/col (`A[p0,*]=0`, `A[p0,p0]=1`, `b[p0]=0`).
  *Alternative:* soft mean‑zero constraint with tiny weight preserves SPD.
* **Preconditioning:** prefer **IC0**; optionally **amgcl** if available.

### 9.3 Mirroring & Orientation

* Detect island handedness by comparing signed 3D area vs signed UV area (per triangle).
  Aggregate per island (majority vote); emit sidecar booleans `flip_u`, `flip_v`.
* Flip the corresponding component(s) of `(x,y)` **prior** to normalization.

### 9.4 Multithreading

* Parallelize per tile → per island (TBB).
  Avoid write races by writing to island‑local buffers then compositing.
* Ensure deterministic reductions when `--deterministic` is set (static partitioner, ordered merges).

### 9.5 UDIM & UV Sets

* Detect `<UDIM>` in paths; else infer from UV ranges.
  If single map + multi‑tile UVs → require `--replicate-udims`.

### 9.6 I/O Details

* Normals are **linear**. If sRGB tag detected → warn and proceed as linear.
* Accept PNG/TIF/EXR; EXR recommended for float precision.

---

## 10) API Surface (for future GUI)

```cpp
// n2d::NormalizationMode
enum class NormalizationMode { Auto, XYZ, XY, None };

struct BakeParams {
  std::filesystem::path mesh_path;
  std::string material_selector;      // name or index string
  std::string uv_set;                 // name or index string
  std::string normal_pattern;         // supports <UDIM>
  std::string output_pattern;         // supports <UDIM>
  bool y_is_down = false;             // +Y default
  bool export_sidecars = false;
  std::optional<std::filesystem::path> cache_dir;
  int threads = (int)std::thread::hardware_concurrency();
  float amplitude_mm = 1.0f;
  float max_slope = 10.0f;
  NormalizationMode normalization = NormalizationMode::Auto;
  bool deterministic = false;
  std::optional<std::filesystem::path> debug_dumps_dir;
};

struct Result {
  std::vector<std::filesystem::path> outputs;  // EXR paths
  std::vector<std::string> log_lines;
};

tl::expected<Result, N2DError> Bake(const BakeParams& p);
```

Keep all state inside the call; thread‑safe.

---

## 11) Determinism & Performance

* `--deterministic` enforces: static partitioning; stable reduction order; `/fp:strict` and FMA‑stable math where applicable; optional single‑threaded solver if bitwise identity is mandatory.
* Performance guidance: 4K map with \~10 islands should complete in **seconds** with IC0 or amgcl and/or coarse‑to‑fine initialization on an 8‑core CPU. Document expected ranges.

---

## 12) Acceptance Criteria

* **Correctness:** unit tests pass; no bevel at island borders; synthetic golden RMS threshold met.
* **Robustness:** supports OBJ/FBX/GLTF; normal maps PNG/TIF/EXR; UDIM 1001–1010; handles 2/3/4‑ch normals; overlap diagnostics.
* **UX:** `inspect` reports are clear; `bake` errors actionable; `--interactive` resolves common ambiguities.
* **Determinism:** identical hashes across runs when `--deterministic`; `--threads 1` and default produce identical outputs in practice for most scenes (document caveats).
* **Output:** EXR written; metadata truthful (`metric`, `units`).

---

## 13) Phase‑wise "Done" Checks (Quick)

* **P0:** VS builds; `--help` works.
* **P1:** Materials & UV sets; UDIMs; overlaps.
* **P2:** Normal maps load (2/3/4‑ch); UDIM expansion; sRGB warnings.
* **P3:** Chart masks/flags exported; coverage sane.
* **P4:** Divergence centered; flat → near‑zero.
* **P5:** Poisson passes golden RMS; border remains flat.
* **P6:** UDIM end‑to‑end; multi‑thread determinism.
* **P7:** CLI polish; EXR metadata; JSON schema validated.
* **P8:** Metric map optional path proven.

---

## 15) CI & Container (Docker) Setup

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

**Common CI pitfalls:**

* “invalid builtin-baseline”: the baseline SHA in the manifest/config does not exist in the container’s vcpkg clone. Fetch/checkout that SHA before bootstrap.
* Running `vcpkg install` manually: unnecessary in manifest mode—CMake + toolchain will drive installation.
* Missing binary cache: enable `VCPKG_FEATURE_FLAGS=binarycaching` and set a persistent `VCPKG_DEFAULT_BINARY_CACHE` to speed up builds.

---

## 14) Glossary

* **Neumann (zero‑flux) boundary:** derivative normal to boundary is zero; implemented by reusing interior edge flux in divergence at the boundary.
* **Staggered‑flux divergence:** difference of face‑centered fluxes; more stable pairing with 5‑point Laplacian than naive 4‑neighbor differences.
* **BC5/ATI2N:** 2‑channel compressed normal format; requires reconstructing Z from X,Y.
* **IC0:** Incomplete Cholesky preconditioner; accelerates CG for SPD systems.
