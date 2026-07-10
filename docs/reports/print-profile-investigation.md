# Documentation and print-profile audit

Audit date: 2026-07-10. This report records the evidence used to correct the
README and investigate two print-profile concerns. It describes the current
implementation; no production algorithm or default was changed.

> Historical note: the VTP validation issue identified in this report was fixed
> in the subsequent low-risk workflow tranche. The two print-profile algorithm
> concerns remain documented limitations pending a separate geometry task.

## Baseline and commands

Repository state before the documentation patch: `8fb6668`. The workspace
already contained an unrelated untracked `poetry.lock`; it was not read as
authoritative metadata and was not modified.

The environment initially had no Python project dependencies. After
`python -m pip install -e '.[dev]'`, the first test collection failed because
PyMeshLab could not load `libGL.so.1`. Installing the same Ubuntu packages as CI
resolved it:

```console
sudo apt-get update
sudo apt-get install -y libgl1 libopengl0
python -m pytest -q -ra
```

Environment and successful result:

| item | audited value |
|---|---|
| platform | Linux 6.8.0-1052-azure, x86-64, glibc 2.39 (Ubuntu 24.04 package base) |
| Python | 3.12.1, GCC 13.3.0 |
| NumPy / SciPy | 2.5.1 / 1.18.0 |
| pydicom / SimpleITK | 3.0.2 / 2.5.5 |
| VTK / Trimesh | 9.6.2 / 4.12.2 |
| PyMeshLab | 2025.7.post1 |
| pytest | 9.1.1 |
| OpenGL packages | `libgl1 1.7.0-1build1`, `libopengl0 1.7.0-1build1` |
| result | 160 passed in 36.83 s; 0 skipped; 398 upstream NumPy-shape deprecation warnings |

CLI help was captured for all six subcommands, and every major README command
was parsed with `dicom_surface.cli.build_parser()`: `list`, `presets`, `convert`,
`merge`, `validate`, and `repair` all matched the parser.

No documentation linter or link-check configuration exists in the repository.

### Implementation map used for the audit

| concern | implementation |
|---|---|
| CLI entry point and all subcommands | `src/dicom_surface/__main__.py`; parser and `list`/`presets`/`convert`/`merge`/`validate`/`repair` handlers in `cli.py` |
| Tissue presets and print profiles | `presets.py`; CLI composition in `_resolve_print_profile`; runtime composition in `pipeline.build_mask` |
| Series discovery/grouping/ranking | `series.discover`, `Series.unusable_reason`, `rank`, and `select` |
| DICOM volume loading and HU reliance | `volume.load`, `has_calibrated_hu`, and `warnings_for` |
| Kernel rounding and morphology | `geometry.py`; threshold/filter/grid/thin-selector/resampling functions in `segment.py` |
| Conversion pipeline/provenance | `pipeline.build_mask` and `pipeline.convert` |
| Flying Edges, component filtering, smoothing, decimation, normals, writers | `surface.py` |
| File validation and welding-dependent metrics | `validate.py`; optional external repair in `repair.py` |
| FFT/ICP registration and metrics | `registration.py` |
| Patient comparison, gates, occupancy union, merge provenance | `merge.py` |
| Metadata/dependency/OS support | `pyproject.toml` and `.github/workflows/ci.yml` |

Relevant regression evidence is in `tests/test_geometry.py`, `test_series.py`,
`test_segment.py`, `test_surface.py`, `test_pipeline.py`, `test_merge.py`, and
`test_printing.py`. Git history was inspected for measurements whose scripts and
fixtures are no longer present.

## Claim-audit summary

| Topic | What the README said | What the code/evidence does | Documentation action |
|---|---|---|---|
| Print-profile table | Separate `thicken` and minimum-feature values, including nonexistent resin 0.4 mm thickening | `PrintProfile` has only `closing_mm`, `min_island_mm3`, and `min_feature_mm`; `T/2` is derived as the selective dilation radius | Removed `thicken`; renamed the target and explained derivation and CLI precedence |
| Topology-safe resampling | Resampling was topology-safe and marching cubes always manifold | `vtkFlyingEdges3D` at 0.5 is used; VTK warns of degenerate triangles; tests and history show structures/tunnels/genus can change | Distinguished output validity from preservation of mask/anatomical topology and retained post-write validation |
| Smoothing table | Implied deviation from a wholly unsmoothed pipeline | The zero post-smooth row still follows 20 initial smoothing iterations, component selection, and 600k decimation; benchmark artifacts are absent | Named preceding stages and marked the table historical/non-reproducible |
| Identity gates | Mixed “conflicting demographics”, “only differing IDs refuse”, and “without complaint” | Name, birth-date, or sex conflict refuses; ID/issuer/name rules are ordered; missing/de-identified data warns; `--force` has bounded scope | Added an ordered decision table and replaced “without complaint” |
| Privacy | Patient-field values never surface | Those five values are discarded, but UID, description, paths, and other linkable/free-text metadata can surface | Kept the tested patient-value statement and added residual metadata risk |
| Minimum wall | Claimed exact/no-wall-thinner-than target and printed-model guarantee | Integer morphology acts on a resampled binary mask; discrete support differs from `2r×spacing`; no final-mesh thickness test exists | Reframed as an intended mask-space feature target, with no final STL/manufacturing guarantee |
| Through-plane sampling | A 0.6 mm wall could not exist with 0.8 mm slices | Reliability depends on orientation, slice response, reconstruction, and partial volume; sub-pitch anatomy is not reliably inferred in that direction | Replaced the categorical statement |
| Sharp kernels | Used unexplained `[BHUY]r?[6-9]x` shorthand | Literal case-insensitive regex is `(?:^|[^0-9])(?:[BHUY]r?|BONE|EDGE|LUNG)\s*_?([6-9]\d)`; no manufacturer normalization | Printed the literal rule, verified `Hr68`/`B70f`, and documented heuristic errors |
| Linux/OpenGL | Desktop Linux always had requirements and both packages only affected decimation | Minimal audit image failed PyMeshLab import without `libGL.so.1`; CI history ties `libopengl0` to the meshing plugin | Used bounded distro-specific wording and separated failure modes |
| STL welding | Exact-coordinate weld and rhetorical 100%-boundary absolute | Trimesh processing rounds at digits derived from `tol.merge` (`1e-8` here); default processing does not remove duplicate faces or fix winding | Documented version-dependent approximate welding and measurement semantics |
| DICOM frames | Studies rarely shared a usable frame | `FrameOfReferenceUID` is never read; registration always runs | Stated implementation behavior without a prevalence claim |
| Decimation target | MeshLab always hit the target exactly and preserved every mesh | Filter constraints can prevent a target; tests cover specific synthetics; historical skull assets are absent | Qualified as a target/observation and added a current-version torus check |
| DICOM support | Broadly implied DICOM-series compatibility | Classic single-frame path is tested; enhanced multi-frame/mosaics are unsupported; compression and several pixel/geometry cases are untested or delegated | Added a compatibility/limitations table |
| Registration | Mixed method narrative and unexplained overlap/Dice/RMS | 2 mm 3D FFT translation, two-pass trimmed point-to-plane ICP, 4 mm correspondence tolerance, 1.5 mm shared-FOV Dice grid | Added precise metric and sampling definitions |
| Registration gates | “Calibrated” four positive pairs and one impostor | Positives are repeat studies/reconstructions of one skull; one negative; thresholds were selected on that set; artifacts/protocols are missing | Relabelled exploratory empirical selection and scoped the table |
| Quantitative claims | Presented development results with broad precision | Real-data scripts/scans/versions/hardware are absent; synthetic tests reproduce narrower properties | Labelled historical values, removed implied bounds, and added an evidence section |
| Byte identity | Universal byte-for-byte default output | Test asserts mask equality only; manual historical checks are not retained and did not cover every format/platform | Limited claim to mask-stage identity |
| Self-intersections | Listed among validation results without explaining opt-in status | CLI checks only with `--self-intersections`; edge validity/watertightness do not exclude intersections | Documented opt-in scope and stage coverage |
| CLI terminology/examples | Examples were mostly correct but implied equal end-to-end support for every output extension | Parser accepts all shown forms; measured inter-slice distance can differ from `SliceThickness`; Trimesh validation rejects VTP | Parser-checked examples, corrected spacing terms, and documented the VTP `--no-validate` workaround |
| Merge-profile rationale | Complete per-scan transform treated as distributive | Only dilation distributes; closing, island removal, grid conversion, and selector do not | Corrected rationale and linked the order experiment below |

## Concern 1: thin features near thick anatomy

### Exact operation examined

For feature target `T`, after the effective closing and any printability-grid
conversion, production executes:

```text
r         = T / 2
core      = binary opening(mask, ceil(r / spacing), sitkBall)
reachable = binary dilation(core, ceil(r / spacing), sitkBall)
thin      = mask AND NOT reachable
output    = mask OR binary dilation(thin, ceil(r / spacing), sitkBall)
```

The opening and dilation are 3-D SimpleITK ball operations with foreground value
1. A profile runs after the preset's threshold, first island filter, median,
optional opening, and effective closing; it runs before the second island filter,
padding, 0.5-isovalue Flying Edges extraction, initial smoothing, surface-shell
filter, optional PyMeshLab decimation, and post-smoothing.

The resin/FDM targets are 0.6/1.2 mm; their closing extents are 3.2/4.8 mm and
island floors 100/200 mm³. The focused selector cases below use the FDM target
on its exact selected 0.3 mm isotropic grid (`r=2` voxels). A second source grid
of 0.4×0.4×1.2 mm was converted through production
`printability_grid_mm`/`to_printability_grid` to the same 0.3 mm work grid.
Inputs are deliberately post-closing masks so the selector is isolated.

### Diagnostic method

A temporary, uncommitted `audit_print_experiments.py` imported production
`segment`/`surface` functions and was run as:

```console
python -m py_compile audit_print_experiments.py
python audit_print_experiments.py
```

It generated a long and short attached fin, long plate bridge, short bridge,
near-parallel membrane tethered to a wall, one- and two-voxel oblique walls,
curved thin shell attached to a thick shell, and a skull-like plate/rim. It
recorded selector coverage, added mask volume, analytical opposing-axis run
thickness where defined, total volume/bounds, final-mesh nearest-vertex symmetric
surface distance, edge validity, components/genus, and PyMeshLab intersecting-face
selection. Each mesh received the `bone` surface stage: 20 initial iterations,
600k target (a no-op because these meshes were already smaller), and 25
post-iterations. Cross-sections were visually inspected before the script was
deleted.

`selected` below is the fraction of the labelled target feature in production's
actual `thin` selector. Thickness is direct opposing-axis support at a controlled
probe, not `mask \ opening`.

| post-closing case | native grid | selected | probe thickness in→out | total volume growth | final surface p95 | observation |
|---|---:|---:|---:|---:|---:|---|
| long fin | 0.3³ | 93.9% | 0.3→1.5 mm | 7.5% | 0.66 mm | Fin visibly thickened except its core-reach attachment band |
| short fin (<0.6 mm reach) | 0.3³ | **0.0%** | **0.3→0.3 mm** | **0.0%** | 0 | Visibly unchanged |
| long bridge | 0.3³ | 85.2% | 0.3→1.5 mm | 2.4% | 0.59 mm | Central bridge thickened |
| short bridge | 0.3³ | **0.0%** | **0.3→0.3 mm** | **0.0%** | 0 | Visibly unchanged between thick blocks |
| parallel membrane | 0.3³ | 95.1% | 0.3→5.7 mm on the sampled joined run | 20.2% | 1.16 mm | Membrane grew and joined the nearby wall |
| one-/two-voxel oblique walls | 0.3³ | 95.7% / 93.9% | analytical input width 0.44/0.84 mm | 90.6% / 77.2% | not comparable / 0.72 mm | Both visibly grew; one-voxel component filtering makes its whole-mesh distance unsuitable |
| curved attached shell | 0.3³ | 79.2% | radial probe only | 55.4% | not comparable | Most thin hemisphere grew; attachment region remained in core reach |
| skull plate + bulky rim | 0.3³ | 96.9% | 0.3→1.5 mm | 31.0% | 0.60 mm | Representative long plate grew; two intersecting faces appeared after the surface stage |
| long fin / bridge | 0.4×0.4×1.2 → 0.3³ | 91.0% / 84.0% | 0.3→1.5 mm | 6.4% / 2.3% | 0.66 / 0.59 mm | Same practical result after anisotropic-source conversion |
| short fin / bridge | 0.4×0.4×1.2 → 0.3³ | **0.0% / 0.0%** | bridge 0.3→0.3 mm | **0.0% / 0.0%** | 0 | Same reach blind spot |
| curved shell | 0.4×0.4×1.2 → 0.3³ | 50.1% | radial probe only | 54.0% | not comparable | Partial-volume resampling changed which shell material was selected |
| skull plate + rim | 0.4×0.4×1.2 → 0.3³ | 0.0% | 1.5→1.5 mm | 0.0% | 0 | Resampling already made the controlled plate thicker than the target |

All final meshes were edge-watertight, with zero boundary/non-manifold edges and
one retained component. That did not imply intersection-free output: several
sharp input phantoms already intersected after smoothing, and the isotropic
skull-plate output changed from 0 to 2 PyMeshLab-selected intersecting faces.

The global post-profile `thin_fraction` for both unchanged short cases was about
1%, so `thin_material_warning` returned no warning under its 5% threshold. The
miss can therefore be silent even though the localized feature is four times
thinner than the requested target.

### Assessment and recommendation

The mechanism is confirmed, not hypothetical: a short feature inside the
dilated core's reach is intentionally classified as reachable and cannot seed
growth. Representative long fins, bridges, oblique walls, and the isotropic
skull-like plate were thickened, so the selector still avoids the much larger
global inflation of uniform dilation. Practical severity is **moderate but not
yet quantified on patient anatomy**: the failure is large and silent on controlled
short fin/bridge geometry, but no orbital/septal patient fixture or physical
print was available.

Recommendation: document the limitation now, add focused permanent tests in a
separate implementation task, and treat a better detector/final validator as a
future implementation fix rather than leaving the current “minimum” implication
unqualified. Candidate approaches include a local-thickness/maximal-inscribed-
sphere transform, medial distance analysis, directional opposing-surface distance,
constrained growth from an independently detected thin set, or final-mesh
thickness validation. Regression risks are global inflation from corner speckle,
closing holes/tunnels, merging nearby anatomy, high memory use on 300M-voxel
grids, anisotropic/partial-volume instability, and extra self-intersections.

## Concern 2: profile each scan before union

### Exact current order

Registration always uses anatomical masks. Only after identity and registration
gates pass does `merge` rebuild each mask with the profile. Thus the complete
current transform is:

```text
P(X) = threshold → island floor → median → optional preset opening
       → max(preset, profile) closing → optional printability-grid conversion
       → opened-core selection and selective dilation → island floor

current = max(linear_resample(antialias(P(A))),
              linear_resample(antialias(P(B))))
```

The max is on float occupancy fields, followed by 0.5-isovalue meshing. There is
no post-union morphology. Only the final dilation primitive distributes over a
binary union; `P` as a whole does not.

### Diagnostic method and results

The same temporary script compared current order with diagnostic-only
`P(threshold(max(resample(A), resample(B)), 0.5))` on a 0.3 mm common grid using
the actual FDM closing, island, feature, grid, occupancy, and surface-stage
values. Tissue-preset median/opening were no-ops so the processing-order effect
was isolated; every individual input exceeded the 200 mm³ island floor. Final
surfaces received the same `bone` smoothing/target stage as above.

| case | mask difference (% of anatomical union) | current vs alternative local thickness | final surface p95 / max | topology/visual result |
|---|---:|---:|---:|---|
| Identical observations with partial fields | 0% | 1.5 / 1.5 mm | 0 / 0 | Identical |
| Complementary halves of a uniformly thin wall | 0% | 1.5 / 1.5 mm | 0 / 0 | Identical; dilation distributivity applies because both selectors chose the wall |
| Two complementary thin layers whose five-layer union is already thick | **91.8%** | **2.7 / 1.5 mm** | **0.66 / 0.80 mm** | Per-scan result visibly over-thick; 14 intersecting faces versus 0 |
| Same thin plate observed with a two-voxel offset and one-layer overlap | **95.0%** | **2.7 / 1.5 mm** | **0.67 / 0.80 mm** | Same visible over-thickening; 4 intersecting faces versus 0 |
| Thick rim in A plus adjoining plate in B | 2.7% | 2.1 / 2.1 mm | 0.38 / 0.54 mm | Small localized differences at plate/rim junctions; both valid |
| Partial-overlap shell with one-voxel registration error | 0.1% | not sampled | 0.007 / 0.21 mm | Negligible at printer scale in this case |
| 3.6 mm inter-scan seam under the 4.8 mm FDM closing extent | 9.1% | gap open / seam joined | 16.59 / 16.80 mm | Current mask stayed two components and final largest-shell filtering discarded one; fusion-first added 241 mm³ and joined them |
| Half-voxel-offset oblique fractional boundary | 0% | not sampled | 0 / 0 | Identical in this constructed case; both had 12 intersecting faces after smoothing |

All compared final meshes were edge-watertight with zero boundary/non-manifold
edges. The seam case's final one-component metrics hide an important difference:
current order had two mask components before `largest_component`, so half of the
geometry was dropped, whereas fusion-first closing joined and retained both. That
joining could repair a true seam or incorrectly bridge unrelated surfaces.

### Assessment and recommendation

The current order is not uniformly worse: four cases were identical, the
registration-error shell differed negligibly, and per-scan processing avoids
closing across a gap that may be anatomical rather than a seam. The diagnostic
fractional case did not show a current-order advantage, but one case cannot prove
that preserving each scan's post-profile occupancy is valueless.

The overlaid-layer and offset-plate counterexamples are nevertheless material at
printer scale: 1.2 mm of unnecessary mask thickness and about 0.66 mm p95 surface
displacement, plus self-intersections not present fusion-first. Recommendation:
open a future implementation task to compare post-union thin selection or a full
fusion-first profile on representative registered skull pairs. Do not switch
orders solely from these synthetics: fusion-first closing can join seams and
discarding fractional occupancy can reduce boundary fidelity.

Registration inputs and gates should remain anatomical and therefore need not
change. Performance could improve by running expensive print morphology once,
but the fused extent can approach the 800M-voxel merge limit, versus the current
per-scan 300M printability budgets; peak memory could become substantially worse.
A hybrid could keep per-scan closing/occupancy but select thin material after
union. It would need explicit overlap handling and final thickness/intersection
validation.

## Decimation spot check

The README's historical skull numbers were not reproducible, so the audit also
ran a current-version procedural check:

```text
trimesh.creation.torus(major_radius=20, minor_radius=2,
                       major_sections=768, minor_sections=512)
surface.decimate(poly, target)
```

The input had 786,432 faces, genus 1, and volume 1,579.079 mm³. With PyMeshLab
2025.7.post1, targets 600,000 and 250,000 produced exactly those counts, remained
watertight with genus 1 and zero boundary/non-manifold edges, and changed volume
by -0.00165% and -0.01544%. A 1,000,000 target returned the original object
without modification. This is one procedural torus, not a guarantee for
constrained skull or thin-wall meshes.

## Intentional non-changes

The audit intentionally left all of the following unchanged:

- DICOM discovery, orientation grouping, ranking, loading, intensity handling,
  and series-selection rules;
- tissue presets, print-profile values, CLI options/defaults, thresholds,
  morphology extents, island floors, component policies, and voxel budgets;
- thresholding, median/opening/closing order, selective thin-set definition,
  dilation rounding, printability-grid selection, antialiasing, occupancy
  threshold, padding, and field-of-view capping;
- VTK Flying Edges extraction/isovalue, coordinate transforms, normals,
  smoothing iterations/passbands, component filtering, PyMeshLab filter/options,
  decimation targets, and post-smoothing;
- registration lattice, sampling/seed, ICP trimming/tolerance/iterations,
  transform direction, overlap/Dice/RMS definitions, registration thresholds,
  identity logic, `--force` scope, and unmodified-anatomy registration;
- per-scan profile-before-union order, fractional max union, fused grid default
  and limit, output formats, validation behavior, repair behavior, and
  provenance contents.

Only `README.md` and this report were changed.

## Suggested follow-up issues

### 1. Detect or validate short thin features hidden by core reach

- **Observed failure:** a 0.3 mm fin and bridge within the FDM core's 0.6 mm
  dilation reach were selected 0%, stayed 0.3 mm, and did not trigger the 5%
  warning.
- **Reproduction:** on a 0.3 mm isotropic mask, attach a one-voxel fin shorter
  than two voxels to a large block, or connect two blocks with a one-voxel
  cross-section bridge shorter than the reach; run `segment.thicken(mask, 1.2)`.
- **Expected behavior:** either grow the localized feature toward the documented
  target or report that it remains below the target after final meshing.
- **Potential approaches:** independent local-thickness/medial-distance
  detection, directional opposing-surface distance, constrained dilation, and/or
  final-mesh thickness validation.
- **Regression risks:** global surface inflation, closure of anatomical pores,
  component/tunnel changes, high memory/runtime, anisotropy sensitivity, and
  self-intersections.

### 2. Evaluate thin selection after registered occupancy fusion

- **Observed failure:** complementary/offset observations whose union was 1.5 mm
  thick became 2.7 mm with current per-scan FDM processing; p95 surface distance
  was about 0.66 mm and only current order self-intersected.
- **Reproduction:** create two 0.3 mm-grid plates, each 2–3 layers thick, offset
  so their union is five adjacent layers; compare `P(A)∪P(B)` with `P(A∪B)`.
- **Expected behavior:** avoid redundant thickening in overlap without sealing
  an anatomical seam or losing fractional-boundary fidelity.
- **Potential approaches:** fuse occupancy then apply all of `P`; retain per-scan
  closing but select/grow after union; hybrid per-scan plus post-union validation;
  or explicitly resolve overlap before morphology.
- **Regression risks:** an 800M fused-grid memory peak, different island/component
  behavior, seam closure, fractional-occupancy loss, runtime changes, and altered
  intersection/genus outcomes. Registration must continue to use unmodified
  anatomical masks so gate calibration is not weakened.

### 3. Add final print-profile thickness and self-intersection regression checks

- **Observed failure:** current diagnostics measure a global mask-opening proxy,
  not final mesh thickness; focused outputs could remain thin without warning,
  and smoothing/profile order produced PyMeshLab-selected intersecting faces on
  otherwise watertight meshes.
- **Reproduction:** reuse the short-fin/bridge and offset-layer cases above and
  evaluate after the complete surface stage.
- **Expected behavior:** tests should fail when the final feature silently misses
  its intended target or acquires self-intersections.
- **Potential approaches:** controlled opposing-surface distances for synthetic
  tests, a bounded mesh-thickness estimator, and opt-in or sampled intersection
  checks.
- **Regression risks:** expensive/flaky geometric predicates, tolerance and
  version sensitivity, and pressure to alter production thresholds merely to
  satisfy a synthetic oracle.

### 4. Resolved: make VTP output compatible with default validation

- **Observed failure:** `surface.write` accepts `.vtp`, but both default
  post-conversion validation and `dicom-surface validate` call Trimesh, which
  raises `NotImplementedError: file_type 'vtp' not supported` in 4.12.2.
- **Reproduction:** write the synthetic sphere used by surface tests to VTP and
  call `validate.validate(path, self_intersections=False)`.
- **Expected behavior:** every advertised output extension should either complete
  the default conversion/validation workflow or be rejected early with an
  actionable message.
- **Potential approaches:** validate VTP through VTK/PyMeshLab, convert only for
  measurement without rewriting the file, or make `--no-validate` automatic with
  a clear warning for unsupported validator formats.
- **Regression risks:** metric differences between readers/welding rules,
  inconsistent self-intersection support, and extension-specific code paths.

Resolution: VTP is now read and triangulated with VTK before entering the same
Trimesh metric pipeline. The VTP self-intersection path passes the in-memory mesh
to PyMeshLab. Both paths have regression tests; see
`docs/bugs/vtp-validation.md`.
