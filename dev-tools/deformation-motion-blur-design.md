# Deformation (vertex) motion blur — E9 design

Status: **Phase 1 (plumbing) implemented.** Roadmap item E9 — engine-level
per-vertex motion blur for meshes and curves. Adapter-side prerequisites
(A5 step-collection infra) are done; backend acceleration is not started —
meshes carrying a vertex series still render their static `vertices`
(static fallback, same stance as E7 curve data in serialized scenes).

Phase-1 surface (implemented):

- `ExtTriangleMesh::SetVertexMotion(times, stepVerts)` — `motionVertTimes`
  + `motionVertSteps` members. Validation: ≥2 steps, strictly increasing
  finite times, every step buffer has the mesh vertex count (constant
  topology). `GetVertexAtTime(i, t)` lerps between adjacent steps and
  clamps outside the range (MotionSystem convention). Positions only;
  per-step normals are not stored.
- Propagation: `ApplyTransform` transforms every step buffer;
  `CopyExt`/`Copy` keep the series iff `meshVertices` is not overridden
  (same rule as curve data); `Merge` keeps it iff every input has a
  series with identical times — partial presence or differing times
  throw (same convention as UV/AOV mismatches). `Delete` frees steps;
  serialization does not persist them (load clears → static fallback).
- Scene plumbing: `Scene::SetMeshVertexMotion` →
  `ExtMeshCache::SetMeshVertexMotion` (plain `TYPE_EXT_TRIANGLE` meshes
  only — set it on the base shape, not instance/motion wrappers) +
  `GEOMETRY_EDIT` so the next `Preprocess` rebuilds the DataSet.
- Public API: `luxcore::Scene::SetMeshVertexMotion(meshName, times,
  timesCount, verts, vertsCount)` — flat step-major float array —
  and the pyluxcore binding
  `Scene.SetMeshVertexMotion(name, times, [step (N,3) arrays])`.
- Properties: `<prefix>.motion.N.time` + `<prefix>.motion.N.vertices`
  parsed in `CreateInlinedMesh`, covering both `scene.objects.*` inlined
  meshes and `scene.shapes.* type=inlinedmesh`. The object-level
  transform-motion wrapper is now skipped when no step defines
  `.transformation` (vertex-only motion keeps the plain-mesh path —
  avoids pushing static transforms through the motion path).
- Unit test: `vertexmotion_test` (`dev-tools/e9_vertexmotion_test.cpp`,
  wired in `src/luxrays/CMakeLists.txt`) — 22 asserts covering
  validation, lerp/clamp, copy/merge/transform/serialization rules.
  Building it surfaced a pre-existing upstream bug —
  `TriangleMesh::save()` serialized `vertices.Count()` as the triangle
  count, corrupting every mesh where the two differ; fixed separately.

## Why this is an engine task

LuxCore's `MotionSystem` (`include/luxrays/core/geometry/motionsystem.h`)
interpolates `Transform`s only — `InterpolatedTransform` decomposes two
matrices and samples between them. There is no per-vertex time data
anywhere in the mesh model, BVH, or scene properties. Scenes that deform
vertices (character animation, FLIP/cloth caches, GN-evaluated geometry,
hair driven by armature) currently render with zero motion blur.

## Data model

`ExtTriangleMesh` gains an optional vertex time series:

- `motionTimes` — N shutter times (already the convention used by
  `MotionSystem` and by the `motion.N.*` scene properties).
- `motionVerts[N]` — N vertex buffers, one per step. Positions only;
  normals are recomputed at shading time (or per-step normals if a cheap
  cross-product normal is insufficient — decide during implementation).
- Same linear-interpolation model as transform motion: at ray time t,
  vertex position is the lerp between the surrounding steps.
- Topology is constant across steps — Blender-side export enforces this
  (vertex-count mismatch falls back to the center step, same policy as
  the A5/point-cloud instancing paths).

Memory cost: N× the vertex array. For hair (strands.cpp) the curve
control points get the same treatment — `curveCPs` becomes
`curveCPs[N]`; tessellated fallback gets `motionVerts` the same way.

## Per-backend path

| Backend | Mechanism |
|---|---|
| Metal HWRT | `MTLAccelerationStructureMotionTriangleGeometryDescriptor` and `MTLAccelerationStructureMotionCurveGeometryDescriptor` — per-keyframe vertex buffers, HW interpolates. Gate on `device.supportsMotionBlur`; falls back to static otherwise. Extends the E7 curve-AS path directly. |
| Embree (CPU) | `rtcSetGeometryTimeStepCount` + `RTC_BUFFER_TYPE_VERTEX` per timestep — native multi-segment motion blur for meshes and curves. |
| OpenCL (SW) | No HW motion support: build the BVH over each triangle's *swept* AABB (union of all step positions), and interpolate verts to `ray.time` inside the leaf intersection test. Same asymptotic traversal, wider bounds. The .cl kernels already carry `ray.time` for transform motion. |
| OptiX/CUDA | `OptixMotionGeometryDesc` vertex buffers — E8-era refresh handles this; out of scope for the first implementation. |

## Scene/API surface

- `Scene::DefineMesh`/`DefineMeshExt` gains an optional
  `(times, vertsPerStep)` argument, or a `motion.N.vertices` property —
  properties version preferred since everything else is property-driven
  (`motion.N.time`, `motion.N.transformation`).
- `Mesh::ApplyTransform`/`Merge`/`CopyExt` must propagate the vertex
  series under the same rules the E7 curve data uses (kept iff vertices
  are not overridden; merge keeps motion iff all inputs carry it).
- Serialization: same stance as E7 curve data — v1 does not persist
  vertex series; deserialized scenes render static.

## BlendLuxCore side (adapter, mostly done)

- The A5 step loop in `export/motion_blur.py` already re-evaluates the
  depsgraph at every shutter step — deformation export adds
  `foreach_get("co")` on the evaluated mesh per step (vectorized, same
  fast path used for hair point collection).
- Gate: `enable_motion_blur` on the object; topology check per step.
- Hair: `depsgraph` particle/curve positions per step — same collection
  point as A5's dupli matrices.

## Phasing

1. ~~ExtTriangleMesh vertex series + serialization/merge rules + scene
   property plumbing (no backend uses it yet — pure plumbing).~~ **Done**
   — see the Phase-1 surface list above.
2. Metal motion geometry descriptor (primary GPU-first target).
3. OpenCL swept-bound software path.
4. Embree timestep path (CPU parity).
5. BlendLuxCore mesh + hair export.
6. Validation scenes: animated character mesh, GN-deformed geometry,
   armature-driven hair — A/B vs static, plus a divergence-stress scene.

Open question: whether per-step normals are needed for smooth-shaded
deforming surfaces or whether shading-time recomputation suffices —
decide after the first visual tests.
