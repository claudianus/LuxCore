# Path Guiding — online path-space importance sampling

Status: implemented (M1 CPU + M2/M2b/M2c GPU). Learns a per-region incident
radiance field during the render and samples it to reduce variance.

## What and why

In indirectly-lit scenes the useful light arrives from a few directions, but
uniform BSDF sampling wastes most rays. Path guiding builds an approximate
incident-radiance field (a spatial-directional mixture / table) **during**
rendering and samples directions proportional to it — the renderer "learns
where the light is". This is a leading practical variance-reduction method.

## References

- Vorba, Herholz, Křivánek, et al. **On-line Learning of Parametric Mixture
  Models for Light-Transport Simulation.** SIGGRAPH 2014. (Online learning.)
- Müller, Gross, Novák. **Practical Path Guiding for Efficient Light-Transport
  Simulation.** EGSR 2017. (The guiding method our implementation follows.)
- Jensen. **Importance Driven Path Tracing Using the Photon Map.** EGWR 1995.
  (Early guiding idea.)

## Implementation

| Stage | Commit | What |
|---|---|---|
| M1 | `CPU path guiding M1` | Uniform spatial grid, one-sample MIS (CPU). |
| M2b | `GPU path guiding M2b` | Frozen coarse-table sampling on GPU. |
| M2b-2 | `GPU online training` | Record buffers + drain + swap: GPU trains its table online. |
| M2c | `indirect training + adaptive mixture` | Indirect-radiance training data + depth/count adaptive mixture (CPU+GPU). |
| Rounds | `frozen training rounds` | Cache finalized after N training rounds. |

- Stability fixes: `Fix GPU+guiding exit-GC crash` (drain-path file dumps +
  redundant zeroing removed).
- Experiments (absolute smoothing, stratified ubin) are **default-off**.

## Test scenes / validation

- Indirect-heavy interiors (Cornell, classroom) — guiding reduces variance vs
  unguided at equal spp; verified no bias (guiding only changes sampling, not
  the estimator).

## Platforms

CPU (M1) and GPU (M2b/M2c: OpenCL + Metal).
