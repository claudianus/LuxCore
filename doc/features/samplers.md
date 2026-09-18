# Samplers — PMJ02 + blue-noise dithered Sobol

Status: implemented. Two low-discrepancy sampler improvements on top of the
existing Sobol infrastructure.

## What and why

Convergence rate is bounded by how evenly samples cover the domain. Two
orthogonal upgrades were added:

1. **PMJ02** — a progressive multi-jittered sequence with a provably better
   discrepancy bound than classic jittering, especially at low sample counts.
2. **Blue-noise dithered Sobol** — scrambles the Sobol sequence with a
   screen-space blue-noise mask so residual integration error appears as
   high-frequency noise (perceptually uniform) rather than structured banding.

## References

- Christensen, Kensler, Kilpatrick. **Progressive Multi-Jittered Sample
  Sequences.** EGSR 2018. (PMJ02 / PMJ02bn.)
- Heitz, Belcour, Ostroumoukhov, Coeurjolly, Iehl. **A Low-Discrepancy Sampler
  that Distributes Monte Carlo Errors as a Blue Noise in Screen Space.**
  SIGGRAPH 2019. (Blue-noise error distribution.)
- Ahmed, Wonka. **Screen-Space Blue-Noise Diffusion of Monte Carlo Sampling
  Error via Hierarchical Ordering of Pixels.** SIGGRAPH Asia 2020.
- Sobol'. **On the Distribution of Points in a Cube and the Approximate
  Evaluation of Integrals.** 1967. (Base sequence.)

## Implementation

- `P1-4: PMJ02 sampler` — the PMJ02 sequence generator.
- `SOTA P0: blue-noise dithered Sobol` — Sobol + screen-space blue-noise mask;
  unbiased.
- Exposed in Blender (`Expose Sobol blue-noise dithering in the UI`).

## Test scenes / validation

- Convergence at equal spp vs the previous sampler; the dither mask must not
  bias the estimator (verified it only permutes error, not the mean). PMJ02
  statistical (chi-square) uniformity is a roadmap validation item.

## Platforms

CPU, OpenCL GPU, Metal GPU (sampler tables are device-resident).
