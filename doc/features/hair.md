# Hair — Chiang BSDF + strand vertex AOVs

Status: implemented (host + device). `hairmat` material and per-vertex strand
AOVs that feed Cycles "Hair Info" outputs.

## What and why

Two pieces make production hair shading work:

1. **A hair BSDF.** A fibre BCSDF with longitudinal roughness and a
   Melanin/eumelanin-pheomelanin absorption model gives the characteristic
   glints + colour of real hair — far better than a diffuse/glossy
   approximation.
2. **Per-strand varying attributes.** Cycles' Hair Info node needs `strand u`
   (parametric coordinate along the fibre) and a per-strand random value.
   These are stored as **vertex AOVs** on the tessellated strand mesh and
   readable in shading.

## References

- Chiang, Bitterli, Tappan, Burley. **A Practical and Controllable Hair and Fur
  Model for Production Path Tracing.** Computer Graphics Forum 35(2), 2016.
  (The `hairmat` model — longitudinal/azimuthal roughness + melanin.)
- Huang, Weidlich, et al. **A Microfacet-Based Hair Scattering Model.** EGSR
  2022. (Newer model, listed as future work.)

## Implementation

- `hairmat` material — Chiang model parameters (melanin, roughness, shift,
  ior). Host (CPU) and device (OpenCL/Metal) BSDF.
- `StrandsShape` writes `strand-u` and `strand-random` as vertex AOVs; Cycles
  Hair Info → these AOVs via the node reader.
- `cyHairFile` out-of-bounds read fixed.
- Hair geometry is tessellated to triangles (ribbon/solid modes); native Metal
  curve primitives are a roadmap item (E7) — current cost is tessellation
  overhead, not a wrong default.

## Test scenes / validation

- `scenes/strands/`, `scenes/cornell/hairmat-test.scn`, `scenes/strands/strandu-test.scn`.
- Renders on CPU and Metal GPU.

## Platforms

CPU, OpenCL GPU, Metal GPU.
