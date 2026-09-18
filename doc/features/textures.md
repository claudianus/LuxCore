# New / extended textures — blackbody (texturable) + whitenoise

## Blackbody texture — texturable temperature

Status: implemented (CPU + GPU + spectral). `blackbody.temperature` accepts a
texture, not only a constant.

**What/why.** Emission driven by a spatial field — fire / volume temperature,
procedural heat — needs the Planckian colour evaluated per shading point.
Previously `.temperature` was a compile-time constant, so a density grid's
temperature channel could not drive the colour. Now it can.

**How.** The RGB path uses a 160-entry temperature→RGB lookup table baked from
`TemperatureToWhitePoint` (the full Planck→XYZ integral is too costly per
shading point); the same normalized table is embedded in the kernel source.
The spectral path still evaluates the true Planckian SPD per wavelength
(`Spectral_BlackbodyEval` / `BlackbodySPD`), so spectral renders keep exact
physics rather than an RGB approximation.

- Constant `.temperature = <K>` still works and is byte-identical (exact
  integral, not the LUT).
- `blackbody.normalize = 0|1` — 1 gives normalized chromaticity, 0 the raw
  physical scale.
- **References:** Planck's law (blackbody SPD); CIE colour-matching functions.
  `scenes/cornell/bb-test.scn` demonstrates a density grid driving fire
  emission temperature.
- **Validation:** deterministic surface render is **bit-identical** CPU↔GPU
  (max|diff| = 0); volume fire render identical on CPU / OpenCL / Metal /
  spectral.

## Whitenoise texture — deterministic 3D-seed noise

Status: implemented (CPU + GPU). Cycles White Noise equivalent.

**What/why.** A hash of a 3D position seed to a deterministic [0,1) value /
colour — uncorrelated per-point noise used for randomization (sprinkle
variation, anti-repetition). Cycles parity.

**How.** `whitenoise.texture` = the vector seed (default `position`);
`whitenoise.seed` = an offset. The seed is hashed by bit-reinterpretation
(`as_uint` on GPU / union on CPU) + an integer avalanche mix so CPU and GPU
produce **identical** values even for negative / large world coordinates — a
plain arithmetic float→uint cast is UB for negatives and diverged between the
two (found and fixed during validation).

- `whitenoise.scn` float + colour outputs.
- **Validation:** CPU↔GPU identical distribution on a Cornell-box render.

## Platforms

Both textures: CPU, OpenCL GPU, Metal GPU (cl2msl-compatible kernel code).
