# Spectral Rendering — hero-wavelength spectral transport

Status: implemented (CPU + GPU). True per-wavelength sampling instead of a
fixed RGB triple, enabling physically-correct dispersion.

## What and why

Ordinary renderers carry one RGB triple per path, which cannot reproduce
wavelength-dependent effects (prism dispersion, thin-film interference, real
material spectra). Hero-wavelength sampling carries a small set of wavelengths
anchored at a randomly sampled **hero wavelength**; each path evaluates the
scene's spectral response at its wavelengths and the result is projected back
to RGB through colour-matching functions. LuxCore samples 3 wavelengths per
path (hero + 2 rotated).

## References

- Wilkie, Nawaz, Droske, Weidlich, Hanika. **Hero Wavelength Spectral
  Sampling.** Computer Graphics Forum (EGSR) 33(4), 2014. (The hero-λ method.)
- Radziszewski, Boryczko, Sun. **An Improved Technique for Full Spectral
  Rendering.** Journal of WSCG 2009. (Companion-wavelength idea.)
- CIE 1931 colour-matching functions — used to integrate SPD -> XYZ -> RGB.

## Implementation

| Commit | What |
|---|---|
| `Spectral transport S1+S2 (CPU)` | Hero-λ state on the path; SPD evaluation at path wavelengths; CPU. |
| `A2: GPU hero-wavelength spectral transport` | Same transport on OpenCL/Metal PATHOCL engines. |
| `spectral GPU hardening` | NaN/OOB guards (added after a WindowServer watchdog panic). |
| `export spectral flag for RTPATHCPU` | Blender UI exposes the spectral option. |

- Spectral tables (CIE, Planck SPD, Smits RGB->SPD) are shared between CPU and
  the `.cl` kernel source (duplicated constant tables — the codebase has no
  `.h` include path in kernel compilation).
- Glass dispersion is the canonical visible result; emission/absorption media
  also evaluate at path wavelengths.

## Test scenes / validation

- `scenes/cornell/cornell-spectral*.scn`, prism dispersion.
- CPU/GPU output compared for parity.

## Platforms

CPU, OpenCL GPU, Metal GPU.
