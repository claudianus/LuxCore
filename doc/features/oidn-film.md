# Film hardware image pipeline + OIDN denoising

Status: implemented (partial). The film tone-map / output image pipeline runs
on the render device (incl. Metal); Intel Open Image Denoise is integrated but
runs on CPU on Apple (see gap).

## What and why

Two hardware-acceleration items around the film/output stage:

1. **Hardware film pipeline** — the image pipeline (tonemap, gamma, plugin
   chain) executes as GPU kernels on the same device as the render, instead of
   downloading the film to the CPU each output.
2. **OIDN** — Intel's machine-learning denoiser removes residual Monte Carlo
   noise; production renders rely on it for usable output at moderate spp.

## Implementation

- `Film hardware image pipeline: run on the render engine's device` — the
  pipeline kernels are compiled and run on the render device (OpenCL or
  Metal), so film output stays on-GPU.
- `OIDN: use the Metal device on Apple Silicon, fall back to CPU` — the code
  requests `oidn::DeviceType::Metal` on Apple.

### Known gap (honest status)

The vendored OIDN 2.5.0 ships only the **CPU device binary** on this platform —
`device_metal` is not built, so the request silently falls back to CPU
denoising. Rebuilding OIDN with its Metal device is roadmap item E1. On other
platforms (CUDA/ROCm/SYCL) OIDN GPU may already be available from stock builds.

## References

- Intel Open Image Denoise — https://www.openimagedenoise.org (ML denoiser).

## Platforms

Film HW pipeline: OpenCL + Metal. OIDN: CPU here; GPU where a device build
exists.
