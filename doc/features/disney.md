# Disney transmission — integrated BTDF in the principled shader

Status: implemented. The Disney/principled material gained a real dielectric
transmission lobe, replacing a mix-with-a-transmissive-material hack.

## What and why

The Disney "principled" BRDF is the de-facto production uber-shader (Blender's
Principled BSDF is based on it). Its `transmission` parameter blends between
the opaque lobe stack and a dielectric BTDF (refraction with roughness /
thin-film). Previously LuxCore's Disney material had no integrated
transmission, so the adapter faked it by mixing materials — fragile and
slow. Now the transmission lobe is a first-class part of the BSDF.

## References

- Burley. **Extending the Disney BRDF.** Journal of Computer Graphics
  Techniques 4(1), 2015. (The full principled model incl. transmission.)
- Burley & Walt Disney Animation Studios. **Physically-Based Shading at
  Disney.** SIGGRAPH 2012 course notes. (Original model.)
- Walter, Marschner, Li, Torrance. **Microfacet Models for Refraction through
  Rough Surfaces.** EGSR 2007. (The rough-dielectric BTDF used.)

## Implementation

- `disney` material `transmission` parameter → rough-dielectric BTDF weight;
  thin-film interference supported. Host + device eval.
- The Blender adapter maps Principled v2 `transmission` directly (drops the
  old mix hack): `Principled v2: integrated Disney transmission + thin-film`.

## Test scenes / validation

- `scenes/cornell/cornell-disney.scn` — transmissive principled material.
- CPU / Metal GPU renders.

## Platforms

CPU, OpenCL GPU, Metal GPU.
