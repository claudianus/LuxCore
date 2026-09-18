/*
 * Copyright (C) Andrew Helmer 2020.
 * Licensed under MIT Open-Source License.
 *
 * Minimal utility subset vendored from pmj-cpp
 * (https://github.com/Andrew-Helmer/pmj-cpp) for LuxCoreRender:
 * only what the non-best-candidate PMJ02 generator needs
 * (Point, UniformRand, UniformInt, plus a seed hook so renders stay
 * reproducible instead of pulling from std::random_device).
 */
#ifndef LUXCORE_PMJ_UTIL_H_
#define LUXCORE_PMJ_UTIL_H_

namespace pmj {

typedef struct {
  double x;
  double y;
} Point;

// Fixed-seed hook: call once before GenerateSamples for reproducible sets.
// (Upstream uses a thread_local random_device; LuxCore wants determinism.)
void SetSeed(unsigned int seed);

// Gets a random double between any two numbers. Thread-safe after SetSeed.
double UniformRand(double min = 0.0, double max = 1.0);
// Generates a random int in the given range. Thread-safe after SetSeed.
int UniformInt(int min, int max);

}  // namespace pmj

#endif  // LUXCORE_PMJ_UTIL_H_
