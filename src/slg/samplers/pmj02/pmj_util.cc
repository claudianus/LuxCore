/*
 * Copyright (C) Andrew Helmer 2020.
 * Licensed under MIT Open-Source License.
 *
 * Minimal utility subset vendored from pmj-cpp
 * (https://github.com/Andrew-Helmer/pmj-cpp) for LuxCoreRender.
 */
#include "pmj_util.h"

#include <random>

namespace pmj {

namespace {
// Seeded once per generated set (see SetSeed); upstream used random_device.
thread_local static std::default_random_engine gen(1u);
}  // namespace

void SetSeed(unsigned int seed) {
  gen.seed(seed);
}

double UniformRand(double min, double max) {
  thread_local static std::uniform_real_distribution<double> uniform;

  std::uniform_real_distribution<double>::param_type param(min, max);

  return uniform(gen, param);
}

int UniformInt(int min, int max) {
  thread_local static std::uniform_int_distribution<int> uniform;

  std::uniform_int_distribution<int>::param_type param(min, max);

  return uniform(gen, param);
}

}  // namespace pmj
