/*
 * Copyright (C) Andrew Helmer 2020.
 * Licensed under MIT Open-Source License.
 *
 * PMJ(0,2) sequence generator (Christensen et al. 2018, Pharr 2019 fast
 * generation). Only the plain PMJ02 entry point is vendored; the
 * blue-noise/balancing experiments from pmj-cpp are omitted.
 * Vendored from pmj-cpp (https://github.com/Andrew-Helmer/pmj-cpp)
 * for LuxCoreRender (include paths flattened, seed hook added).
 */
#ifndef LUXCORE_PMJ02_H_
#define LUXCORE_PMJ02_H_

#include <memory>

#include "pmj_util.h"

namespace pmj {

// Generates progressive multi-jittered (0,2) samples WITHOUT blue noise
// properties. Takes in a number of samples.
std::unique_ptr<Point[]> GetPMJ02Samples(const int num_samples);

}  // namespace pmj

#endif  // LUXCORE_PMJ02_H_
