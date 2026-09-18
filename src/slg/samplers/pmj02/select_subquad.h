/*
 * Copyright (C) Andrew Helmer 2020.
 * Licensed under MIT Open-Source License.
 *
 * Subquadrant selection for the PMJ02 generator. Only the SwapXOrY variant
 * is needed for PMJ(0,2); the rest of pmj-cpp's select_subquad is omitted.
 * Vendored from pmj-cpp (https://github.com/Andrew-Helmer/pmj-cpp)
 * for LuxCoreRender (include paths flattened).
 */
#ifndef LUXCORE_SELECT_SUBQUAD_H_
#define LUXCORE_SELECT_SUBQUAD_H_

#include <utility>
#include <vector>

#include "pmj_util.h"

namespace pmj {
  typedef std::vector<std::pair<int, int>> (*subquad_fn)(
      const Point samples[], const int dim);

  /*
   * This will randomly choose once to swap X or swap Y, and will always swap
   * X or Y for all subquadrants. For PMJ02, this ensures that the next set of
   * samples are themselves a (0,2) sequence.
   *
   * Credit goes to Simon Brown for discovering this method with his Rust
   * implementation: https://github.com/sjb3d/pmj
   */
  std::vector<std::pair<int, int>> GetSubQuadrantsSwapXOrY(
    const Point samples[],
    const int dim);
}  // namespace pmj

#endif  // LUXCORE_SELECT_SUBQUAD_H_
