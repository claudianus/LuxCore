/*
 * Copyright (C) Andrew Helmer 2020.
 * Licensed under MIT Open-Source License.
 *
 * SwapXOrY subquadrant selection only; see select_subquad.h.
 * Vendored from pmj-cpp (https://github.com/Andrew-Helmer/pmj-cpp)
 * for LuxCoreRender (include paths flattened).
 */
#include "select_subquad.h"

#include <utility>
#include <vector>

#include "pmj_util.h"

namespace pmj {

std::vector<std::pair<int, int>> GetSubQuadrantsSwapXOrY(
    const Point samples[],
    const int dim) {
  const int quad_dim = dim / 2;
  const int n = quad_dim*quad_dim;

  std::vector<std::pair<int, int>> choices(n);

  const bool swap_x = UniformRand() < 0.5;

  for (int i = 0; i < n; i++) {
    const Point& sample = samples[i];
    int x_pos = static_cast<int>(sample.x * dim);
    int y_pos = static_cast<int>(sample.y * dim);

    if (swap_x) x_pos = x_pos ^ 1;
    else
      y_pos = y_pos ^ 1;

    choices[i] = {x_pos, y_pos};
  }

  return choices;
}

}  // namespace pmj
