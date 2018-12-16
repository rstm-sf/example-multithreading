// Copyright (C) 2018 Rustam Sayfutdinov (rstm.sf@gmail.com)
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A_ PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef EXAMPLE_UTILS_H_
#define EXAMPLE_UTILS_H_

#include <vector>

namespace ex_m_thr {

template <typename T = float>
std::vector<T> generate_square_block_matrix(
    std::size_t nrows, std::size_t nblocks) {
  std::vector<T> mat(nrows * nrows);

  const std::size_t offset = nrows / nblocks;
  const std::size_t balance = nrows - offset * nblocks;

  std::size_t at, to;
  std::size_t start {0};
  for (std::size_t k = 0; k < nblocks; ++k) {
    at = start;
    to = start + offset;
    if (k < balance)
      ++to;

    for (std::size_t i = at; i < to; ++i)
      for (std::size_t j = at; j < to; ++j)
        if (i == j) {
          mat[i * nrows + j] = static_cast<T>(offset + 100);
        } else {
          mat[i * nrows + j] = static_cast<T>(1);
        }

    start = to;
  }

  return mat;
}

} // namespace ex_m_thr

#endif // EXAMPLE_UTILS_H_