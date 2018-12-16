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

#ifndef EXAMPLE_BLOCK_JACOBI_H_
#define EXAMPLE_BLOCK_JACOBI_H_

#include <thread>
#include <vector>

namespace ex_m_thr {

template <typename T = float>
class BlockJacobi {
public:
  BlockJacobi<T>();
  ~BlockJacobi<T>();
  BlockJacobi<T>(const BlockJacobi<T>&);
  BlockJacobi<T>(BlockJacobi<T>&&);
  BlockJacobi<T>& operator=(const BlockJacobi<T>&);
  BlockJacobi<T>& operator=(BlockJacobi<T>&&);

  BlockJacobi<T>(std::size_t nbs, std::size_t nrows, const std::vector<T>& A);

  std::vector<T> step_solution_gauss_seidel(
    const std::vector<T>& lhs, const std::vector<T>& rhs);
  std::vector<T> times(const std::vector<T>& rhs) const;

private:
  void step_solution_gauss_seidel_thr(
    const std::vector<T>& lhs,
    const std::vector<T>& rhs,
    std::vector<T>& lhs_new,
    std::uint32_t thr_id);

  const std::size_t nrows_;
  std::vector<T> A_;

  const std::size_t nblocks_;
  std::vector<std::size_t> offsets_;
};

template <typename T>
BlockJacobi<T>::BlockJacobi() = default;

template <typename T>
BlockJacobi<T>::~BlockJacobi() = default;

template <typename T>
BlockJacobi<T>::BlockJacobi(const BlockJacobi<T>&) = default;

template <typename T>
BlockJacobi<T>::BlockJacobi(BlockJacobi<T>&&) = default;

template <typename T>
BlockJacobi<T>& BlockJacobi<T>::operator=(const BlockJacobi<T>&) = default;

template <typename T>
BlockJacobi<T>& BlockJacobi<T>::operator=(BlockJacobi<T>&&) = default;

template <typename T>
BlockJacobi<T>::BlockJacobi(
    std::size_t nblocks, std::size_t nrows, const std::vector<T>& A)
  : nblocks_(nblocks), nrows_(nrows), A_(A) {
  std::size_t offset = nrows_ / nblocks_;
  std::size_t balance = nrows_ - offset * nblocks_;

  std::size_t start {0};
  offsets_.reserve(nblocks_ + 1);
  offsets_.push_back(start);
  for (std::size_t i = 0; i < balance; ++i) {
    start += offset + 1;
    offsets_.push_back(start);
  }
  for (std::size_t i = balance; i < nblocks_; ++i) {
    start += offset;
    offsets_.push_back(start);
  }
}

template <typename T>
std::vector<T> BlockJacobi<T>::step_solution_gauss_seidel(
    const std::vector<T>& lhs,
    const std::vector<T>& rhs) {
  std::vector<T> lhs_new(rhs);

  std::vector<std::thread> threads;
  for (std::size_t k = 0; k < nblocks_; ++k) {
    std::thread thr([&lhs, &rhs, &lhs_new, k, this] {
      step_solution_gauss_seidel_thr(lhs, rhs, lhs_new, k);
    }); 
    threads.emplace_back(std::move(thr));
  }

  for(auto& thr : threads)
    thr.join();

  return lhs_new;
}

template <typename T>
std::vector<T> BlockJacobi<T>::times(const std::vector<T>& rhs) const {
  std::vector<T> result;
  result.reserve(rhs.size());
  for (std::size_t k = 0; k < nblocks_; ++k) {
    const std::size_t at = offsets_[k];
    const std::size_t to = offsets_[k + 1];

    for(std::size_t i = at; i < to; ++i) {
      T r {0.0};
      for (std::size_t j = at; j < to; ++j)
        r += A_[i * nrows_ + j] * rhs[j];
      result.push_back(r);
    }
  }

  return result;
}

template <typename T>
void BlockJacobi<T>::step_solution_gauss_seidel_thr(
    const std::vector<T>& lhs,
    const std::vector<T>& rhs,
    std::vector<T>& lhs_new,
    std::uint32_t thr_id) {
  const std::size_t at = offsets_[thr_id];
  const std::size_t to = offsets_[thr_id + 1];

  for (std::size_t i = at; i < to; ++i) {
    for (std::size_t j = at; j < i; ++j)
      lhs_new[i] -= A_[i * nrows_ + j] * lhs_new[j];

    for (std::size_t j = i + 1; j < to; ++j)
      lhs_new[i] -= A_[i * nrows_ + j] * lhs[j];

    lhs_new[i] /= A_[i * nrows_ + i];
  }
}

} // namespace ex_m_thr

#endif // EXAMPLE_BLOCK_JACOBI_H_