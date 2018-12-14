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

#ifndef EXAMPLE_LINEAR_SYSTEM_H_
#define EXAMPLE_LINEAR_SYSTEM_H_

#include <cmath>
#include <numeric>
#include <vector>

namespace ex_m_thr {

template <typename T = float>
class LinearSystem {
public:
  LinearSystem<T>();
  ~LinearSystem<T>();
  LinearSystem<T>(const LinearSystem<T>&);
  LinearSystem<T>(LinearSystem<T>&&);
  LinearSystem<T>& operator=(const LinearSystem<T>&);
  LinearSystem<T>& operator=(LinearSystem<T>&&);

  std::vector<T> solution() const;

  std::size_t nsteps() const;

  void solve();

private:
  std::vector<T> step_gauss_seidel();
  std::vector<T> step_sor(T w);
  bool is_convergence(std::vector<T>& lhs_new);

  const std::size_t max_steps_;
  const T accuracy_;

  std::vector<T> step_norms_;

  std::size_t nrows_;
  std::size_t ncols_;

  std::vector<T> A_;
  std::vector<T> lhs_;
  std::vector<T> rhs_;
};

// Пример из https://s-mat-pcs.oulu.fi/~mpa/matreng/eem5_4-1.htm
template <typename T>
LinearSystem<T>::LinearSystem()
  : max_steps_(100),
    accuracy_(1.0e-6),
    nrows_(3),
    ncols_(nrows_),
    A_({4.0,  1.0, -1.0, 2.0,  7.0,  1.0, 1.0, -3.0, 12.0}),
    lhs_(nrows_),
    rhs_({3.0, 19.0, 31.0}) {
  step_norms_.reserve(max_steps_);
};

template <typename T>
LinearSystem<T>::~LinearSystem() = default;

template <typename T>
LinearSystem<T>::LinearSystem(const LinearSystem<T>&) = default;

template <typename T>
LinearSystem<T>::LinearSystem(LinearSystem<T>&&) = default;

template <typename T>
LinearSystem<T>& LinearSystem<T>::operator=(const LinearSystem<T>&) = default;

template <typename T>
LinearSystem<T>& LinearSystem<T>::operator=(LinearSystem<T>&&) = default;

template <typename T>
std::vector<T> LinearSystem<T>::solution() const { return lhs_; }

template <typename T>
std::size_t LinearSystem<T>::nsteps() const { step_norms_.size(); }

template <typename T>
void LinearSystem<T>::solve() {
  std::vector<T> lhs_new;

  for (std::size_t i = 0; i < max_steps_; ++i) {
    lhs_new = step_gauss_seidel();
    bool is_stop = is_convergence(lhs_new);
    lhs_ = lhs_new;

    if (is_stop)
      break;
  }
}

template <typename T>
std::vector<T> LinearSystem<T>::step_gauss_seidel() {
  std::vector<T> lhs_new(rhs_);

  for (std::size_t i = 0; i < nrows_; ++i) {
    for (std::size_t j = 0; j < i; ++j)
      lhs_new[i] -= A_[i * nrows_ + j] * lhs_new[j];

    for (std::size_t j = i + 1; j < ncols_; ++j)
      lhs_new[i] -= A_[i * nrows_ + j] * lhs_[j];

    lhs_new[i] /= A_[i * nrows_ + i];
  }

  return lhs_new;
}

template <typename T>
std::vector<T> LinearSystem<T>::step_sor(T w) {
  std::vector<T> lhs_new(rhs_);

  for (std::size_t i = 0; i < nrows_; ++i) {
    for (std::size_t j = 0; j < i; ++j)
      lhs_new[i] -= A_[i * nrows_ + j] * lhs_new[j];

    for (std::size_t j = i + 1; j < ncols_; ++j)
      lhs_new[i] -= A_[i * nrows_ + j] * lhs_[j];

    lhs_new[i] /= A_[i * nrows_ + i];
  }

  for (std::size_t i = 0; i < nrows_; ++i)
    lhs_new[i] = (static_cast<T>(1.0) - w) * lhs_[i] + w * lhs_new[i];

  return lhs_new;
}

template <typename T>
bool LinearSystem<T>::is_convergence(std::vector<T>& lhs_new) {
  T dd {0.0};
  for (std::size_t i = 0; i < nrows_; ++i) {
    T d = lhs_new[i] - lhs_[i];
    dd += d * d;
  }
  T xx = std::inner_product(
    std::begin(lhs_), std::end(lhs_), std::begin(lhs_), 0.0);
  T norm = std::sqrt(dd / xx);

  step_norms_.push_back(norm);

  return norm <= accuracy_;
}

} // namespace ex_m_thr

#endif // EXAMPLE_LINEAR_SYSTEM_H_