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
#include <initializer_list>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace ex_m_thr {

enum class Method {
  GaussSeidel,
  SOR
};

template <typename T = float>
class LinearSystem {
public:
  LinearSystem<T>();
  virtual ~LinearSystem<T>();
  LinearSystem<T>(const LinearSystem<T>&);
  LinearSystem<T>(LinearSystem<T>&&);
  LinearSystem<T>& operator=(const LinearSystem<T>&);
  LinearSystem<T>& operator=(LinearSystem<T>&&);

  LinearSystem<T>(
    std::size_t max_steps,
    T accuracy,
    std::size_t nrows,
    std::initializer_list<T> A,
    std::initializer_list<T> rhs);

  std::vector<T> solution() const;

  std::size_t nsteps() const;
  std::vector<T> r_residual_norms() const;

  virtual void solve(Method method = Method::GaussSeidel);

protected:
  std::vector<T> step_solution_gauss_seidel();
  std::vector<T> step_solution_sor(T w = 0.5);
  bool is_convergence();

  const std::size_t max_steps_;
  const T accuracy_;

  std::vector<T> r_residual_norms_;

  const std::size_t nrows_;
  const std::size_t ncols_;

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
  r_residual_norms_.reserve(max_steps_);
};

template <typename T>
LinearSystem<T>::LinearSystem(
    std::size_t max_steps,
    T accuracy,
    std::size_t nrows,
    std::initializer_list<T> A,
    std::initializer_list<T> rhs)
  : max_steps_(max_steps),
    accuracy_(accuracy),
    nrows_(nrows),
    ncols_(nrows),
    A_(A_),
    lhs_(nrows),
    rhs_(rhs) {
  r_residual_norms_.reserve(max_steps_);
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
std::size_t LinearSystem<T>::nsteps() const { r_residual_norms_.size(); }

template <typename T>
std::vector<T> LinearSystem<T>::r_residual_norms() const {
  return r_residual_norms_;
}

template <typename T>
void LinearSystem<T>::solve(Method method) {
  for (std::size_t i = 0; i < max_steps_; ++i) {
    switch (method) {
      case Method::GaussSeidel: lhs_ = step_solution_gauss_seidel(); break;
      case Method::SOR:         lhs_ = step_solution_sor(); break;
      default:
        throw std::runtime_error("LinearSystem::Solve: Undefined method!");
    }

    if (is_convergence())
      break;
  }
}

template <typename T>
std::vector<T> LinearSystem<T>::step_solution_gauss_seidel() {
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
std::vector<T> LinearSystem<T>::step_solution_sor(T w) {
  std::vector<T> lhs_new(rhs_);

  for (std::size_t i = 0; i < nrows_; ++i) {
    for (std::size_t j = 0; j < i; ++j)
      lhs_new[i] -= A_[i * nrows_ + j] * lhs_new[j];

    for (std::size_t j = i + 1; j < ncols_; ++j)
      lhs_new[i] -= A_[i * nrows_ + j] * lhs_[j];

    // SOR
    lhs_new[i] = lhs_[i] + w * (lhs_new[i] / A_[i * nrows_ + i] - lhs_[i]);
  }

  return lhs_new;
}

template <typename T>
bool LinearSystem<T>::is_convergence() {
  T rr {0.0};
  T xx {0.0};
  for (std::size_t i = 0; i < nrows_; ++i) {
    T residual {0.0};
    // ncols == nrows
    for (std::size_t j = 0; j < nrows_; ++j)
      residual += A_[i * nrows_ + j] * lhs_[j];
    residual -= rhs_[i];

    rr += residual * residual;
    xx += lhs_[i] * lhs_[i];
  }

  T r_residual_norm = std::sqrt(rr / xx);
  r_residual_norms_.push_back(r_residual_norm);

  return r_residual_norm <= accuracy_;
}

} // namespace ex_m_thr

#endif // EXAMPLE_LINEAR_SYSTEM_H_