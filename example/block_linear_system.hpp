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

#ifndef EXAMPLE_BLOCK_LINEAR_SYSTEM_H_
#define EXAMPLE_BLOCK_LINEAR_SYSTEM_H_

#include <cmath>
#include <initializer_list>
#include <stdexcept>
#include <vector>

#include "block_jacobi.hpp"
#include "linear_system.hpp"

namespace ex_m_thr {

template <typename T = float>
class BlockLinearSystem : public LinearSystem<T> {
public:
  BlockLinearSystem<T>();
  virtual ~BlockLinearSystem<T>();
  BlockLinearSystem<T>(const BlockLinearSystem<T>&);
  BlockLinearSystem<T>(BlockLinearSystem<T>&&);
  BlockLinearSystem<T>& operator=(const BlockLinearSystem<T>&);
  BlockLinearSystem<T>& operator=(BlockLinearSystem<T>&&);

  BlockLinearSystem<T>(
    std::size_t nblocks,
    std::size_t max_steps,
    T accuracy,
    std::size_t nrows,
    const std::vector<T>& A,
    const std::vector<T>& rhs);

  BlockLinearSystem<T>(
    std::size_t nblocks,
    std::size_t max_steps,
    T accuracy,
    std::size_t nrows,
    std::initializer_list<T> A,
    std::initializer_list<T> rhs);

private:
  virtual std::vector<T> step_solution_gauss_seidel() override;

  BlockJacobi<T> preconditioner_;
};

template <typename T>
BlockLinearSystem<T>::BlockLinearSystem()
  : LinearSystem<T>(), preconditioner_(2, this->nrows_, this->A_) {};

template <typename T>
BlockLinearSystem<T>::~BlockLinearSystem() = default;

template <typename T>
BlockLinearSystem<T>::BlockLinearSystem(const BlockLinearSystem<T>&) = default;

template <typename T>
BlockLinearSystem<T>::BlockLinearSystem(BlockLinearSystem<T>&&) = default;

template <typename T>
BlockLinearSystem<T>& BlockLinearSystem<T>::operator=(const BlockLinearSystem<T>&) = default;

template <typename T>
BlockLinearSystem<T>& BlockLinearSystem<T>::operator=(BlockLinearSystem<T>&&) = default;

template <typename T>
BlockLinearSystem<T>::BlockLinearSystem(
    std::size_t nblocks,
    std::size_t max_steps,
    T accuracy,
    std::size_t nrows,
    const std::vector<T>& A,
    const std::vector<T>& rhs)
  : LinearSystem<T>(max_steps, accuracy, nrows, A, rhs),
    preconditioner_(nblocks, this->nrows_, this->A_) {};

template <typename T>
BlockLinearSystem<T>::BlockLinearSystem(
    std::size_t nblocks,
    std::size_t max_steps,
    T accuracy,
    std::size_t nrows,
    std::initializer_list<T> A,
    std::initializer_list<T> rhs)
  : LinearSystem<T>(max_steps, accuracy, nrows, A, rhs),
    preconditioner_(nblocks, this->nrows_, this->A_) {};

template <typename T>
std::vector<T> BlockLinearSystem<T>::step_solution_gauss_seidel() {
  return preconditioner_.step_solution_gauss_seidel(this->lhs_, this->rhs_);
}

} // namespace ex_m_thr

#endif // EXAMPLE_BLOCK_LINEAR_SYSTEM_H_