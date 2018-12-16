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
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <cmath>

#include <gtest/gtest.h>

#include "linear_system.hpp"

class LinearSystemTests : public ::testing::Test {};

TEST_F(LinearSystemTests, test_default) {
  ex_m_thr::LinearSystem ls;
  std::vector<float> exact_solution({1.0, 2.0, 3.0});

  ls.solve();

  float dd {0.0};
  for (std::size_t i = 0; i < exact_solution.size(); ++i) {
    float d = ls.solution()[i] - exact_solution[i];
    dd += d * d;
  }

  EXPECT_TRUE(std::sqrt(dd) < 1.0e-6);
}

TEST_F(LinearSystemTests, ex) {
  std::size_t max_steps{100};
  float accuracy {1.0e-6};
  std::size_t nrows {4};
  std::vector<float> A({
    10.0, -1.0,  0.0,  0.0,
    -1.0, 11.0,  0.0,  0.0,
     0.0,  0.0, 10.0, -1.0,
     0.0,  0.0, -1.0,  8.0
  });
  std::vector<float> rhs({6.0, 25.0, -11.0, 15.0});

  ex_m_thr::LinearSystem ls(max_steps, accuracy, nrows, A, rhs);

  std::vector<float> exact_solution({0.83486235, 2.3486238, -0.92405063, 1.7594937});

  ls.solve();

  float dd {0.0};
  for (std::size_t i = 0; i < exact_solution.size(); ++i) {
    float d = ls.solution()[i] - exact_solution[i];
    dd += d * d;
  }

  EXPECT_TRUE(std::sqrt(dd) < 1.0e-5);
}