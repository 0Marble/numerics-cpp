#include <cmath>
#include <iostream>

#include "numerics.h"

using namespace Numerics;

int main() {
  auto f = LambdaBasedFunction([](float x) { return std::exp(x); });

  float answer = std::exp(5.0f) - std::exp(-5.0f), eps = 0.0001f;
  float computed_answer =
      SimpsonMethodIntegrator().integrate(f, -5.0f, 5.0f, eps);

  if (std::abs(answer - computed_answer) > 10.0f * eps) {
    std::cerr << "Expected " << answer << " +- " << eps << ", got "
              << computed_answer << "\n";
    return 1;
  }

  return 0;
}