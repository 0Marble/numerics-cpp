#include <cmath>
#include <iostream>

#include "numerics.h"

using namespace Numerics;

int main() {
  auto f =
      LambdaBasedFunction([](float x) { return 1.0 + 4.0 / (x * x + 1.0f); });
  auto g = LambdaBasedFunction([](float x) { return x * x * x; });

  float computed_answer, answer = 1.344, eps = 0.001;
  auto res = SecantMethodRootFinder().root(f, g, 0.5, 2.0, eps);

  switch (res.index()) {
    case 0:
      computed_answer = std::get<0>(res);
      break;
    case 1:
      std::cerr << std::get<1>(res) << "\n";
      return 1;
    default:
      break;
  }

  if (std::abs(computed_answer - answer) > eps) {
    std::cerr << "Expected " << answer << " +- " << eps << ", got "
              << computed_answer << "\n";
    return 1;
  }
  return 0;
}