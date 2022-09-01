#ifndef __CONCRETE_FUNCTIONS_H
#define __CONCRETE_FUNCTIONS_H

#include <functional>

namespace Numerics {
class Function {
 public:
  virtual float operator()(float x) const = 0;
};

class LambdaBasedFunction : public Function {
 private:
  std::function<float(float)> f;

 public:
  LambdaBasedFunction() = delete;
  LambdaBasedFunction(const std::function<float(float)>& f);

  virtual float operator()(float x) const override;
};
}  // namespace Numerics

#endif