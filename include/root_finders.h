#ifndef __CONCRETE_ROOT_FINDERS_H
#define __CONCRETE_ROOT_FINDERS_H

#include <variant>

#include "error.h"
#include "functions.h"

namespace Numerics {
class RootFinder {
 public:
  virtual std::variant<float, Error> root(const Function& f, const Function& g,
                                          float from, float to,
                                          float eps) const = 0;
};

class SecantMethodRootFinder : public RootFinder {
 public:
  virtual std::variant<float, Error> root(const Function& f, const Function& g,
                                          float from, float to,
                                          float eps) const override;
};

}  // namespace Numerics

#endif