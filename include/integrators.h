#ifndef __CONCRETE_INTEGRATORS_H
#define __CONCRETE_INTEGRATORS_H

#include "functions.h"

namespace Numerics {
class Integrator {
 public:
  virtual float integrate(const Function& f, float from, float to,
                          float eps) const = 0;
};

class SimpsonMethodIntegrator : public Integrator {
 public:
  SimpsonMethodIntegrator();
  virtual float integrate(const Function& f, float from, float to,
                          float eps) const override;
  float integrate_step(const Function& f, float from, float to, uint32_t& n,
                       std::vector<float>& computed_points) const;
};
}  // namespace Numerics

#endif