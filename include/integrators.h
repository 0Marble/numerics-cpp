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
 private:
  uint32_t start_node_count;

 public:
  SimpsonMethodIntegrator() = delete;
  SimpsonMethodIntegrator(uint32_t start_node_count);
  virtual float integrate(const Function& f, float from, float to,
                          float eps) const override;
  float integrate_step(const Function& f, float from, float to, uint32_t& n,
                       std::vector<float>& computed_points) const;
};
}  // namespace Numerics

#endif