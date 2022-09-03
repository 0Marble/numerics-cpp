
#include "numerics.h"

namespace Numerics {

std::ostream& operator<<(std::ostream& out, const Error& e) {
  switch (e) {
    case Error::BadRange:
      out << "Error: BadRange";
      break;
    default:
      break;
  }
  return out;
}

LambdaBasedFunction::LambdaBasedFunction(const std::function<float(float)>& f)
    : f(f) {}
float LambdaBasedFunction::operator()(float x) const { return f(x); }

std::variant<float, Error> SecantMethodRootFinder::root(const Function& f,
                                                        const Function& g,
                                                        float from, float to,
                                                        float eps) const {
  auto [min, max] = std::minmax(from, to);
  float a = min, b = max;
  float f_a = f(a) - g(a), f_b = f(b) - g(b);
  if (f_a == 0.0) return a;
  if (f_b == 0.0) return b;

  while (true) {
    if (f_a * f_b > 0.0) return Error::BadRange;

    float c = (a * f_b - b * f_a) / (f_b - f_a);
    float f_c = f(c) - g(c);
    if (f_c == 0.0) return c;

    float midpoint = (b - a) / 2.0 + a;

    bool first_derivative_positive = f_a > 0.0f;
    bool second_derivative_positive =
        f(midpoint) - g(midpoint) > (f_a - f_b) / 2.0f;

    if (first_derivative_positive xor second_derivative_positive) {
      a = c;
      f_a = f_c;

      if (f_c * (f(c + eps) - g(c + eps)) < 0) {
        return c;
      }
    } else {
      b = c;
      f_b = f_c;

      if (f_c * (f(c - eps) - g(c - eps)) < 0) {
        return c;
      }
    }
  }
}

SimpsonMethodIntegrator::SimpsonMethodIntegrator() {}

float SimpsonMethodIntegrator::integrate_step(
    const Function& f, float from, float to, uint32_t& n,
    std::vector<float>& computed_points) const {
  if (computed_points.size() < 3) {
    float step = (to - from) / 3.0;
    computed_points.reserve(3);

    computed_points.push_back(f(from));
    computed_points.push_back(f(to));
    computed_points.push_back(f((to + from) / 2.0));

    n = 2;
    return (2.0 * computed_points[0] + 2.0 * computed_points[1] +
            4.0 * computed_points[2]) *
           (to - from) / 6.0;
  }

  float sum = 0.0;
  float step = (to - from) / float(n);
  computed_points.reserve(n + n + 1);
  for (uint32_t i = 0; i < n; i++) {
    float x = step * (float(i) + 0.5) + from;
    float y = f(x);
    computed_points.push_back(y);
    sum += 4.0 * y;
  }

  for (uint32_t i = 1; i < n; i++) {
    sum += 2.0 * computed_points[i];
  }

  sum += computed_points[0] + computed_points[n];
  n += n;
  float new_step = (to - from) / float(n);
  return sum * new_step / 3.0;
}

float SimpsonMethodIntegrator::integrate(const Function& f, float from,
                                         float to, float eps) const {
  std::vector<float> computed_points;
  uint32_t n = 0;
  float prev_int = integrate_step(f, from, to, n, computed_points);
  while (true) {
    float cur_int = integrate_step(f, from, to, n, computed_points);
    if (std::abs(prev_int - cur_int) < eps * 15.0) {
      return cur_int;
    }
    prev_int = cur_int;
  }
}

}  // namespace Numerics
