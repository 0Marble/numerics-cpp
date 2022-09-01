#ifndef __ERROR_H
#define __ERROR_H

#include <iostream>

namespace Numerics {
enum class Error : uint8_t { BadRange };

std::ostream& operator<<(std::ostream& out, const Error& e);

}  // namespace Numerics

#endif