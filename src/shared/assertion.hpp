/* Copyright (C) 2023  Enrico Degregori, Wilton Jaciel Loch
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef SRC_SHARED_ASSERTION_HPP_
#define SRC_SHARED_ASSERTION_HPP_

#include <iostream>

#ifndef NDEBUG
  #include <cassert>
  #define YAOP_ASSERT(condition)                                              \
  {                                                                           \
      if (!(condition)) {                                                     \
          std::cerr << "Assertion failed at " << __FILE__ << ":" << __LINE__; \
          std::cerr << " inside " << __FUNCTION__ << std::endl;               \
          std::cerr << "Condition: " << #condition;                           \
          abort();                                                            \
      }                                                                       \
  }
#else
  #define YAOP_ASSERT(condition) (condition)
#endif

#endif  // SRC_SHARED_ASSERTION_HPP_
