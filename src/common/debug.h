#ifndef FUROO_COMMON_DEBUG_H
#define FUROO_COMMON_DEBUG_H

#include <cmath>
#include <iostream>
#include <sstream>

#ifdef FUROO_DEBUG
#define CERR(A) std::cerr << "[INFO]: " << A << std::endl
#else
#define CERR(A)                                                                \
  do {                                                                         \
  } while (0)
#endif

#define LOG_ERR(A) CERR(A)

#define LOG_INFO(A) CERR(A)

#define LOG_LOCATION "[" << __FILE__ << ", line " << __LINE__ << "]"

#define UNUSED_VARIABLE(x) ((void)x)

#define THROW(A, B)                                                            \
  {                                                                            \
    if (!(A)) {                                                                \
      std::ostringstream message;                                              \
      message << LOG_LOCATION;                                                 \
      std::cerr << (B) << std::endl;                                           \
      throw(message.str());                                                    \
    }                                                                          \
  }

#define NOT_TESTED()                                                           \
  { std::cout << "Code not tested! " << LOG_LOCATION << std::endl; }

#define NOT_IMPLEMENTED()                                                      \
  { std::cout << "Function not implemented! " << LOG_LOCATION << std::endl; }

#define ASSERT_FATAL(A)                                                        \
  {                                                                            \
    if (!(A)) {                                                                \
      std::cout << "Assertion failed in " << LOG_LOCATION << std::endl;        \
      exit(1);                                                                 \
    }                                                                          \
  }

#define ASSERT(A)                                                              \
  {                                                                            \
    if (!(A))                                                                  \
      std::cout << "Assertion failed in " << LOG_LOCATION << std::endl;        \
  }

#define CHECK_IN_BETWEEN(A, B, C)                                              \
  {                                                                            \
    if (!((A) >= (B) && (A) <= (C)))                                           \
      std::cout << "Assertion failed in " << LOG_LOCATION << std::endl;        \
  }

#define CHECK_FLOAT_EQUAL(A, B)                                                \
  {                                                                            \
    if (fabs((A) - (B)) < 1e-8)                                                \
      std::cout << LOG_LOCATION << " " << std::endl;                           \
  }

#define LOG std::cout << LOG_LOCATION << " "

#define PRINT(A) std::cout << A << std::endl;

#define DUMP_VECTOR(V)                                                         \
  {                                                                            \
    std::cout << "VECTOR in " << LOG_LOCATION << std::endl;                    \
    for (size_t i = 0; i < V.size(); ++i)                                      \
      std::cout << V[i] << " ";                                                \
    std::cout << std::endl;                                                    \
  }

#define DUMP_MATRIX(M)                                                         \
  {                                                                            \
    std::cout << "MATRIX in " << LOG_LOCATION << std::endl;                    \
    for (int i = 0; i < M.size(); ++i) {                                       \
      for (int j = 0; j < M[i].size(); ++j)                                    \
        std::cout << M[i][j] << " ";                                           \
    }                                                                          \
    std::cout << std::endl;                                                    \
  }

namespace furoo {

// Condition to stop processing
inline void concatenate(std::ostringstream &s) {}

template <typename H, typename... T>
void concatenate(std::ostringstream &s, H p, T... t) {
  s << p;
  concatenate(s, t...);
}

template <class... Args> std::string concat(const Args &... args) {
  std::ostringstream s;
  concatenate(s, args...);
  return s.str();
}

inline void printBits(unsigned int n) {
  for (int i = 31; i >= 0; i--)
    if ((1 << i) & n)
      std::cout << '1';
    else
      std::cout << '0';
}
} // namespace furoo
#endif
