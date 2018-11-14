#ifndef SMALL_HELPERS_HH
#define SMALL_HELPERS_HH

#include <string>
#include <sstream>
#include <cstdio>
#include <iostream>

using std::cerr;
using std::endl;

template< class a_type >
std::string to_string(a_type zahl) {
  std::stringstream fluss;
  fluss << zahl;
  return fluss.str();
};

// std::string to_string(double zahl) {
//   std::stringstream fluss;
//   fluss << zahl;
//   return fluss.str();
// };

FILE* fopen_with_check(std::string fname, const char* attr);

#endif
