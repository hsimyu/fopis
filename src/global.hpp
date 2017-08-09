#ifndef __TDPIC_GLOBAL_H_INCLUDED__
#define __TDPIC_GLOBAL_H_INCLUDED__
#include <iostream>
#include <map>
#include <math.h>
#include <string>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <boost/format.hpp>
#include <boost/multi_array.hpp>

using std::cout;
using std::endl;
using boost::format;

using tdArray = boost::multi_array<double, 3>;
using ObjectDefinedMap = boost::multi_array<bool, 3>;
using ObjectNodes = std::map< unsigned int, std::array<unsigned int, 3> >;

//! constants
constexpr double e = 1.60217733e-19;
constexpr double eps0 = 8.854187817e-12;
constexpr double mu0 = 4.0 * M_PI * 1.0e-7;
constexpr double me = 9.1093818872e-31;
constexpr double c = 2.99792458e8;
constexpr int MAX_PARTICLE_NUM = 100000;

//! Enumerators
enum class AXIS {x, y, z};
enum class AXIS_SIDE {low, up};

#endif
