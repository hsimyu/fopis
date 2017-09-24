#ifndef __TDPIC_GLOBAL_H_INCLUDED__
#define __TDPIC_GLOBAL_H_INCLUDED__
#include <iostream>
#include <map>
#define _USE_MATH_DEFINES
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

class Particle;
using ParticleArray = std::vector<std::vector<Particle>>;

//! constants
constexpr double e = 1.60217733e-19;
constexpr double eps0 = 8.854187817e-12;
constexpr double mu0 = 4.0 * M_PI * 1.0e-7;
constexpr double me = 9.1093818872e-31;
constexpr double c = 2.99792458e8;
constexpr int MAX_PARTICLE_NUM = 100000000;

//! Enumerators
enum class AXIS {x, y, z};
enum class AXIS_SIDE {low, up};

//! 物体情報用の構造体
struct ObjectInfo_t {
    std::string name;
    std::string file_name;
    std::string surface_type;
    unsigned int history_width;
    double potential_fix;
    std::vector<std::string> emit_particle_names;
    std::map<int, std::string> materials;
};

#endif
