#ifndef __TDPIC_SPACECRAFT_H_INCLUDED__
#define __TDPIC_SPACECRAFT_H_INCLUDED__
#include "global.hpp"

//! @class: Spacecraft
class Spacecraft {
private:
    static unsigned int num_of_spacecraft;
    std::string name;
    objectArray object_map;
    void construct(const size_t, const size_t, const size_t);

public:
    Spacecraft(const size_t nx, const size_t ny, const size_t nz) :
        name("spacecraft" + std::to_string(num_of_spacecraft)), object_map(boost::extents[0][0][0]) {
        construct(nx, ny, nz);
    }

    Spacecraft(const size_t nx, const size_t ny, const size_t nz, const std::string _name) : name(_name), object_map(boost::extents[0][0][0]) {
        construct(nx, ny, nz);
    }

    void setName(const std::string _name) { name = _name; }
    std::string getName() const { return name; }
    friend std::ostream& operator<<(std::ostream&, const Spacecraft&);
};
#endif