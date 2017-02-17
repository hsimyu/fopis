#include "position.hpp"
#include "particle.hpp"

Position::Position(const Particle& p){
    this->setXYZ(p.x, p.y, p.z);
}

Velocity::Velocity(const Particle& p){
    this->set(p.vx, p.vy, p.vz);
}

std::ostream& operator<<(std::ostream& ost, Position const& pos) {
    ost << format("ipos = %d, %d, %d") % pos.i % pos.j % pos.k << endl;
    ost << format("dpos = %d, %d, %d") % pos.x % pos.y % pos.z << endl;
    return ost;
}
