#include "position.hpp"
#include "particle.hpp"

Position::Position(const Particle& p){
    this->setXYZ(p.x, p.y, p.z);
}

Velocity::Velocity(const Particle& p){
    this->set(p.vx, p.vy, p.vz);
}

std::ostream& operator<<(std::ostream& ost, const Position& pos) {
    ost << format("ipos = %d, %d, %d\n") % pos.i % pos.j % pos.k;
    ost << format("dpos = %s, %s, %s\n") % pos.x % pos.y % pos.z;
    ost << format("dx1 = %e, dx2 = %e\n") % pos.dx1 % pos.dx2;
    ost << format("dy1 = %e, dy2 = %e\n") % pos.dy1 % pos.dy2;
    ost << format("dz1 = %e, dz2 = %e") % pos.dz1 % pos.dz2 << endl;
    return ost;
}

std::ostream& operator<<(std::ostream& ost, const Velocity& vel) {
    ost << format("velocity = %s, %s, %s\n") % vel.vx % vel.vy % vel.vz << endl;
    return ost;
}