#include "global.hpp"
#include "particle.hpp"
#include "environment.hpp"

void Particle::setPosition(Position const& pos){
    x = pos.x;
    y = pos.y;
    z = pos.z;
}

double Particle::getEnergy(void) const {
    return 0.5 * (vx*vx + vy*vy + vz*vz) * Environment::ptype[typeId].getMass();
}

// util
std::ostream& operator<<(std::ostream& ost, Particle const& p){
    ost << "[  Particle  ]" << endl;
    ost << "         id: " << p.typeId << endl;
    ost << "   position: " << format("%f, %f, %f") % p.x % p.y % p.z << endl;
    ost << "   velocity: " << format("%f, %f, %f") % p.vx % p.vy % p.vz << endl;
    ost << "    isValid: " << format("%s") % p.isValid << endl;
    return ost;
}
