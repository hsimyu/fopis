#include "global.hpp"
#include "position.hpp"
#include "particle.hpp"
#include "environment.hpp"
#include "utils.hpp"

void Particle::setPosition(Position const& pos){
    x = pos.x;
    y = pos.y;
    z = pos.z;
}

void Particle::setVelocity(Velocity const& v){
    vx = v.vx;
    vy = v.vy;
    vz = v.vz;
}

Position&& Particle::getPosition(void) const {
    return Position{x, y, z};
}

Position&& Particle::getOldPosition(void) const {
    return Position{x - vx, y - vy, z - vz};
}

Position&& Particle::getNewPosition(void) const {
    return Position{x + vx, y + vy, z + vz};
}

double Particle::getEnergy(void) const {
    return 0.5 * (vx*vx + vy*vy + vz*vz) * Environment::ptype[typeId].getMass();
}

void Particle::generateNewPosition(const double min_x, const double max_x, const double min_y, const double max_y, const double min_z, const double max_z)
{
    this->setPosition( Environment::ptype[typeId].generateNewPosition(min_x, max_x, min_y, max_y, min_z, max_z) );
}

void Particle::generateNewVelocity(void)
{
    this->setVelocity( Environment::ptype[typeId].generateNewVelocity() );
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
