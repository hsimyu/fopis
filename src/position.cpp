#include "position.hpp"
#include "particle.hpp"

Position::Position(const Particle& p){
    this->setXYZ(p.x, p.y, p.z);
}

Velocity::Velocity(const Particle& p){
    this->set(p.vx, p.vy, p.vz);
}
