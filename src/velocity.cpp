#include "particle.hpp"

Velocity::Velocity(){}
Velocity::~Velocity(){}

void Velocity::set(double _vx, double _vy, double _vz){
    vx = _vx;
    vy = _vy;
    vz = _vz;
}

double Velocity::getVX() const { return vx; }
double Velocity::getVY() const { return vy; }
double Velocity::getVZ() const { return vz; }
