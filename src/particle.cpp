#include "tdpic.h"
Particle::Particle(){}
Particle::~Particle(){}

void Particle::setPosition(const Position& pos){
    x = pos.getX();
    y = pos.getY();
    z = pos.getZ();
}
void Particle::setPosition(double _x, double _y, double _z){
    x = _x;
    y = _y;
    z = _z;
}

void Particle::setVelocity(const Velocity& vel){
    vx = vel.getVX();
    vy = vel.getVY();
    vz = vel.getVZ();
}
void Particle::setVelocity(double _vx, double _vy, double _vz){
    vx = _vx;
    vy = _vy;
    vz = _vz;
}

double Particle::getX()  const { return x; }
double Particle::getY()  const { return y; }
double Particle::getZ()  const { return z; }
double Particle::getVX() const { return vx; }
double Particle::getVY() const { return vy; }
double Particle::getVZ() const { return vz; }

// ParticleInfo Class
ParticleInfo::ParticleInfo(){}

void ParticleInfo::setParticleSize(int size){
    sizeOfSuperParticle = size;
}

int ParticleInfo::getParticleSize(){
    return sizeOfSuperParticle;
}
