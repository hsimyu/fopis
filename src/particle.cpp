#include <tdpic.h>

Particle::Particle(){}

// Copy Constructer
Particle::Particle(const Particle& p){
    x = p.x;
    y = p.y;
    z = p.z;
    vx = p.vx;
    vy = p.vy;
    vz = p.vz;
    typeId = p.typeId;
}

Particle::~Particle(){
    std::cout << "Particle Destructor is called" << std::endl;
}

void Particle::setPosition(const Position& pos){
    x = pos.x;
    y = pos.y;
    z = pos.z;
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

//! 個別計算用
double Particle::getEnergy(void) const {
    return 0.5 * (vx*vx + vy*vy + vz*vz) * Environment::ptype[typeId].getMass();
}

//! まとめて計算する時用
double Particle::getSquaredMagnitudeOfVelocity(void) const {
    return (vx*vx + vy*vy + vz*vz);
}

//! 位置の更新
void Particle::updatePosition(void) {
    x += vx;
    y += vy;
    z += vz;
}

// util
std::ostream& operator<<(std::ostream& ost, const Particle& p){
    ost << "[  Particle  ]" << std::endl;
    ost << "         id: " << p.typeId << std::endl;
    ost << "   position: " << format("%f, %f, %f") % p.x % p.y % p.z << std::endl;
    ost << "   velocity: " << format("%f, %f, %f") % p.vx % p.vy % p.vz << std::endl;
    return ost;
}
