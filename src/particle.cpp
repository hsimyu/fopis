#include <tdpic.h>

Particle::Particle(void){
    isValid = 1;
}

// Copy Constructer
Particle::Particle(Particle const& p){
    x = p.x;
    y = p.y;
    z = p.z;
    vx = p.vx;
    vy = p.vy;
    vz = p.vz;
    typeId = p.typeId;
    isValid = p.isValid;
}

Particle& Particle::operator=(Particle const& rhs){
    x = rhs.x;
    y = rhs.y;
    z = rhs.z;
    vx = rhs.vx;
    vy = rhs.vy;
    vz = rhs.vz;
    typeId = rhs.typeId;
    isValid = rhs.isValid;
}

Particle::~Particle(){}

void Particle::setPosition(Position const& pos){
    x = pos.x;
    y = pos.y;
    z = pos.z;
}

void Particle::setPosition(const double _x, const double _y, const double _z){
    x = _x;
    y = _y;
    z = _z;
}

void Particle::setVelocity(Velocity const& vel){
    vx = vel.getVX();
    vy = vel.getVY();
    vz = vel.getVZ();
}

void Particle::setVelocity(const double _vx, const double _vy, const double _vz){
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
std::ostream& operator<<(std::ostream& ost, Particle const& p){
    ost << "[  Particle  ]" << std::endl;
    ost << "         id: " << p.typeId << std::endl;
    ost << "   position: " << format("%f, %f, %f") % p.x % p.y % p.z << std::endl;
    ost << "   velocity: " << format("%f, %f, %f") % p.vx % p.vy % p.vz << std::endl;
    ost << "    isValid: " << format("%s") % p.isValid << std::endl;
    return ost;
}
