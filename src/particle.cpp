#include <tdpic.h>

Particle::Particle(){
    // std::cout << "Particle Constructer is called" << std::endl;
}
Particle::~Particle(){
    // std::cout << "Particle Destructor is called" << std::endl;
}

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

void Particle::setTypeId(const int _id) { typeId = _id; }
int Particle::getTypeId(void) const { return typeId; }

//! 個別計算用
double Particle::getEnergy(void) const {
    return 0.5 * (vx*vx + vy*vy + vz*vz) * Environment::ptype[typeId].getMass();
}

//! まとめて計算する時用
double Particle::getSquaredMagnitudeOfVelocity(void) const {
    return (vx*vx + vy*vy + vz*vz);
}
