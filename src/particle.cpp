#include <tdpic.h>

Particle::Particle(){
    // std::cout << "Particle Constructer is called" << std::endl;
}
Particle::~Particle(){
    // std::cout << "Particle Destructor is called" << std::endl;
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

double Particle::getX()  const { return x; }
double Particle::getY()  const { return y; }
double Particle::getZ()  const { return z; }
double Particle::getVX() const { return vx; }
double Particle::getVY() const { return vy; }
double Particle::getVZ() const { return vz; }

void Particle::setVX(const double _vx) { vx = _vx; }
void Particle::setVY(const double _vy) { vy = _vy; }
void Particle::setVZ(const double _vz) { vz = _vz; }

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

//! 位置の更新
void Particle::updatePosition(void) {
    x += vx;
    y += vy;
    z += vz;
}
