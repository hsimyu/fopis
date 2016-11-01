#include "tdpic.h"
#include <iostream>

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

// ParticleType Class
ParticleType::ParticleType(){}

// setter
void ParticleType::setId(int _id){ id = _id; }
void ParticleType::setCharge(double _charge){ charge = _charge; }
void ParticleType::setMass(double _mass){ mass = _mass; }
void ParticleType::setDensity(double _density){ density = _density; }
void ParticleType::setTemperature(double _temp){ temperature = _temp; }
void ParticleType::setName(std::string _name){ name = _name; }
void ParticleType::setType(std::string _type){ type = _type; }
void ParticleType::setSize(int _size){ size = _size; }
void ParticleType::setTotalNumber(int _num){ totalNumber = _num; }
void ParticleType::setPcell(int _pcell){ particle_per_cell = _pcell; }

// getter
int ParticleType::getId() const { return id; }
std::string ParticleType::getType() const { return type; }
std::string ParticleType::getName() const { return name; }
double ParticleType::getCharge() const { return charge; }
double ParticleType::getMass() const { return mass; }
double ParticleType::getDensity() const { return density; }
double ParticleType::getTemperature() const { return temperature; }
int ParticleType::getPcell() const { return particle_per_cell; }
int ParticleType::getSize() const { return size; }
int ParticleType::getTotalNumber() const { return totalNumber; }

int ParticleType::calcTotalNumber(const Environment* env, int nr){
    totalNumber = env->cell_x * env->cell_y * env->cell_z * nr;
    return totalNumber;
}

// util
std::ostream& operator<<(std::ostream& ost, const ParticleType& ptype){
    ost << "[Particle  : " << ptype.getName() << "]" << std::endl;
    ost << "         id: " << ptype.getId() << std::endl;
    ost << "       type: " << ptype.getType() << std::endl;
    ost << "       mass: " << ptype.getMass() << std::endl;
    ost << "     charge: " << ptype.getCharge() << std::endl;
    ost << "    density: " << ptype.getDensity() << std::endl;
    ost << "temperature: " << ptype.getTemperature() << std::endl;
    ost << "       size: " << ptype.getSize() << std::endl;
    ost << "   per_cell: " << ptype.getPcell() << std::endl;
    ost << "totalNumber: " << ptype.getTotalNumber() << std::endl;
    return ost;
}
