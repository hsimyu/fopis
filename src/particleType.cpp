#include "global.hpp"
#include "environment.hpp"
#include "particle.hpp"
#include "utils.hpp"

// ParticleType Class
ParticleType::ParticleType(){}

// setter
void ParticleType::setId(int _id){ id = _id; }
void ParticleType::setCharge(double _charge){ charge = _charge; }
void ParticleType::setMass(double _mass){ mass = _mass; }
void ParticleType::setDensity(double _density){ density = _density; }
void ParticleType::setTemperature(double _temp){
    //! eV形式で入力されると仮定
    //! 内部的にはkB Teの値で持つ => eをかけて保存
    temperature = e * _temp;
}
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
double ParticleType::getTemperature() const { return temperature / e; }
double ParticleType::getTrueTemperature() const { return temperature; }
int ParticleType::getPcell() const { return particle_per_cell; }
int ParticleType::getSize() const { return size; }
int ParticleType::getTotalNumber() const { return totalNumber; }

int ParticleType::calcSize(void){
    size = static_cast<int>(pow(Environment::dx, 3) * density / static_cast<double>(particle_per_cell));
    return size;
}

// initializer for ambient plasma
int ParticleType::calcTotalNumber(void){
    totalNumber = Environment::cell_x * Environment::cell_y * Environment::cell_z * particle_per_cell;
    return totalNumber;
}

double ParticleType::calcDeviation() const {
    return Utils::Normalizer::normalizeVelocity(sqrt(temperature / (mass * me)));
}

std::string ParticleType::calcMemory() const {
    static constexpr double memory_per_particle = 8.0*6 + 4.0*2;

    double pmem = (this->getTotalNumber() > Environment::max_particle_num) ?
        this->getTotalNumber() * memory_per_particle : Environment::max_particle_num * memory_per_particle;

    return Utils::prettyMemoryString(pmem);
}

// util
std::ostream& operator<<(std::ostream& ost, const ParticleType& ptype){
    ost << "[ParticleType  : " << ptype.getName() << "]" << endl;
    ost << "             id: " << ptype.getId() << endl;
    ost << "           type: " << ptype.getType() << endl;
    ost << "           mass: " << ptype.getMass() << endl;
    ost << "         charge: " << ptype.getCharge() << "e" << endl;
    ost << "        density: " << ptype.getDensity() << "/m^3" << endl;
    ost << "    temperature: " << ptype.getTemperature() << "eV" << endl;
    ost << "           size: " << ptype.getSize() << endl;
    ost << "       per_cell: " << ptype.getPcell() << endl;
    ost << "    totalNumber: " << ptype.getTotalNumber() << endl;
    ost << "         memory: " << ptype.calcMemory() << endl << endl;
    return ost;
}
