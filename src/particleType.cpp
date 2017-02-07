#include "global.hpp"
#include "environment.hpp"
#include "particle.hpp"
#include "utils.hpp"

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
