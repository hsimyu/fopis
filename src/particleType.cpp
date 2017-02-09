#include "global.hpp"
#include "environment.hpp"
#include "position.hpp"
#include "particle.hpp"
#include "grid.hpp"
#include "utils.hpp"
#include <random>

ParticleType::ParticleType(void) :
    mt_x(10684930 + MPIw::Environment::rank),
    mt_y(99881 + MPIw::Environment::rank),
    mt_z(861200045 + MPIw::Environment::rank),
    mt_vx(930 + MPIw::Environment::rank),
    mt_vy(98076621 + MPIw::Environment::rank),
    mt_vz(7662566 + MPIw::Environment::rank) {}

int ParticleType::calcSize(void){
    size = static_cast<int>(pow(Environment::dx, 3) * density / static_cast<double>(particle_per_cell));
    return size;
}

// initializer for ambient plasma
int ParticleType::calcTotalNumber(void){
    totalNumber = Environment::cell_x * Environment::cell_y * Environment::cell_z * particle_per_cell;
    return totalNumber;
}

double ParticleType::calcThermalVelocity() const {
    return Utils::Normalizer::normalizeVelocity(sqrt(2.0f * temperature / (mass * me)));
}

Position ParticleType::generateNewPosition(
        const double min_x, const double max_x,
        const double min_y, const double max_y,
        const double min_z, const double max_z) {

    std::uniform_real_distribution<> dist_x(min_x, max_x);
    std::uniform_real_distribution<> dist_y(min_y, max_y);
    std::uniform_real_distribution<> dist_z(min_z, max_z);

    Position p(dist_x(mt_x), dist_y(mt_y), dist_z(mt_z));
    return p;
}

Velocity ParticleType::generateNewVelocity(void) {
    const double deviation = this->calcDeviation();
    std::normal_distribution<> dist_vx(0.0, deviation);
    std::normal_distribution<> dist_vy(0.0, deviation);
    std::normal_distribution<> dist_vz(0.0, deviation);

    Velocity v(dist_vx(mt_vx), dist_vy(mt_vy), dist_vz(mt_vz));
    return v;
}

std::vector<double> ParticleType::calcFlux(Grid const& g) const {
    //
    // const double ddv = 1e-3;
    // const double nvx = 12.0f/ddv;
    //
    // double fluxx1 = 0.0f;
    // double fluxx2 = 0.0f;
    //
    // //! -6.0から6.0まで積分
    // //! x * exp(-x^2)?
    // //! @note: これはv = 0を境目にして、どちらの方向へのフラックスが多いのかを計算する
    // //! 完全に等方的 (drift速度がゼロ) の場合は0.5ずつ
    // for(int i = 0; i < nvx - 1; ++i) {
    //     double arg1 = i*ddv - 6.0f;
    //     double argx1 = pow(arg1, 2);
    //
    //     if(arg1 > 0.0f) {
    //         fluxx1 += arg1 * exp(-argx1) * ddv;
    //     } else {
    //         fluxx2 += arg1 * exp(-argx1) * ddv;
    //     }
    // }
    //
    // fluxx1 *= this->calcThermalVelocity()/sqrt(M_PI);
    // fluxx2 *= -1.0f * this->calcThermalVelocity()/sqrt(M_PI);

    std::vector<double> flux(6);
    //! 等方的なfluxを仮定
    const double flux0 = 0.5f * this->calcThermalVelocity()/sqrt(M_PI) * Utils::Normalizer::normalizeDensity(this->getDensity());
    const double areax = pow(Utils::Normalizer::normalizeLength(Environment::dx), 2) * g.getNZ() * g.getNY();
    const double areay = pow(Utils::Normalizer::normalizeLength(Environment::dx), 2) * g.getNZ() * g.getNX();
    const double areaz = pow(Utils::Normalizer::normalizeLength(Environment::dx), 2) * g.getNX() * g.getNY();

    flux[0] = areax*flux0 / this->getSize();
    flux[1] = areax*flux0 / this->getSize();
    flux[2] = areay*flux0 / this->getSize();
    flux[3] = areay*flux0 / this->getSize();
    flux[4] = areaz*flux0 / this->getSize();
    flux[5] = areaz*flux0 / this->getSize();

    return flux;
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
