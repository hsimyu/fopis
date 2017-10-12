#include "global.hpp"
#include "environment.hpp"
#include "particle_type.hpp"
#include "position.hpp"
#include "particle.hpp"
#include "grid.hpp"
#include "utils.hpp"
#include "normalizer.hpp"
#include <random>

ParticleType::ParticleType(void) :
    mt_x(10684930 + MPIw::Environment::rank),
    mt_y(99881 + MPIw::Environment::rank),
    mt_z(861200045 + MPIw::Environment::rank),
    mt_vx(930 + MPIw::Environment::rank),
    mt_vy(98076621 + MPIw::Environment::rank),
    mt_vz(7662566 + MPIw::Environment::rank) {}

void ParticleType::setId(int _id){
    id = _id;
    mt_x.seed(10684930 + MPIw::Environment::rank + 106*_id);
    mt_y.seed(99881 + MPIw::Environment::rank + 171*_id);
    mt_z.seed(861200045 + MPIw::Environment::rank + 95*_id);
    mt_vx.seed(930 + MPIw::Environment::rank + 2*_id);
    mt_vy.seed(98076621 + MPIw::Environment::rank - 50 * _id);
    mt_vz.seed(7662566 + MPIw::Environment::rank + 4 * _id);
}

double ParticleType::calcDebyeLength(void) const {
    return sqrt( temperature * eps0 / (pow(e, 2) * density));
}

int ParticleType::updateSize(void) {
    size = static_cast<int>(pow(Environment::dx, 3) * density / static_cast<double>(particle_per_cell));
    return size;
}

// initializer for ambient plasma
int ParticleType::updateTotalNumber(void){
    //! 担当するセルの数は上側境界にいるかどうかで変わる
    int cellX = (Environment::onHighXedge) ? Environment::cell_x - 1 : Environment::cell_x;
    int cellY = (Environment::onHighYedge) ? Environment::cell_y - 1 : Environment::cell_y;
    int cellZ = (Environment::onHighZedge) ? Environment::cell_z - 1 : Environment::cell_z;
    totalNumber = cellX * cellY * cellZ * particle_per_cell;

    if (totalNumber > Environment::max_particle_num) {
        if (Environment::isRootNode) {
            cout << format("[WARNING] total particle number %d for '%s' exceeds embedded max particle number %d.") % totalNumber % name % Environment::max_particle_num << endl;
        }

        totalNumber = Environment::max_particle_num;
    }
    return totalNumber;
}

double ParticleType::calcThermalVelocity() const {
    return Normalizer::normalizeVelocity(sqrt(2.0 * temperature / (mass * me)));
}

double ParticleType::calcDeviation() const {
    return Normalizer::normalizeVelocity(sqrt(temperature / (mass * me)));
}

double ParticleType::calcPlasmaFrequency(void) const {
    return sqrt(density * pow(e, 2) / (me * eps0));
}

std::string ParticleType::calcMemory() const {
    static constexpr double memory_per_particle = 8.0*6 + 4.0*2;
    static constexpr double memory_buff_coeff = 2.0;

    double pmem = this->getTotalNumber() * memory_per_particle * memory_buff_coeff;

    return Utils::prettyMemoryString(pmem);
}

void ParticleType::proceedGeneratedCounts(const ParticleType::GeneratedCount_t& target_counts) {
    assert(target_counts.size() == generated_counts.size());

    for(size_t i = 0; i < generated_counts.size(); ++i) {
        const size_t discard_number = target_counts[i] - generated_counts[i];

        if (discard_number > 0) {
            switch(i) {
                case 0:
                    mt_x.discard(discard_number);
                    break;
                case 1:
                    mt_y.discard(discard_number);
                    break;
                case 2:
                    mt_z.discard(discard_number);
                    break;
                case 3:
                    mt_vx.discard(discard_number);
                    break;
                case 4:
                    mt_vy.discard(discard_number);
                    break;
                case 5:
                    mt_vz.discard(discard_number);
                    break;
                default:
                    break;
            }
        }

        generated_counts[i] = target_counts[i];
    }
}

// util
void ParticleType::printInfo() const {
    cout << "[ParticleType : " << getName() << "]" << endl;
    cout << "  id: " << getId() << endl;
    cout << "  type: " << getType() << endl;
    cout << "  mass: " << getMass() << endl;
    cout << "  charge: " << getCharge() << "e" << endl;
    cout << "  density: " << getDensity() << "/m^3" << endl;
    cout << "  temperature: " << getTemperature() << "eV" << endl;
    cout << "  size: " << getSize() << endl;
    cout << "  per_cell: " << getPcell() << endl;
    cout << "  total_number in node: " << getTotalNumber() << endl;
    cout << "  total number in all node: " << (getTotalNumber() * Environment::proc_x * Environment::proc_y * Environment::proc_z) << endl;
    cout << "  memory / process: " << calcMemory() << endl;
}

//! 背景プラズマ用
std::vector<double> AmbientParticleType::calcFlux(Grid const& g) const {
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
    const double flux0 = 0.5 * this->calcThermalVelocity()/sqrt(M_PI) * Normalizer::normalizeDensity(this->getDensity());
    const double areax = pow(Normalizer::normalizeLength(Environment::dx), 2) * g.getNZ() * g.getNY();
    const double areay = pow(Normalizer::normalizeLength(Environment::dx), 2) * g.getNZ() * g.getNX();
    const double areaz = pow(Normalizer::normalizeLength(Environment::dx), 2) * g.getNX() * g.getNY();

    flux[0] = areax*flux0 / this->getSize();
    flux[1] = areax*flux0 / this->getSize();
    flux[2] = areay*flux0 / this->getSize();
    flux[3] = areay*flux0 / this->getSize();
    flux[4] = areaz*flux0 / this->getSize();
    flux[5] = areaz*flux0 / this->getSize();

    return flux;
}

//! 粒子生成Factory関数の実体
Particle AmbientParticleType::generateNewParticle(const double min_x, const double max_x, const double min_y, const double max_y, const double min_z, const double max_z) {
    Particle p(id);
    p.setPosition(this->generateNewPosition(min_x, max_x, min_y, max_y, min_z, max_z));
    p.setVelocity(this->generateNewVelocity());

    return p;
}

Particle AmbientParticleType::generateNewParticle(const double min_x, const double max_x, const double min_y, const double max_y, const double min_z, const double max_z, const Velocity& vel) {
    Particle p(id);
    p.setPosition(this->generateNewPosition(min_x, max_x, min_y, max_y, min_z, max_z));
    p.setVelocity(vel);

    return p;
}

Position AmbientParticleType::generateNewPosition(
        const double min_x, const double max_x,
        const double min_y, const double max_y,
        const double min_z, const double max_z) {

    std::uniform_real_distribution<> dist_x(min_x, max_x);
    std::uniform_real_distribution<> dist_y(min_y, max_y);
    std::uniform_real_distribution<> dist_z(min_z, max_z);

    Position p(dist_x(mt_x), dist_y(mt_y), dist_z(mt_z));
    incrementPositionGeneratedCount();
    return p;
}

Velocity AmbientParticleType::generateNewVelocity(void) {
    const double deviation = this->calcDeviation();
    std::normal_distribution<> dist_vx(0.0, deviation);
    std::normal_distribution<> dist_vy(0.0, deviation);
    std::normal_distribution<> dist_vz(0.0, deviation);

    Velocity v(dist_vx(mt_vx), dist_vy(mt_vy), dist_vz(mt_vz));
    incrementVelocityGeneratedCount();
    return v;
}