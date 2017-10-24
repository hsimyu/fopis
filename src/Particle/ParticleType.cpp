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
    cout << "  current density (OML): " << (e * density * sqrt(temperature / (2.0 * M_PI * mass * me))) << '\n';
    cout << "  size: " << getSize() << endl;
    cout << "  per_cell: " << getPcell() << endl;
    cout << "  total_number in node: " << getTotalNumber() << endl;
    cout << "  total number in all node: " << (getTotalNumber() * Environment::proc_x * Environment::proc_y * Environment::proc_z) << endl;
    cout << "  memory / process: " << calcMemory() << endl;
}
