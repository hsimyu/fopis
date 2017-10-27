#include "particle_type.hpp"
#include "position.hpp"
#include "particle.hpp"
#include "normalizer.hpp"
#include <random>

/*
Particle BeamParticleType::generateNewParticle(const Position& relative_emission_position, const std::array<double, 3>& emission_vector) {
    Particle p(id);
    Velocity vel = this->generateNewVelocity(emission_vector);
    p.setVelocity(vel);
    p.setPosition(this->generateNewPosition(relative_emission_position, emission_vector, vel));

    return p;
}

Position BeamParticleType::generateNewPosition(const Position& relative_emission_position, const std::array<double, 3>& emission_vector, const Velocity& vel) {
    std::uniform_real_distribution<> dist_t(0.0, 1.0);
    const double random_timewidth = dist_t(mt_x);
    incrementGeneratedCount(0);

    if (fabs(emission_vector[0]) == 1.0) {
        //! 放出方向と直交する平面上に幅をもたせる
        std::uniform_real_distribution<> dist_y(0.0, Normalizer::normalizeLength(emission_radius));
        std::uniform_real_distribution<> dist_z(0.0, Normalizer::normalizeLength(emission_radius));
        const double random_ywidth = dist_y(mt_y);
        const double random_zwidth = dist_z(mt_z);
        incrementGeneratedCount(1);
        incrementGeneratedCount(2);

        //! 放出位置を -1 velocity しておくことで、粒子更新時に時間が同期するようにする
        return Position{
            relative_emission_position.x + (random_timewidth - 1.0) * vel.vx,
            relative_emission_position.y + (random_timewidth - 1.0) * vel.vy + random_ywidth,
            relative_emission_position.z + (random_timewidth - 1.0) * vel.vz + random_zwidth
        };
    } else if (fabs(emission_vector[1]) == 1.0) {
        std::uniform_real_distribution<> dist_x(0.0, Normalizer::normalizeLength(emission_radius));
        std::uniform_real_distribution<> dist_z(0.0, Normalizer::normalizeLength(emission_radius));
        const double random_xwidth = dist_x(mt_x);
        const double random_zwidth = dist_z(mt_z);
        incrementGeneratedCount(0);
        incrementGeneratedCount(2);

        return Position{
            relative_emission_position.x + (random_timewidth - 1.0) * vel.vx + random_xwidth,
            relative_emission_position.y + (random_timewidth - 1.0) * vel.vy,
            relative_emission_position.z + (random_timewidth - 1.0) * vel.vz + random_zwidth
        };

    } else if (fabs(emission_vector[2]) == 1.0) {
        std::uniform_real_distribution<> dist_x(0.0, Normalizer::normalizeLength(emission_radius));
        std::uniform_real_distribution<> dist_y(0.0, Normalizer::normalizeLength(emission_radius));
        const double random_xwidth = dist_x(mt_x);
        const double random_ywidth = dist_y(mt_y);
        incrementGeneratedCount(0);
        incrementGeneratedCount(1);

        return Position{
            relative_emission_position.x + (random_timewidth - 1.0) * vel.vx + random_xwidth,
            relative_emission_position.y + (random_timewidth - 1.0) * vel.vy + random_ywidth,
            relative_emission_position.z + (random_timewidth - 1.0) * vel.vz
        };
    } else {
        std::string error_message = (format("[ERROR] At BeamParticleType::generateNewPosition: Now 'emission_vector' must be [+-1,0,0] or [0,+-1,0] or [0,0,+-1].")).str();
        throw std::invalid_argument(error_message);
    }
}

Velocity BeamParticleType::generateNewVelocity(const std::array<double, 3>& emission_vector) {
    const double deviation = this->calcDeviation();
    std::normal_distribution<> dist_vx(0.0, deviation);
    std::normal_distribution<> dist_vy(0.0, deviation);
    std::normal_distribution<> dist_vz(0.0, deviation);

    double emission_normal = sqrt(pow(emission_vector[0], 2) + pow(emission_vector[1], 2) + pow(emission_vector[2], 2));
    double velocity_coeff = getAcceleration() / emission_normal;

    incrementVelocityGeneratedCount();

    return Velocity{
        dist_vx(mt_vx) + velocity_coeff * emission_vector[0],
        dist_vy(mt_vy) + velocity_coeff * emission_vector[1],
        dist_vz(mt_vz) + velocity_coeff * emission_vector[2]
    };
}

double BeamParticleType::getAcceleration() const {
    return sqrt(2.0 * accel_potential / mass);
}
*/

double PhotoElectronParticleType::getEmissionAmount() const {
    return current_density / (fabs(charge) * size);
}

void PhotoElectronParticleType::printInfo() const {
    ParticleType::printInfo();
    cout << "  current_density: " << Normalizer::unnormalizeCurrentDensity(this->getCurrentDensity()) << " A/m^2\n";
    cout << "  emiss_amount: " << this->getEmissionAmount() << " particles / step" << endl;
}
