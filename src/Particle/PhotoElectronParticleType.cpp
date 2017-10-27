#include "particle_type.hpp"
#include "position.hpp"
#include "particle.hpp"
#include "normalizer.hpp"
#include <random>

PhotoElectronParticleType::PhotoElectronParticleType() : EmissionParticleType{}, angle_gen{MPIw::Environment::rank} {}

Particle PhotoElectronParticleType::generateNewParticle(const Position& relative_emission_position, const std::array<double, 3>& emission_vector) {
    Particle p(id);
    Velocity vel = this->generateNewVelocity(emission_vector);
    p.setVelocity(vel);
    p.setPosition(this->generateNewPosition(relative_emission_position, emission_vector, vel));

    return p;
}

Position PhotoElectronParticleType::generateNewPosition(const Position& relative_emission_position, const std::array<double, 3>& emission_vector, const Velocity& vel) {
    std::uniform_real_distribution<> dist_t(0.0, 1.0);
    const double random_timewidth = dist_t(mt_x);
    incrementGeneratedCount(0);

    if (fabs(emission_vector[0]) == 1.0) {
        //! 放出方向と直交する平面1グリッド分に均等に粒子を分布させる
        std::uniform_real_distribution<> dist_y(0.0, 1.0);
        std::uniform_real_distribution<> dist_z(0.0, 1.0);
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
        std::uniform_real_distribution<> dist_x(0.0, 1.0);
        std::uniform_real_distribution<> dist_z(0.0, 1.0);
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
        std::uniform_real_distribution<> dist_x(0.0, 1.0);
        std::uniform_real_distribution<> dist_y(0.0, 1.0);
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
        std::string error_message = (format("[ERROR] PhotoElectronParticleType::generateNewPosition: Now 'emission_vector' must be [+-1,0,0] or [0,+-1,0] or [0,0,+-1].")).str();
        throw std::invalid_argument(error_message);
    }
}

Velocity PhotoElectronParticleType::generateNewVelocity(const std::array<double, 3>& emission_vector) {
    const double thermal_velocity = this->calcThermalVelocity();
    const double depression_angle = angle_gen.genDepressionAngle();
    const double azimuth_angle = angle_gen.genAzimuthAngle();

    if (fabs(emission_vector[0]) == 1.0) {
        Velocity v(
            thermal_velocity * std::cos(depression_angle) * emission_vector[0],
            thermal_velocity * std::sin(depression_angle) * std::cos(azimuth_angle),
            thermal_velocity * std::sin(depression_angle) * std::sin(azimuth_angle)
        );
        return v;
    } else if (fabs(emission_vector[1]) == 1.0) {
        Velocity v(
            thermal_velocity * std::sin(depression_angle) * std::sin(azimuth_angle),
            thermal_velocity * std::cos(depression_angle) * emission_vector[1],
            thermal_velocity * std::sin(depression_angle) * std::cos(azimuth_angle)
        );
        return v;
    } else if (fabs(emission_vector[2]) == 1.0) {
        Velocity v(
            thermal_velocity * std::sin(depression_angle) * std::cos(azimuth_angle),
            thermal_velocity * std::sin(depression_angle) * std::sin(azimuth_angle),
            thermal_velocity * std::cos(depression_angle) * emission_vector[2]
        );
        return v;
    } else {
        std::string error_message = (format("[ERROR] PhotoElectronParticleType::generateNewVelocity: Now 'emission_vector' must be [+-1,0,0] or [0,+-1,0] or [0,0,+-1].")).str();
        throw std::invalid_argument(error_message);
    }
}

//! @note: PhotoElectronの場合はCurrentDensityを返す
double PhotoElectronParticleType::getEmissionAmount() const {
    return current_density / (fabs(charge) * size);
}

void PhotoElectronParticleType::printInfo() const {
    ParticleType::printInfo();
    cout << "  current_density: " << Normalizer::unnormalizeCurrentDensity(this->getCurrentDensity()) << " A/m^2\n";
    cout << "  emiss_amount: " << this->getEmissionAmount() << " particles / (m^2 * step)" << endl;
}
