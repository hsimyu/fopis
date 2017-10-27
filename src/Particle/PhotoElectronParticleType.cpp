#include "particle_type.hpp"
#include "position.hpp"
#include "particle.hpp"
#include "normalizer.hpp"
#include <random>

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
    const double deviation = this->calcDeviation();
    const double thermal_velocity = this->calcThermalVelocity();

    if (fabs(emission_vector[0]) == 1.0) {
        std::normal_distribution<> dist_vx(0.0, thermal_velocity);
        std::normal_distribution<> dist_vy(0.0, deviation);
        std::normal_distribution<> dist_vz(0.0, deviation);

        Velocity v(dist_vx(mt_vx), dist_vy(mt_vy), dist_vz(mt_vz));
        incrementVelocityGeneratedCount();

        double sign = std::copysign(1.0, emission_vector[0]);
        while ((v.vx * sign) < 0.0) {
            v.vx = dist_vx(mt_vx);
            v.vy = dist_vy(mt_vy);
            v.vz = dist_vz(mt_vz);
            incrementVelocityGeneratedCount();
        }
        return v;
    } else if (fabs(emission_vector[1]) == 1.0) {
        std::normal_distribution<> dist_vx(0.0, deviation);
        std::normal_distribution<> dist_vy(0.0, thermal_velocity);
        std::normal_distribution<> dist_vz(0.0, deviation);

        Velocity v(dist_vx(mt_vx), dist_vy(mt_vy), dist_vz(mt_vz));
        incrementVelocityGeneratedCount();

        double sign = std::copysign(1.0, emission_vector[1]);
        while ((v.vy * sign) < 0.0) {
            v.vx = dist_vx(mt_vx);
            v.vy = dist_vy(mt_vy);
            v.vz = dist_vz(mt_vz);
            incrementVelocityGeneratedCount();
        }
        return v;
    } else if (fabs(emission_vector[2]) == 1.0) {
        std::normal_distribution<> dist_vx(0.0, deviation);
        std::normal_distribution<> dist_vy(0.0, deviation);
        std::normal_distribution<> dist_vz(0.0, thermal_velocity);

        Velocity v(dist_vx(mt_vx), dist_vy(mt_vy), dist_vz(mt_vz));
        incrementVelocityGeneratedCount();

        double sign = std::copysign(1.0, emission_vector[2]);
        while ((v.vz * sign) < 0.0) {
            v.vx = dist_vx(mt_vx);
            v.vy = dist_vy(mt_vy);
            v.vz = dist_vz(mt_vz);
            incrementVelocityGeneratedCount();
        }
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
