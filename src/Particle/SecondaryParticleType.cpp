#include "particle_type.hpp"
#include "position.hpp"
#include "particle.hpp"
#include "normalizer.hpp"
#include <Spacecraft/Material.hpp>
#include <Spacecraft/Incident.hpp>
#include <random>

SecondaryParticleType::SecondaryParticleType() : EmissionParticleType{},
    angle_gen{MPIw::Environment::rank},
    mt_ptype_gen(884751 + 2 * MPIw::Environment::rank),
    mt_true_see_gen(85 + 7 * MPIw::Environment::rank) {}

std::vector<Particle> SecondaryParticleType::generateNewParticles(const IncidentInfo_t& incident, const MaterialInfo_t& material) {
    std::vector<Particle> parray;

    // Velocity vel = this->generateNewVelocity(emission_vector);
    // p.setVelocity(vel);
    // p.setPosition(this->generateNewPosition(relative_emission_position, emission_vector, vel));

    const auto energy = incident.getIncidentEnergyInElectronVolt();
    const auto angle = incident.getIncidentAngle();
    const double r = this->getElasticBackscatterCoeff(material, energy, angle);
    const double eta = this->getInelasticBackscatterCoeff(material, energy, angle);
    double delta = this->getTrueSecondaryCoeff(material, energy, angle);

    if (Environment::isRootNode) {
        cout << format("incident energy = %s\n") % energy;
        cout << format("incident angle = %s\n") % ((180.0 / M_PI) * angle);
        cout << format("r = %s\n") % r;
        cout << format("eta = %s\n") % eta;
        cout << format("delta = %s") % delta << endl;
    }

    const double emission_type_random = dist_uniform(mt_ptype_gen);
    cout << format("random = %s") % emission_type_random << endl;
    if (emission_type_random < r) {
        cout << format("This is elastic backscattering.") << endl;
        Particle p(id);
        parray.push_back( std::move(p) );
    } else if (emission_type_random < (r + eta)) {
        cout << format("This is inelastic backscattering.") << endl;
        Particle p(id);
        parray.push_back( std::move(p) );
    } else {
        cout << format("This is true secondary emisssion.") << endl;
        //! deltaを "一時粒子に対する真の二次電子の数"から"吸収された（後方散乱されない）粒子の数に対する二次粒子の数に変換"
        //! By K.Hoshi et al. 2017
        delta /= (1.0 - r - eta);

        unsigned int emission_count = std::floor(delta);
        const double delta_resid = delta - std::floor(delta);

        const double emission_count_random = dist_uniform(mt_true_see_gen);
        if (emission_count_random < delta_resid) {
            emission_count++;
        }

        cout << format("Delta = %s.\n") % delta;
        cout << format("Emission Count = %s.\n") % emission_count << endl;

        // for (unsigned int itr = 0; itr < emission_count; ++itr) {
        // }
    }

    return parray;
}

/*
Position SecondaryParticleType::generateNewPosition(const Position& relative_emission_position, const std::array<double, 3>& emission_vector, const Velocity& vel) {
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
        std::string error_message = (format("[ERROR] SecondaryParticleType::generateNewPosition: Now 'emission_vector' must be [+-1,0,0] or [0,+-1,0] or [0,0,+-1].")).str();
        throw std::invalid_argument(error_message);
    }
}

Velocity SecondaryParticleType::generateNewVelocity(const std::array<double, 3>& emission_vector) {
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
        std::string error_message = (format("[ERROR] SecondaryParticleType::generateNewVelocity: Now 'emission_vector' must be [+-1,0,0] or [0,+-1,0] or [0,0,+-1].")).str();
        throw std::invalid_argument(error_message);
    }
}*/

void SecondaryParticleType::printInfo() const {
    ParticleType::printInfo();
}

double SecondaryParticleType::getTrueSecondaryCoeff(const MaterialInfo_t& material, const double incident_energy, const double incident_angle) const {
    return this->computeWhippleTrueSecondaryCoeff(material, incident_energy, incident_angle);
};

double SecondaryParticleType::getElasticBackscatterCoeff(const MaterialInfo_t& material, const double incident_energy, const double incident_angle) const {
    return this->computeWhippleElasticCoeff(material, incident_energy, incident_angle);
};

double SecondaryParticleType::getInelasticBackscatterCoeff(const MaterialInfo_t& material, const double incident_energy, const double incident_angle) const {
    return this->computeWhippleInelasticCoeff(material, incident_energy, incident_angle);
};

//! 各係数の実装
double SecondaryParticleType::computeWhippleTrueSecondaryCoeff(const MaterialInfo_t& material, const double incident_energy, const double incident_angle) const {
    const double delta = (1.114 * material.delta_max / std::cos(incident_angle)) * std::pow(material.epsi_max / incident_energy, 0.35) * (1.0 - std::exp(-2.28 * std::cos(incident_angle) * std::pow(incident_energy / material.epsi_max, 1.35)));
    return delta;
};

double SecondaryParticleType::computeWhippleInelasticCoeff(const MaterialInfo_t& material, const double incident_energy, const double incident_angle) const {
    //! 10keV以上向けの係数
    const double eta1 = 1.0 - std::pow(2.0 / incident_energy, 0.037 * material.atomic_number);

    //! 10keV以下向けの係数
    const double eta2 = 0.1 * std::exp(-incident_energy / 5000.0);

    //! from Whipple, Eq. (3.17)
    double eta0 = 0.0;
    if (incident_energy > 50.0 && incident_energy < 1000.0) {
        eta0 = std::log(incident_energy / 50.0) / std::log(20.0);
    } else if (incident_energy > 1000.0) {
        eta0 = 1.0;
    }

    //! empirical factor for angular dependent obtained by Laframboise and Kamitsuma, 1983
    return eta0 * (eta1 + eta2) * (7.37 * std::pow(material.atomic_number, -0.56875)) * (1.0 - std::cos(incident_angle));
};

double SecondaryParticleType::computeWhippleElasticCoeff(const MaterialInfo_t& material, const double incident_energy, const double incident_angle) const {
    constexpr double rydberg_factor = 1.0 / 13.54;
    return rydberg_factor * std::pow(material.fermi_energy, 4) / (16.0 * (std::pow(incident_energy + material.fermi_energy, 3) + std::pow(material.fermi_energy, 3)));
}