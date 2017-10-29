#include "particle_type.hpp"
#include "position.hpp"
#include "particle.hpp"
#include "normalizer.hpp"
#include <Spacecraft/Material.hpp>
#include <Spacecraft/Incident.hpp>
#include <random>

SecondaryParticleType::SecondaryParticleType() : EmissionParticleType{},
    truesee_angle_gen{MPIw::Environment::rank},
    inelastic_angle_gen{851035553 + 5 * MPIw::Environment::rank},
    mt_ptype_gen(884751 + 2 * MPIw::Environment::rank),
    mt_inelastic_gen(985739508 + 13 * MPIw::Environment::rank),
    mt_true_see_gen(85 + 7 * MPIw::Environment::rank) {}

std::vector<Particle> SecondaryParticleType::generateNewParticles(const IncidentInfo_t& incident, const MaterialInfo_t& material) {
    std::vector<Particle> parray;

    const auto energy = incident.getIncidentEnergyInElectronVolt();
    const auto angle = incident.getIncidentAngle();
    const double r = (incident.getCharge() < 0.0) ? this->getElasticBackscatterCoeff(material, energy, angle) : 0.0;
    const double eta = (incident.getCharge() < 0.0) ? this->getInelasticBackscatterCoeff(material, energy, angle) : 0.0;
    double delta = this->getTrueSecondaryCoeff(material, energy, angle);

    const double emission_type_random = dist_uniform(mt_ptype_gen);

    //! 確率的二次電子放出の実装
    //! K. Hoshi et al. IEEE Trans. on Plasma Science, 2017
    if (emission_type_random < r) {
        parray.push_back( this->generateElasticBackscatterParticle(incident, material) );
    } else if (emission_type_random < (r + eta)) {
        parray.push_back( this->generateInelasticBackscatterParticle(incident, material) );
    } else {
        //! deltaを "一時粒子に対する真の二次電子の数"から"吸収された（後方散乱されない）粒子の数に対する二次粒子の数に変換"
        delta /= (1.0 - r - eta);

        unsigned int emission_count = std::floor(delta);
        const double delta_resid = delta - std::floor(delta);

        const double emission_count_random = dist_uniform(mt_true_see_gen);
        if (emission_count_random < delta_resid) {
            emission_count++;
        }

        for (unsigned int itr = 0; itr < emission_count; ++itr) {
            parray.push_back( this->generateTrueSecondaryParticle(incident, material) );
        }
    }

    return parray;
}

Particle SecondaryParticleType::generateElasticBackscatterParticle(const IncidentInfo_t& incident, const MaterialInfo_t& material) {
    Particle p(id);
    const double remain_time = incident.getRemainingTime();
    const Position pos = incident.getPosition();

    if (incident.isXsurfaceIncident()) {
        p.setVelocity(-incident.getVx(), incident.getVy(), incident.getVz());
        p.setPosition(pos.x + remain_time * p.vx, pos.y + remain_time * p.vy, pos.z + remain_time * p.vz);
        return p;
    } else if (incident.isYsurfaceIncident()) {
        p.setVelocity(incident.getVx(), -incident.getVy(), incident.getVz());
        p.setPosition(pos.x + remain_time * p.vx, pos.y + remain_time * p.vy, pos.z + remain_time * p.vz);
        return p;
    } else {
        p.setVelocity(incident.getVx(), incident.getVy(), -incident.getVz());
        p.setPosition(pos.x + remain_time * p.vx, pos.y + remain_time * p.vy, pos.z + remain_time * p.vz);
        return p;
    }
}

Particle SecondaryParticleType::generateInelasticBackscatterParticle(const IncidentInfo_t& incident, const MaterialInfo_t& material) {
    Particle p(id);
    Velocity v;
    const double emission_energy = dist_uniform(mt_inelastic_gen) * incident.getIncidentEnergy();
    const double emission_velocity = Normalizer::convertEnergyToVelocity(emission_energy, mass);
    const double depression_angle = inelastic_angle_gen.genDepressionAngle();
    const double azimuth_angle = inelastic_angle_gen.genAzimuthAngle();

    if (incident.isXsurfaceIncident()) {
        const double inverse_sign = (incident.getVx() > 0.0) ? -1.0 : 1.0;
        v.set(
            emission_velocity * std::cos(depression_angle) * inverse_sign,
            emission_velocity * std::sin(depression_angle) * std::cos(azimuth_angle),
            emission_velocity * std::sin(depression_angle) * std::sin(azimuth_angle)
        );
    } else if (incident.isYsurfaceIncident()) {
        const double inverse_sign = (incident.getVy() > 0.0) ? -1.0 : 1.0;
        v.set(
            emission_velocity * std::sin(depression_angle) * std::sin(azimuth_angle),
            emission_velocity * std::cos(depression_angle) * inverse_sign,
            emission_velocity * std::sin(depression_angle) * std::cos(azimuth_angle)
        );
    } else {
        const double inverse_sign = (incident.getVz() > 0.0) ? -1.0 : 1.0;
        v.set(
            emission_velocity * std::sin(depression_angle) * std::cos(azimuth_angle),
            emission_velocity * std::sin(depression_angle) * std::sin(azimuth_angle),
            emission_velocity * std::cos(depression_angle) * inverse_sign
        );
    }
    p.setVelocity(v);

    const double remain_time = incident.getRemainingTime();
    const Position pos = incident.getPosition();
    p.setPosition(pos.x + remain_time * p.vx, pos.y + remain_time * p.vy, pos.z + remain_time * p.vz);

    return p;
}

Particle SecondaryParticleType::generateTrueSecondaryParticle(const IncidentInfo_t& incident, const MaterialInfo_t& material) {
    Particle p(id);
    Velocity v;

    //! 二次電子の放出エネルギーは2eVの熱速度とする
    const double emission_energy = Normalizer::normalizeEnergy(temperature);
    const double emission_velocity = Normalizer::convertEnergyToVelocity(emission_energy, mass);
    const double depression_angle = truesee_angle_gen.genDepressionAngle();
    const double azimuth_angle = truesee_angle_gen.genAzimuthAngle();

    if (incident.isXsurfaceIncident()) {
        const double inverse_sign = (incident.getVx() > 0.0) ? -1.0 : 1.0;
        v.set(
            emission_velocity * std::cos(depression_angle) * inverse_sign,
            emission_velocity * std::sin(depression_angle) * std::cos(azimuth_angle),
            emission_velocity * std::sin(depression_angle) * std::sin(azimuth_angle)
        );
    } else if (incident.isYsurfaceIncident()) {
        const double inverse_sign = (incident.getVy() > 0.0) ? -1.0 : 1.0;
        v.set(
            emission_velocity * std::sin(depression_angle) * std::sin(azimuth_angle),
            emission_velocity * std::cos(depression_angle) * inverse_sign,
            emission_velocity * std::sin(depression_angle) * std::cos(azimuth_angle)
        );
    } else {
        const double inverse_sign = (incident.getVz() > 0.0) ? -1.0 : 1.0;
        v.set(
            emission_velocity * std::sin(depression_angle) * std::cos(azimuth_angle),
            emission_velocity * std::sin(depression_angle) * std::sin(azimuth_angle),
            emission_velocity * std::cos(depression_angle) * inverse_sign
        );
    }
    p.setVelocity(v);

    const double remain_time = incident.getRemainingTime();
    const Position pos = incident.getPosition();
    p.setPosition(pos.x + remain_time * p.vx, pos.y + remain_time * p.vy, pos.z + remain_time * p.vz);

    return p;
}

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