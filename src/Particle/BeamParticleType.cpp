#include "particle_type.hpp"
#include "position.hpp"
#include "particle.hpp"
#include "normalizer.hpp"
#include <random>

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

    if (fabs(emission_vector[0]) == 1.0) {
        //! 放出方向と直交する平面上に幅をもたせる
        std::uniform_real_distribution<> dist_y(0.0, Normalizer::normalizeLength(emission_radius));
        std::uniform_real_distribution<> dist_z(0.0, Normalizer::normalizeLength(emission_radius));
        const double random_ywidth = dist_y(mt_y);
        const double random_zwidth = dist_z(mt_z);

        return Position{
            relative_emission_position.x + random_timewidth * vel.vx,
            relative_emission_position.y + random_timewidth * vel.vy + random_y_width,
            relative_emission_position.z + random_timewidth * vel.vz + random_z_width
        };
    } else if (fabs(emission_vector[1]) == 1.0) {
        std::uniform_real_distribution<> dist_x(0.0, Normalizer::normalizeLength(emission_radius));
        std::uniform_real_distribution<> dist_z(0.0, Normalizer::normalizeLength(emission_radius));
        const double random_xwidth = dist_x(mt_x);
        const double random_zwidth = dist_z(mt_z);

        return Position{
            relative_emission_position.x + random_timewidth * vel.vx + random_x_width,
            relative_emission_position.y + random_timewidth * vel.vy,
            relative_emission_position.z + random_timewidth * vel.vz + random_z_width
        };

    } else if (fabs(emission_vector[2]) == 1.0) {
        std::uniform_real_distribution<> dist_x(0.0, Normalizer::normalizeLength(emission_radius));
        std::uniform_real_distribution<> dist_y(0.0, Normalizer::normalizeLength(emission_radius));
        const double random_xwidth = dist_x(mt_x);
        const double random_ywidth = dist_y(mt_y);

        return Position{
            relative_emission_position.x + random_timewidth * vel.vx + random_x_width,
            relative_emission_position.y + random_timewidth * vel.vy + random_y_width,
            relative_emission_position.z + random_timewidth * vel.vz
        };
    } else {
        std::string error_message = (format("[ERROR] Now 'emission_vector' must be [+-1,0,0] or [0,+-1,0] or [0,0,+-1].")
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

    return Velocity{
        dist_vx(mt_vx) + velocity_coeff * emission_vector[0],
        dist_vy(mt_vy) + velocity_coeff * emission_vector[1],
        dist_vz(mt_vz) + velocity_coeff * emission_vector[2]
    };
}

double BeamParticleType::getEmissionAmount() const {
    return beam_current / (fabs(charge) * size);
}

void BeamParticleType::printInfo() const {
    ParticleType::printInfo();
    cout << "  emission_type: " << emission_type << endl;
    cout << "  beam_potential: " << Normalizer::unnormalizePotential(accel_potential) << "V" << endl;
    cout << "  beam_current: " << Normalizer::unnormalizeCurrent(beam_current) << "A" << endl;
    cout << "  emiss_amount: " << getEmissionAmount() << " particles / step" << endl;
    cout << "  beam_divergence: " << beam_divergence << "A" << endl;

    if (getAcceleration() > 1.0) {
        cout << "[WARNING] Acceleration for this particle type exceeds 1.0: " << getAcceleration() << "." << endl;
    }
}
