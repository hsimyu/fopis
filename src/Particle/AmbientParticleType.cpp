#include "particle_type.hpp"
#include "position.hpp"
#include "particle.hpp"
#include "grid.hpp"
#include "normalizer.hpp"

//! 背景プラズマ用
std::vector<double> AmbientParticleType::calcFlux() const {
    const double ddv = 1e-3;
    const double nvx = 12.0 / ddv;

    std::vector<double> flux(6);
    const auto thermal_velocity = this->calcThermalVelocity();
    const auto normed_density = Normalizer::normalizeDensity(this->getDensity());
    const auto size = this->getSize();

    for(int axis = 0; axis < 3; ++axis) {
        double drift_velocity_to_thermal_velocity;
        switch(axis) {
            case 0:
                drift_velocity_to_thermal_velocity = drift_velocity.vx / thermal_velocity;
                break;
            case 1:
                drift_velocity_to_thermal_velocity = drift_velocity.vy / thermal_velocity;
                break;
            case 2:
                drift_velocity_to_thermal_velocity = drift_velocity.vz / thermal_velocity;
                break;
        }

        //! @note: これはv = 0を境目にして、どちらの方向へのフラックスが多いのかを計算する
        //! x * exp(-x^2) = \int x f(x) = 1方向への速度の平均速度 = |\hat{v_x}|
        //! 暫定的に-6.0から6.0まで積分した割合
        //! 完全に等方的 (drift速度がゼロ) の場合は0.5ずつになる
        for(int i = 0; i < nvx - 1; ++i) {
            double arg1 = i * ddv + (drift_velocity_to_thermal_velocity - 6.0);
            double argx1 = pow(arg1 - drift_velocity_to_thermal_velocity, 2);
        
            if(arg1 > 0.0) {
                flux[2 * axis] += arg1 * exp(-argx1) * ddv;
            } else {
                flux[2 * axis + 1] += arg1 * exp(-argx1) * ddv;
            }
        }

        flux[2 * axis] *= (thermal_velocity / sqrt(M_PI)) * normed_density / size;
        flux[2 * axis + 1] *= -1.0 * (thermal_velocity / sqrt(M_PI)) * normed_density / size;
    }

    //! injection axis が false に設定されている場合は
    //! その方向からの注入なしとする
    for(int axis = 0; axis < flux.size(); ++axis) {
        if (!injection_axis[axis]) flux[axis] = 0.0;
    }

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
    std::normal_distribution<> dist_vx(drift_velocity.vx, deviation);
    std::normal_distribution<> dist_vy(drift_velocity.vy, deviation);
    std::normal_distribution<> dist_vz(drift_velocity.vz, deviation);

    Velocity v(dist_vx(mt_vx), dist_vy(mt_vy), dist_vz(mt_vz));
    incrementVelocityGeneratedCount();
    return v;
}

Velocity AmbientParticleType::generateNewInjectionVelocity(AXIS axis, AXIS_SIDE low_or_up) {
    const double deviation = this->calcDeviation();
    const double thermal_velocity = this->calcThermalVelocity();

    double sign = (low_or_up == AXIS_SIDE::low) ? 1.0 : -1.0;

    switch(axis) {
        case AXIS::x:
            {
                std::normal_distribution<> dist_vx(drift_velocity.vx, thermal_velocity);
                std::normal_distribution<> dist_vy(drift_velocity.vy, deviation);
                std::normal_distribution<> dist_vz(drift_velocity.vz, deviation);

                Velocity v(dist_vx(mt_vx), dist_vy(mt_vy), dist_vz(mt_vz));
                while ((v.vx * sign) < 0.0) {
                    v.vx = dist_vx(mt_vx);
                    v.vy = dist_vy(mt_vy);
                    v.vz = dist_vz(mt_vz);
                    incrementVelocityGeneratedCount();
                }
                return v;
            }
        case AXIS::y:
            {
                std::normal_distribution<> dist_vx(drift_velocity.vx, deviation);
                std::normal_distribution<> dist_vy(drift_velocity.vy, thermal_velocity);
                std::normal_distribution<> dist_vz(drift_velocity.vz, deviation);

                Velocity v(dist_vx(mt_vx), dist_vy(mt_vy), dist_vz(mt_vz));

                while ((v.vy * sign) < 0.0) {
                    v.vx = dist_vx(mt_vx);
                    v.vy = dist_vy(mt_vy);
                    v.vz = dist_vz(mt_vz);
                    incrementVelocityGeneratedCount();
                }
                return v;
            }
        case AXIS::z:
            {
                std::normal_distribution<> dist_vx(drift_velocity.vx, deviation);
                std::normal_distribution<> dist_vy(drift_velocity.vy, deviation);
                std::normal_distribution<> dist_vz(drift_velocity.vz, thermal_velocity);

                Velocity v(dist_vx(mt_vx), dist_vy(mt_vy), dist_vz(mt_vz));

                while ((v.vz * sign) < 0.0) {
                    v.vx = dist_vx(mt_vx);
                    v.vy = dist_vy(mt_vy);
                    v.vz = dist_vz(mt_vz);
                    incrementVelocityGeneratedCount();
                }
                return v;
            }
    }
    throw std::logic_error("[ERROR] Something wrong in AmbientParticleType::generateNewInjectionVelocity().");
}

void AmbientParticleType::printInfo() const {
    ParticleType::printInfo();

    const auto velocity_unnorm = Normalizer::unnormalizeVelocity(1.0);
    cout << format("  drift_velocity: %s km/s, %s km/s, %s km/s\n") %
        (drift_velocity.vx * velocity_unnorm * 1e-3) %
        (drift_velocity.vy * velocity_unnorm * 1e-3) %
        (drift_velocity.vz * velocity_unnorm * 1e-3);
    cout << format("  injection_axis: -x: %s, +x: %s, -y: %s, +y: %s, -z: %s, +z: %s\n") %
        (injection_axis[0] ? "ON" : "OFF") %
        (injection_axis[1] ? "ON" : "OFF") %
        (injection_axis[2] ? "ON" : "OFF") %
        (injection_axis[3] ? "ON" : "OFF") %
        (injection_axis[4] ? "ON" : "OFF") %
        (injection_axis[5] ? "ON" : "OFF") << endl;
}
