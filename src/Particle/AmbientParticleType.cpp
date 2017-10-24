#include "particle_type.hpp"
#include "position.hpp"
#include "particle.hpp"
#include "grid.hpp"
#include "normalizer.hpp"

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
    const double areax = pow(Normalizer::normalizeLength(Environment::dx), 2) * (g.getNZ() - 1) * (g.getNY() - 1);
    const double areay = pow(Normalizer::normalizeLength(Environment::dx), 2) * (g.getNZ() - 1) * (g.getNX() - 1);
    const double areaz = pow(Normalizer::normalizeLength(Environment::dx), 2) * (g.getNX() - 1) * (g.getNY() - 1);

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

Velocity AmbientParticleType::generateNewInjectionVelocity(AXIS axis, AXIS_SIDE low_or_up) {
    const double deviation = this->calcDeviation();
    std::normal_distribution<> dist_vx(0.0, deviation);
    std::normal_distribution<> dist_vy(0.0, deviation);
    std::normal_distribution<> dist_vz(0.0, deviation);

    double sign = (low_or_up == AXIS_SIDE::low) ? 1.0 : -1.0;

    Velocity v(dist_vx(mt_vx), dist_vy(mt_vy), dist_vz(mt_vz));
    switch(axis) {
        case AXIS::x:
            while ((v.vx * sign) < 0.0) {
                v.vx = dist_vx(mt_vx);
                v.vy = dist_vy(mt_vy);
                v.vz = dist_vz(mt_vz);
                incrementVelocityGeneratedCount();
            }
            return v;
        case AXIS::y:
            while ((v.vy * sign) < 0.0) {
                v.vx = dist_vx(mt_vx);
                v.vy = dist_vy(mt_vy);
                v.vz = dist_vz(mt_vz);
                incrementVelocityGeneratedCount();
            }
            return v;
        case AXIS::z:
            while ((v.vz * sign) < 0.0) {
                v.vx = dist_vx(mt_vx);
                v.vy = dist_vy(mt_vy);
                v.vz = dist_vz(mt_vz);
                incrementVelocityGeneratedCount();
            }
            return v;
    }
    throw std::logic_error("[ERROR] Something wrong in AmbientParticleType::generateNewInjectionVelocity().");
}