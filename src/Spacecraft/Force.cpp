#include "particle.hpp"
#include "spacecraft.hpp"
#include "mpiw.hpp"
#include "normalizer.hpp"

//! 推力計算関連
void Spacecraft::updateMaxwellTensorForce(const tdArray& exref, const tdArray& eyref, const tdArray& ezref) {
    Force coulomb{0.0, 0.0, 0.0};

    for(const auto& node : capacity_matrix_relation) {
        const auto cmat_number = node.first;
        const auto& pos = node.second;

        const double t_ex = exref[pos.i][pos.j][pos.k];
        const double t_ey = eyref[pos.i][pos.j][pos.k];
        const double t_ez = ezref[pos.i][pos.j][pos.k];
        const double t_e2 = 0.5 * (t_ex * t_ex + t_ey * t_ey + t_ez * t_ez);

        const auto srf = this->getSurfaceVector(cmat_number);

        coulomb.fx += srf[0] * (t_ex * t_ex - t_e2) + srf[1] * (t_ey * t_ex) + srf[2] * (t_ez * t_ex);
        coulomb.fy += srf[0] * (t_ex * t_ey) + srf[1] * (t_ey * t_ey - t_e2) + srf[2] * (t_ez * t_ey);
        coulomb.fz += srf[0] * (t_ex * t_ez) + srf[1] * (t_ey * t_ez) + srf[2] * (t_ez * t_ez - t_e2);
    }

    force += coulomb * Normalizer::eps0;

    force.fx = MPIw::Environment::Comms[name].sum(force.fx);
    force.fy = MPIw::Environment::Comms[name].sum(force.fy);
    force.fz = MPIw::Environment::Comms[name].sum(force.fz);
}

void Spacecraft::resetForce() {
    force.fx = 0.0;
    force.fy = 0.0;
    force.fz = 0.0;
}

void Spacecraft::accumulateIncidentForce(const Particle& p) {
    auto mass = p.getMass();
    force.fx += p.vx * mass;
    force.fy += p.vy * mass;
    force.fz += p.vz * mass;
}

void Spacecraft::subtractEmissionForce(const Particle& p) {
    auto mass = p.getMass();
    force.fx -= p.vx * mass;
    force.fy -= p.vy * mass;
    force.fz -= p.vz * mass;
}