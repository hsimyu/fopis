#include "global.hpp"
#include "position.hpp"
#include "particle.hpp"
#include "environment.hpp"
#include "utils.hpp"

void Particle::setPosition(Position const& pos){
    x = pos.x;
    y = pos.y;
    z = pos.z;
}

void Particle::setVelocity(Velocity const& v){
    vx = v.vx;
    vy = v.vy;
    vz = v.vz;
}

Position&& Particle::getPosition(void) const {
    return Position{x, y, z};
}

Position&& Particle::getOldPosition(void) const {
    return Position{x - vx, y - vy, z - vz};
}

Position&& Particle::getNewPosition(void) const {
    return Position{x + vx, y + vy, z + vz};
}

//! 電流配分
void Particle::distributeCurrentAtOldPosition(const double q_per_dt, tdArray& jx, tdArray& jy, tdArray& jz) const {
    const Position old_pos = getPosition();
    const Position new_pos = getNewPosition();
    const Position ref_pos = new_pos.getReferencePosition(old_pos);

    // charge flux の計算
    const double flux1[3] = {
        q_per_dt * (ref_pos.x - old_pos.x),
        q_per_dt * (ref_pos.y - old_pos.y),
        q_per_dt * (ref_pos.z - old_pos.z),
    };

    //! 1次の shape factor の計算
    //! 中点における shape factor を取る
    //! グリッドをまたいだ時 ref_position の dx1 は 1.0 となるので、
    const double average_w1[3] = {
        0.5 * (ref_pos.x + old_pos.x) - (static_cast<double>(old_pos.i) - 1.0),
        0.5 * (ref_pos.y + old_pos.y) - (static_cast<double>(old_pos.j) - 1.0),
        0.5 * (ref_pos.z + old_pos.z) - (static_cast<double>(old_pos.k) - 1.0),
    };

    const int pos_i = old_pos.i;
    const int pos_j = old_pos.j;
    const int pos_k = old_pos.k;

    // Current Conserving Weighting, l = 1 (old cell)
    jx[pos_i][pos_j    ][pos_k    ] += flux1[0] * (1.0 - average_w1[1]) * (1.0 - average_w1[2]);
    jx[pos_i][pos_j + 1][pos_k    ] += flux1[0] * average_w1[1] * (1.0 - average_w1[2]);
    jx[pos_i][pos_j    ][pos_k + 1] += flux1[0] * (1.0 - average_w1[1]) * average_w1[2];
    jx[pos_i][pos_j + 1][pos_k + 1] += flux1[0] * average_w1[1] * average_w1[2];

    jy[pos_i    ][pos_j][pos_k    ] += flux1[1] * (1.0 - average_w1[0]) * (1.0 - average_w1[2]);
    jy[pos_i + 1][pos_j][pos_k    ] += flux1[1] * average_w1[0] * (1.0 - average_w1[2]);
    jy[pos_i    ][pos_j][pos_k + 1] += flux1[1] * (1.0 - average_w1[0]) * average_w1[2];
    jy[pos_i + 1][pos_j][pos_k + 1] += flux1[1] * average_w1[0] * average_w1[2];

    jz[pos_i    ][pos_j    ][pos_k] += flux1[2] * (1.0 - average_w1[0]) * (1.0 - average_w1[1]);
    jz[pos_i + 1][pos_j    ][pos_k] += flux1[2] * average_w1[0] * (1.0 - average_w1[1]);
    jz[pos_i    ][pos_j + 1][pos_k] += flux1[2] * (1.0 - average_w1[0]) * average_w1[1];
    jz[pos_i + 1][pos_j + 1][pos_k] += flux1[2] * average_w1[0] * average_w1[1];
}

//! 新しい位置における電流配分
void Particle::distributeCurrentAtNewPosition(const double q_per_dt, tdArray& jx, tdArray& jy, tdArray& jz) const {
    const Position old_pos = getOldPosition();
    const Position new_pos = getPosition();
    const Position ref_pos = new_pos.getReferencePosition(old_pos);

    // charge flux の計算
    const double flux2[3] = {
        q_per_dt * (new_pos.x - ref_pos.x),
        q_per_dt * (new_pos.y - ref_pos.y),
        q_per_dt * (new_pos.z - ref_pos.z),
    };

    const double average_w2[3] = {
        0.5 * (ref_pos.x + new_pos.x) - (static_cast<double>(new_pos.i) - 1.0),
        0.5 * (ref_pos.y + new_pos.y) - (static_cast<double>(new_pos.j) - 1.0),
        0.5 * (ref_pos.z + new_pos.z) - (static_cast<double>(new_pos.k) - 1.0),
    };

    const int pos_i = new_pos.i;
    const int pos_j = new_pos.j;
    const int pos_k = new_pos.k;

    // Current Conserving Weighting, l = 2(new cell)
    jx[pos_i][pos_j    ][pos_k    ] += flux2[0] * (1.0 - average_w2[1]) * (1.0 - average_w2[2]);
    jx[pos_i][pos_j + 1][pos_k    ] += flux2[0] * average_w2[1] * (1.0 - average_w2[2]);
    jx[pos_i][pos_j    ][pos_k + 1] += flux2[0] * (1.0 - average_w2[1]) * average_w2[2];
    jx[pos_i][pos_j + 1][pos_k + 1] += flux2[0] * average_w2[1] * average_w2[2];

    jy[pos_i    ][pos_j][pos_k    ] += flux2[1] * (1.0 - average_w2[0]) * (1.0 - average_w2[2]);
    jy[pos_i + 1][pos_j][pos_k    ] += flux2[1] * average_w2[0] * (1.0 - average_w2[2]);
    jy[pos_i    ][pos_j][pos_k + 1] += flux2[1] * (1.0 - average_w2[0]) * average_w2[2];
    jy[pos_i + 1][pos_j][pos_k + 1] += flux2[1] * average_w2[0] * average_w2[2];

    jz[pos_i    ][pos_j    ][pos_k] += flux2[2] * (1.0 - average_w2[0]) * (1.0 - average_w2[1]);
    jz[pos_i + 1][pos_j    ][pos_k] += flux2[2] * average_w2[0] * (1.0 - average_w2[1]);
    jz[pos_i    ][pos_j + 1][pos_k] += flux2[2] * (1.0 - average_w2[0]) * average_w2[1];
    jz[pos_i + 1][pos_j + 1][pos_k] += flux2[2] * average_w2[0] * average_w2[1];
}

double Particle::getEnergy(void) const {
    return 0.5 * (vx*vx + vy*vy + vz*vz) * Environment::ptype[typeId].getMass();
}

void Particle::generateNewPosition(const double min_x, const double max_x, const double min_y, const double max_y, const double min_z, const double max_z)
{
    this->setPosition( Environment::ptype[typeId].generateNewPosition(min_x, max_x, min_y, max_y, min_z, max_z) );
}

void Particle::generateNewVelocity(void)
{
    this->setVelocity( Environment::ptype[typeId].generateNewVelocity() );
}

// util
std::ostream& operator<<(std::ostream& ost, Particle const& p){
    ost << "[  Particle  ]" << endl;
    ost << "         id: " << p.typeId << endl;
    ost << "   position: " << format("%f, %f, %f") % p.x % p.y % p.z << endl;
    ost << "   velocity: " << format("%f, %f, %f") % p.vx % p.vy % p.vz << endl;
    ost << "    isValid: " << format("%s") % p.isValid << endl;
    return ost;
}
