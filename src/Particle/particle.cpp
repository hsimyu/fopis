#include "global.hpp"
#include "position.hpp"
#include "particle.hpp"
#include "environment.hpp"

std::shared_ptr<ParticleType> Particle::getParticleTypePtr() const {
    return Environment::getParticleType(typeId);
}

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

Position Particle::getPosition(void) const {
    return Position{x, y, z};
}

Position Particle::getOldPosition(void) const {
    return Position{x - vx, y - vy, z - vz};
}

Position Particle::getNextPosition(void) const {
    return Position{x + vx, y + vy, z + vz};
}

Particle Particle::getOldPositionParticle() const {
    return Particle{typeId, x - vx, y - vy, z - vz, vx, vy, vz};
}

Particle Particle::getNextPositionParticle() const {
    return Particle{typeId, x + vx, y + vy, z + vz, vx, vy, vz};
}

Velocity Particle::getVelocity(void) const {
    return Velocity{vx, vy, vz};
}

double Particle::getXMoveRatio() const {
    auto pos = this->getPosition();

    if (vx == 0.0) return 0.0;

    if (vx > 0.0) {
        return pos.dx1 / std::fabs(vx);
    } else {
        return pos.dx2 / std::fabs(vx);
    }
}

double Particle::getYMoveRatio() const {
    auto pos = this->getPosition();

    if (vy == 0.0) return 0.0;

    if (vy > 0.0) {
        return pos.dy1 / std::fabs(vy);
    } else {
        return pos.dy2 / std::fabs(vy);
    }
}

double Particle::getZMoveRatio() const {
    auto pos = this->getPosition();

    if (vz == 0.0) return 0.0;

    if (vz > 0.0) {
        return pos.dz1 / std::fabs(vz);
    } else {
        return pos.dz2 / std::fabs(vz);
    }
}

std::vector<Position> Particle::computeCrossPoints() const {
    std::vector<Position> cross_points;

    const auto pos = getPosition();
    const auto old_pos = getOldPosition();

    //! 次の表面までの距離(速度で正規化)
    const double mvx = 1.0 - getXMoveRatio();
    const double mvy = 1.0 - getYMoveRatio();
    const double mvz = 1.0 - getZMoveRatio();

    unsigned int count = 0;

    Position xcross_point;
    bool move_along_x = (pos.i != old_pos.i);
    if (move_along_x) {
        if (vx > 0.0) {
            xcross_point.setXYZ(std::ceil(old_pos.x), old_pos.y + vy * mvx, old_pos.z + vz * mvx);
        } else {
            xcross_point.setXYZ(std::floor(old_pos.x), old_pos.y + vy * mvx, old_pos.z + vz * mvx);
        }
        ++count;
    }

    Position ycross_point;
    bool move_along_y = (pos.j != old_pos.j);
    if (move_along_y) {
        if (vy > 0.0) {
            ycross_point.setXYZ(old_pos.x + vx * mvy, std::ceil(old_pos.y), old_pos.z + vz * mvy);
        } else {
            ycross_point.setXYZ(old_pos.x + vx * mvy, std::floor(old_pos.y), old_pos.z + vz * mvy);
        }
        ++count;
    }

    Position zcross_point;
    bool move_along_z = (pos.k != old_pos.k);
    if (move_along_z) {
        if (vz > 0.0) {
            zcross_point.setXYZ(old_pos.x + vx * mvz, old_pos.y + vy * mvz, std::ceil(old_pos.z));
        } else {
            zcross_point.setXYZ(old_pos.x + vx * mvz, old_pos.y + vy * mvz, std::floor(old_pos.z));
        }
        ++count;
    }

    if (count == 3) {
        if (mvx <= mvy) {
            if (mvx <= mvz) {
                cross_points.push_back( std::move(xcross_point) );

                if (mvy <= mvz) {
                    cross_points.push_back( std::move(ycross_point) );
                    cross_points.push_back( std::move(zcross_point) );
                } else {
                    cross_points.push_back( std::move(zcross_point) );
                    cross_points.push_back( std::move(ycross_point) );
                }
            } else {
                cross_points.push_back( std::move(zcross_point) );
                cross_points.push_back( std::move(xcross_point) );
                cross_points.push_back( std::move(ycross_point) );
            }
        } else {
            if (mvy <= mvz) {
                cross_points.push_back( std::move(ycross_point) );

                if (mvx <= mvz) {
                    cross_points.push_back( std::move(xcross_point) );
                    cross_points.push_back( std::move(zcross_point) );
                } else {
                    cross_points.push_back( std::move(zcross_point) );
                    cross_points.push_back( std::move(xcross_point) );
                }
            } else {
                cross_points.push_back( std::move(zcross_point) );
                cross_points.push_back( std::move(ycross_point) );
                cross_points.push_back( std::move(xcross_point) );
            }
        }
    } else if (count == 2) {
        if (move_along_x) {
            if (move_along_y) {
                // xy
                if (mvx <= mvy) {
                    cross_points.push_back( std::move(xcross_point) );
                    cross_points.push_back( std::move(ycross_point) );
                } else {
                    cross_points.push_back( std::move(ycross_point) );
                    cross_points.push_back( std::move(xcross_point) );
                }
            } else {
                // xz
                if (mvx <= mvx) {
                    cross_points.push_back( std::move(xcross_point) );
                    cross_points.push_back( std::move(zcross_point) );
                } else {
                    cross_points.push_back( std::move(zcross_point) );
                    cross_points.push_back( std::move(xcross_point) );
                }
            }
        } else {
            // yz
            if (mvy <= mvz) {
                cross_points.push_back( std::move(ycross_point) );
                cross_points.push_back( std::move(zcross_point) );
            } else {
                cross_points.push_back( std::move(zcross_point) );
                cross_points.push_back( std::move(ycross_point) );
            }
        }
    } else if (count == 1) {
        if (move_along_x) {
            cross_points.push_back( std::move(xcross_point) );
        } else if (move_along_y) {
            cross_points.push_back( std::move(ycross_point) );
        } else {
            cross_points.push_back( std::move(zcross_point) );
        }
    }

    return cross_points;
}

Position Particle::getNextXCrossPoint() const {
    Position pos = getPosition();

    if (vx > 0.0) {
        const double mvx = pos.dx2 / vx;
        pos.setXYZ(std::ceil(pos.x), pos.y + vy * mvx, pos.z + vz * mvx);
        return pos;
    } else if (vx < 0.0) {
        const double mvx = pos.dx1 / (-vx);
        pos.setXYZ(std::floor(pos.x), pos.y + vy * mvx, pos.z + vz * mvx);
        return pos;
    } else {
        std::cerr << "[WARNING] Vx of this particle is equal to 0.0 at Particle::getNextXCrossPoint()." << endl;
        return pos;
    }
}

Position Particle::getNextYCrossPoint() const {
    Position pos = getPosition();

    if (vy > 0.0) {
        const double mvy = pos.dy2 / vy;
        pos.setXYZ(pos.x + vx * mvy, std::ceil(pos.y), pos.z + vz * mvy);
        return pos;
    } else if (vy < 0.0) {
        const double mvy = pos.dy1 / (-vy);
        pos.setXYZ(pos.x + vx * mvy, std::floor(pos.y), pos.z + vz * mvy);
        return pos;
    } else {
        std::cerr << "[WARNING] Vy of this particle is equal to 0.0 at Particle::getNextYCrossPoint()." << endl;
        return pos;
    }
}

Position Particle::getNextZCrossPoint() const {
    Position pos = getPosition();

    if (vz > 0.0) {
        const double mvz = pos.dz2 / vz;
        pos.setXYZ(pos.x + vx * mvz, pos.y + vy * mvz, std::ceil(pos.z));
        return pos;
    } else if (vz < 0.0) {
        const double mvz = pos.dz1 / (-vz);
        pos.setXYZ(pos.x + vx * mvz, pos.y + vy * mvz, std::floor(pos.z));
        return pos;
    } else {
        std::cerr << "[WARNING] Vz of this particle is equal to 0.0 at Particle::getNeztZCrossPoint()." << endl;
        return pos;
    }
}

//! 電流配分
void Particle::distributeCurrentAtOldPosition(const double q_per_dt, tdArray& jx, tdArray& jy, tdArray& jz) const {
    const Position old_pos = getPosition();
    const Position new_pos = getNextPosition();
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
    return 0.5 * (vx*vx + vy*vy + vz*vz) * (this->getParticleTypePtr()->getMass());
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
