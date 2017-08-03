#include "spacecraft.hpp"
#include "position.hpp"
#include "particle.hpp"

//! static 変数の実体
unsigned int Spacecraft::num_of_spacecraft = 0;

void Spacecraft::construct(const size_t nx, const size_t ny, const size_t nz) {
    ++num_of_spacecraft;

    // Node ベース, glue cell は要らない?
    objectArray::extent_gen objectExtents;
    object_map.resize(objectExtents[nx][ny][nz]);

    // 初期化
    for(size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t k = 0; k < nz; ++k) {
                object_map[i][j][k] = false;
            }
        }
    }

    //! キャパシタンス行列の要素数
    num_cmat = 0;

    //! テスト定義
    for(size_t i = nx/4; i < 3 * nx/4; ++i) {
        for (size_t j = ny / 4; j < 3 * ny / 4; ++j) {
            for (size_t k = nz / 4; k < 3 * nz / 4; ++k) {
                object_map[i][j][k] = true;
                capacity_matrix_relation[num_cmat][0] = i;
                capacity_matrix_relation[num_cmat][1] = j;
                capacity_matrix_relation[num_cmat][2] = k;
                ++num_cmat;
            }
        }
    }

    //! キャパシタンス行列のサイズを物体サイズに変更
    Cmatrix::extent_gen CmatExtents;
    capacity_matrix.resize(CmatExtents[num_cmat][num_cmat]);

    //! Node ベース, glue cell ありの電荷密度マップ
    tdArray::extent_gen tdExtents;
    charge_map.resize(tdExtents[nx + 2][ny + 2][nz + 2]);
    //! 電荷配列初期化
    for(size_t i = 0; i < nx + 2; ++i) {
        for (size_t j = 0; j < ny + 2; ++j) {
            for (size_t k = 0; k < nz + 2; ++k) {
                charge_map[i][j][k] = 0.0;
            }
        }
    }
}

bool Spacecraft::isIncluded(const Particle& p) const {
    const auto pos = p.getPosition();

    return  object_map[pos.raw_i    ][pos.raw_j    ][pos.raw_k    ] &&
            object_map[pos.raw_i + 1][pos.raw_j    ][pos.raw_k    ] &&
            object_map[pos.raw_i    ][pos.raw_j + 1][pos.raw_k    ] &&
            object_map[pos.raw_i + 1][pos.raw_j + 1][pos.raw_k    ] &&
            object_map[pos.raw_i    ][pos.raw_j    ][pos.raw_k + 1] &&
            object_map[pos.raw_i + 1][pos.raw_j    ][pos.raw_k + 1] &&
            object_map[pos.raw_i    ][pos.raw_j + 1][pos.raw_k + 1] &&
            object_map[pos.raw_i + 1][pos.raw_j + 1][pos.raw_k + 1];
}

void Spacecraft::removeInnerParticle(Particle& p) const {
    if (isIncluded(p)) { p.makeInvalid(); }
}

void Spacecraft::distributeInnerParticleCharge(Particle& p) {
    if (isIncluded(p)) { 
        const auto q = p.getCharge();
        const auto pos = p.getPosition();

        charge_map[pos.i    ][pos.j    ][pos.k    ] += q * pos.dx2 * pos.dy2 * pos.dz2;
        charge_map[pos.i + 1][pos.j    ][pos.k    ] += q * pos.dx1 * pos.dy2 * pos.dz2;
        charge_map[pos.i    ][pos.j + 1][pos.k    ] += q * pos.dx2 * pos.dy1 * pos.dz2;
        charge_map[pos.i + 1][pos.j + 1][pos.k    ] += q * pos.dx1 * pos.dy1 * pos.dz2;
        charge_map[pos.i    ][pos.j    ][pos.k + 1] += q * pos.dx2 * pos.dy2 * pos.dz1;
        charge_map[pos.i + 1][pos.j    ][pos.k + 1] += q * pos.dx1 * pos.dy2 * pos.dz1;
        charge_map[pos.i    ][pos.j + 1][pos.k + 1] += q * pos.dx2 * pos.dy1 * pos.dz1;
        charge_map[pos.i + 1][pos.j + 1][pos.k + 1] += q * pos.dx1 * pos.dy1 * pos.dz1;

        p.makeInvalid();
    }
}

void Spacecraft::applyCharge(tdArray& rho) const {
    //! 電荷分布を場に印加する
    for(size_t i = 0; i < charge_map.shape()[0]; ++i) {
        for (size_t j = 0; j < charge_map.shape()[1]; ++j) {
            for (size_t k = 0; k < charge_map.shape()[2]; ++k) {
                rho[i][j][k] += charge_map[i][j][k];
            }
        }
    }
}

std::ostream& operator<<(std::ostream& ost, const Spacecraft& spc) {
    ost << "name:" << spc.getName() << endl;

    return ost;
}
