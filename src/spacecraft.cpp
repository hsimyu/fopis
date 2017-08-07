#include "spacecraft.hpp"
#include "particle.hpp"
#include "normalizer.hpp"

//! static 変数の実体
unsigned int Spacecraft::num_of_spacecraft = 0;

void Spacecraft::construct(const size_t nx, const size_t ny, const size_t nz) {
    ++num_of_spacecraft;

    // Node ベース, glue cell は要らない?
    objectArray::extent_gen objectExtents;
    object_map.resize(objectExtents[nx][ny][nz]);

    // 初期化
    for(unsigned int i = 0; i < nx; ++i) {
        for (unsigned int j = 0; j < ny; ++j) {
            for (unsigned int k = 0; k < nz; ++k) {
                object_map[i][j][k] = false;
            }
        }
    }

    //! キャパシタンス行列の要素数
    num_cmat = 0;

    //! テスト定義
    for(unsigned int i = 3 * nx / 8; i < 4 * nx / 8; ++i) {
        for (unsigned int j = 3 * ny / 8; j < 4 * ny / 8; ++j) {
            for (unsigned int k = 3 * nz / 8; k < 4 * nz / 8; ++k) {
                //! 自分のノード内の座標である時だけ capacity_matrix_relation を追加する
                //! -> 後々 find() で自動的に処理できる
                object_map[i][j][k] = true;
                // std::map 内に直接オブジェクトを構築する
                capacity_matrix_relation.emplace(std::piecewise_construct, std::make_tuple(num_cmat), std::make_tuple(i, j, k));
                ++num_cmat;
            }
        }
    }

    //! キャパシタンス行列のサイズを物体サイズに変更
    capacity_matrix.resize(num_cmat, num_cmat);

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

Position Spacecraft::getCmatPos(const unsigned int cmat_itr) {
    return capacity_matrix_relation[cmat_itr];
}

auto Spacecraft::getTotalCharge(const tdArray& rho) const {
    double total_charge = 0.0;

    for(size_t cmat_itr = 0; cmat_itr < num_cmat; ++cmat_itr) {
        const auto& pos = capacity_matrix_relation.at(cmat_itr);
        total_charge += rho[pos.i][pos.j][pos.k];
    }

    return total_charge;
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

void Spacecraft::redistributeCharge(tdArray& rho, const tdArray& phi) const {
    cout << "charge before redist: " << getTotalCharge(rho) << endl;

    double charge_before_redist = 0.0;
    for(const auto& one_node : capacity_matrix_relation) {
        const auto j = one_node.first;
        const auto& pos = one_node.second;

        for(size_t i = 0; i < num_cmat; ++i) {
            charge_before_redist += capacity_matrix(i, j) * phi[pos.i][pos.j][pos.k];
        }
    }
    //! ここで sum とる

    const double new_potential = charge_before_redist / total_cmat_value;
    cout << "new potential = " << Normalizer::unnormalizePotential(new_potential) << " V. " << endl;

    for(unsigned int i = 0; i < num_cmat; ++i) {
        double delta_rho = 0.0;

        for(unsigned int j = 0; j < num_cmat; ++j) {
            const auto& pos = capacity_matrix_relation.at(j);
            delta_rho += capacity_matrix(i, j) * (new_potential - phi[pos.i][pos.j][pos.k]);
        }
        //! ここで sum とる

        const auto& target_pos = capacity_matrix_relation.at(i);
        rho[target_pos.i][target_pos.j][target_pos.k] += delta_rho;
    }
    cout << "charge after redist: " << getTotalCharge(rho) << endl;
}

std::ostream& operator<<(std::ostream& ost, const Spacecraft& spc) {
    ost << "    name:" << spc.name << endl;
    ost << "    cmat:" << spc.num_cmat << endl;

    return ost;
}
