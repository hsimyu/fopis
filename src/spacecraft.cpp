#include "spacecraft.hpp"
#include "particle.hpp"
#include "normalizer.hpp"
#include "utils.hpp"

//! static 変数の実体
unsigned int Spacecraft::num_of_spacecraft = 0;

void Spacecraft::construct(const size_t nx, const size_t ny, const size_t nz, const ObjectNodes& nodes, const ObjectNodes& glue_nodes) {
    //! このオブジェクトがプロセス内で有効かどうかを保存しておく
    is_defined_in_this_process = (nodes.size() > 0);
    ++num_of_spacecraft;
    potential = 0.0;
    potential_fix = 0.0;

    if (is_defined_in_this_process) {
        // Node ベース, Glueセルも必要
        ObjectDefinedMap::extent_gen objectExtents;
        object_map.resize(objectExtents[nx + 2][ny + 2][nz + 2]);

        // 初期化
        for(int i = 0; i < nx + 2; ++i) {
            for (int j = 0; j < ny + 2; ++j) {
                for (int k = 0; k < nz + 2; ++k) {
                    object_map[i][j][k] = false;
                }
            }
        }

        //! 判定はGrid側で既に終わっているので必要なし
        for(const auto& node_pair : nodes) {
            const auto cmat_itr = node_pair.first;
            const auto& node_pos = node_pair.second;
            const auto i = node_pos[0];
            const auto j = node_pos[1];
            const auto k = node_pos[2];

            object_map[i][j][k] = true;
            capacity_matrix_relation.emplace(std::piecewise_construct, std::make_tuple(cmat_itr), std::make_tuple(i, j, k));
        }

        //! Glueノードはobject_map側を更新するだけでよい
        for(const auto& node_pair : glue_nodes) {
            const auto& node_pos = node_pair.second;
            const auto i = node_pos[0];
            const auto j = node_pos[1];
            const auto k = node_pos[2];

            object_map[i][j][k] = true;
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
}

Position Spacecraft::getCmatPos(const unsigned int cmat_itr) {
    if (isMyCmat(cmat_itr)) {
        return capacity_matrix_relation[cmat_itr];
    } else {
        throw std::invalid_argument("Invalid Cmat number passed to Spacecraft::getCmatPos().");
    }
}

auto Spacecraft::getTotalCharge(const tdArray& rho) const {
    double total_charge = 0.0;

    for(size_t cmat_itr = 0; cmat_itr < num_cmat; ++cmat_itr) {
        if (isMyCmat(cmat_itr)) {
            const auto& pos = capacity_matrix_relation.at(cmat_itr);
            total_charge += rho[pos.i][pos.j][pos.k];
        }
    }
    total_charge = MPIw::Environment::Comms[name].sum(total_charge);

    return total_charge;
}

void Spacecraft::makeCmatrixInvert(void) {
    //! B行列 -> C行列に変換
    Utils::makeInvert(capacity_matrix);

    total_cmat_value = 0.0;
    //! C_ij の sum を計算して保存しておく
    for(size_t col = 0; col < num_cmat; ++col) {
        for(size_t row = 0; row < num_cmat; ++row) {
            total_cmat_value += capacity_matrix(col, row);
        }
    }
}

bool Spacecraft::isIncluded(const Particle& p) const {
    if (!is_defined_in_this_process) return false;

    const auto pos = p.getPosition();

    return  object_map[pos.i    ][pos.j    ][pos.k    ] &&
            object_map[pos.i + 1][pos.j    ][pos.k    ] &&
            object_map[pos.i    ][pos.j + 1][pos.k    ] &&
            object_map[pos.i + 1][pos.j + 1][pos.k    ] &&
            object_map[pos.i    ][pos.j    ][pos.k + 1] &&
            object_map[pos.i + 1][pos.j    ][pos.k + 1] &&
            object_map[pos.i    ][pos.j + 1][pos.k + 1] &&
            object_map[pos.i + 1][pos.j + 1][pos.k + 1];
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
    for(const auto& one_node : capacity_matrix_relation) {
        const auto& pos = one_node.second;
        rho[pos.i][pos.j][pos.k] += charge_map[pos.i][pos.j][pos.k];
    }
}

void Spacecraft::redistributeCharge(tdArray& rho, const tdArray& phi) {
#ifdef CHARGE_CONSERVATION
    cout << "charge before redist: " << getTotalCharge(rho) << endl;
#endif

    double capacity_times_phi = 0.0;
    //! relationの中には元々内部ノードのPositionしか保存されていないので、
    //! 毎回判定しなくてよい
    for(const auto& one_node : capacity_matrix_relation) {
        const auto j = one_node.first;
        const auto& pos = one_node.second;
        for(size_t i = 0; i < num_cmat; ++i) {
            capacity_times_phi += capacity_matrix(i, j) * phi[pos.i][pos.j][pos.k];
        }
    }
    capacity_times_phi = MPIw::Environment::Comms[name].sum(capacity_times_phi);

    if (potential_fix != 0.0) {
        potential = potential_fix;
    } else {
        potential = capacity_times_phi / total_cmat_value;
    }
    cout << "[" << name << "] potential = " << Normalizer::unnormalizePotential(potential) << " V. " << endl;

    for(unsigned int i = 0; i < num_cmat; ++i) {
        double delta_rho = 0.0;

        for(unsigned int j = 0; j < num_cmat; ++j) {
            if (isMyCmat(j)) {
                const auto& pos = capacity_matrix_relation.at(j);
                delta_rho += capacity_matrix(i, j) * (potential - phi[pos.i][pos.j][pos.k]);
            }
        }
        MPIw::Environment::Comms[name].sum(delta_rho);

        if (isMyCmat(i)) {
            const auto& target_pos = capacity_matrix_relation.at(i);
            rho[target_pos.i][target_pos.j][target_pos.k] += delta_rho;
        }
    }

#ifdef CHARGE_CONSERVATION
    cout << "charge after redist: " << getTotalCharge(rho) << endl;
#endif
}

std::ostream& operator<<(std::ostream& ost, const Spacecraft& spc) {
    ost << "    name:" << spc.name << endl;
    ost << "    cmat:" << spc.num_cmat << endl;

    return ost;
}
