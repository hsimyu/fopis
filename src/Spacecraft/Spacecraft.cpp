#include "spacecraft.hpp"
#include "particle.hpp"
#include "normalizer.hpp"
#include "utils.hpp"
#include <fstream>
#include <stdexcept>
#include <regex>
#include <cassert>
#include <Eigen/Dense>

#define USE_BOOST
#include <simple_vtk.hpp>

//! static 変数の実体
unsigned int Spacecraft::num_of_spacecraft = 0;
void Spacecraft::construct(const size_t nx, const size_t ny, const size_t nz, const ObjectInfo_t& obj_info, const ObjectNodes& nodes, const ObjectNodes& glue_nodes, const ObjectCells& cells, const ObjectNodeTextures& textures, const ObjectConnectivityList& clist) {
    //! このオブジェクトがプロセス内で有効かどうかを保存しておく
    is_defined_in_this_process = (nodes.size() > 0);
    ++num_of_spacecraft;
    potential = 0.0;
    total_charge = 0.0;
    plot_history_width = obj_info.history_width;
    file_name = Utils::extractFileName(obj_info.file_name);
    potential_fix = Normalizer::normalizePotential(obj_info.potential_fix);

    if (is_defined_in_this_process) {
        connected_list = clist;

        //! セル定義マップ
        ObjectDefinedMapBool::extent_gen objectBoolExtents;
        object_cell_map.resize(objectBoolExtents[nx + 1][ny + 1][nz + 1]);

        // セル定義マップを初期化
        for(int i = 0; i < nx + 2; ++i) {
            for (int j = 0; j < ny + 2; ++j) {
                for (int k = 0; k < nz + 2; ++k) {
                    if (i != nx + 1 && j != ny + 1 && k != nz + 1) object_cell_map[i][j][k] = false;
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

            capacity_matrix_relation.emplace(std::piecewise_construct, std::make_tuple(cmat_itr), std::make_tuple(i, j, k));
        }

        //! applyChargeで使うためにGlueノード上のCmatPositionも覚えておく
        for(const auto& node_pair : glue_nodes) {
            const auto cmat_itr = node_pair.first;
            const auto& node_pos = node_pair.second;
            const auto i = node_pos[0];
            const auto j = node_pos[1];
            const auto k = node_pos[2];

            glue_capacity_matrix_relation.emplace(std::piecewise_construct, std::make_tuple(cmat_itr), std::make_tuple(i, j, k));
        }

        //! キャパシタンス行列のサイズを物体サイズに変更
        capacity_matrix.resize(num_cmat, num_cmat);

        tdArray::extent_gen tdExtents;
        for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
            //! Node ベース, glue cell ありの電荷密度マップを生成
            charge_map.emplace_back(tdExtents[nx + 2][ny + 2][nz + 2]);

            //! 電流 in/out を記録する要素を初期化
            current.push_back(0.0);

            const std::string pname = Environment::getParticleType(pid)->getName();

            //! jsonから読み取ったObjectInfo_t内の放出粒子名の指定が定義と一致していなかった場合を排除する
            if (obj_info.emit_particle_info.count(pname) > 0) {
                auto& pinfo = obj_info.emit_particle_info.at(pname);

                //! 相対的な粒子位置を計算
                const auto& rel_pos = Environment::getRelativePositionOnRootWithoutGlue(
                    pinfo.emission_position[0],
                    pinfo.emission_position[1],
                    pinfo.emission_position[2]
                );

                LocalParticleEmissionInfo local_pinfo;
                local_pinfo.relative_emission_position = Position(rel_pos[0], rel_pos[1], rel_pos[2]);
                local_pinfo.emission_vector = {pinfo.emission_vector[0], pinfo.emission_vector[1], pinfo.emission_vector[2]};
                emit_particle_info[ pid ] = local_pinfo;
            }
        }

        //! 電荷配列初期化
        //! @note: 物体の電荷配列マップは総和用の要素(index = 0)を持たない
        for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
            for(size_t i = 0; i < nx + 2; ++i) {
                for (size_t j = 0; j < ny + 2; ++j) {
                    for (size_t k = 0; k < nz + 2; ++k) {
                        charge_map[pid][i][j][k] = 0.0;
                    }
                }
            }
        }

        //! セルベースの物体定義 (物体内部判定用)
        for(const auto& cell_pos : cells) {
            object_cell_map[cell_pos[0]][cell_pos[1]][cell_pos[2]] = true;
        }

        //! 誘電体計算用の静電容量値をノードごとに計算
        if ( this->isDielectricSurface() ) {
            for (const auto& one_node : capacity_matrix_relation) {
                const auto& cmat_number = one_node.first;
                const auto& pos = one_node.second;
                const auto& texture = textures.at(cmat_number);

                //! ノード周りのFaceに割り当てられたTextureの平均値から計算
                double capacitance = 0.0;

                for(const auto& texture_indicies : texture) {
                    if (material_capacitances.count( texture_indicies ) == 0) {
                        //! マッピングに使う物理定数の値を property_list から取り出しておく
                        if (material_property_list.count( material_names[ texture_indicies ] ) > 0 ) {
                            //! 各ノードの capacitance = 真空の誘電率 * 比誘電率でいい?
                            material_capacitances[ texture_indicies ] = Normalizer::normalizeCapacitance(
                                eps0 * material_property_list.at( material_names[ texture_indicies ] ).at("RelativePermittivity")
                            );
                        } else {
                            //! material名の指定がおかしかった場合は落とす
                            std::string error_message = (format("Invalid material name %s is specified to be assigned the texture index %d.") % material_names[ texture_indicies ] % texture_indicies).str();

                            std::cerr << "[ERROR] " << error_message << endl;
                            throw std::invalid_argument(error_message);
                        }
                    }

                    capacitance += material_capacitances[ texture_indicies ];
                }

                capacitance_map[ cmat_number ] = capacitance / texture.size();
            }
        }
    }
}

void Spacecraft::saveWholeNodePositions(const ObjectNodes& whole_nodes) {
    //! RootNodeのみ、全体のノード位置を出力用に保持する
    //! 判定はGrid側で既に終わっているので必要なし
    for(const auto& node_pair : whole_nodes) {
        const auto cmat_itr = node_pair.first;
        const auto& node_pos = node_pair.second;
        const auto i = node_pos[0];
        const auto j = node_pos[1];
        const auto k = node_pos[2];

        whole_capacity_matrix_relation.emplace(std::piecewise_construct, std::make_tuple(cmat_itr), std::make_tuple(i, j, k));
    }
}

Position Spacecraft::getCmatPos(const unsigned int cmat_itr) {
    if (isMyCmat(cmat_itr)) {
        return capacity_matrix_relation[cmat_itr];
    } else {
        throw std::invalid_argument("Invalid Cmat number passed to Spacecraft::getCmatPos().");
    }
}

auto Spacecraft::getTotalCharge(const RhoArray& rho) const {
    double q = 0.0;

    for(size_t cmat_itr = 0; cmat_itr < num_cmat; ++cmat_itr) {
        if (isMyCmat(cmat_itr)) {
            const auto& pos = capacity_matrix_relation.at(cmat_itr);
            q += rho[0][pos.i][pos.j][pos.k];
        }
    }
    q = MPIw::Environment::Comms[name].sum(q);

    return q;
}

void Spacecraft::sumCurrent() {
    for(auto& v : current) {
        v = MPIw::Environment::Comms[name].sum(v);
    }
}

void Spacecraft::resetCurrent() {
    for(auto& v : current) {
        v = 0.0;
    }
}

void Spacecraft::updateTotalCmatValue() {
    total_cmat_value = 0.0;
    //! C_ij の sum を計算して保存しておく
    for(size_t col = 0; col < num_cmat; ++col) {
        for(size_t row = 0; row < num_cmat; ++row) {
            total_cmat_value += capacity_matrix(col, row);
        }
    }
}

void Spacecraft::makeCmatrixInvert(void) {
    //! B行列 -> C行列に変換
    capacity_matrix = capacity_matrix.inverse();
    this->updateTotalCmatValue();
}

bool Spacecraft::isContaining(const Particle& p) const {
    if (!is_defined_in_this_process) return false;

    return this->isContaining(p.getPosition());
}

inline bool Spacecraft::isContaining(const Position& pos) const {
    if (!is_defined_in_this_process) return false;

    if (Environment::isValidPosition(pos)) {
        return object_cell_map[pos.i][pos.j][pos.k];
    } else {
        return false;
    }
}

void Spacecraft::removeInnerParticle(Particle& p) const {
    if (isContaining(p)) { p.makeInvalid(); }
}

void Spacecraft::distributeInnerParticleCharge(Particle& p) {
    if (isContaining(p)) {
        const auto id = p.typeId;
        const auto q = p.getChargeOfSuperParticle();
        auto pos = p.getPosition();
        const auto old_pos = p.getOldPosition();

        bool move_along_x = (pos.i != old_pos.i);
        bool move_along_y = (pos.j != old_pos.j);
        bool move_along_z = (pos.k != old_pos.k);

        if (move_along_x && move_along_y && move_along_z) {
            const double mvx = fabs(p.vx);
            const double mvy = fabs(p.vy);
            const double mvz = fabs(p.vz);

            if (mvx >= mvy) {
                if (mvx >= mvz) {
                    const int di = (p.vx > 0.0) ? 0 : 1;
                    const int dj = (p.vy > 0.0) ? -1 : 1;
                    const int dk = (p.vz > 0.0) ? -1 : 1;
                    pos.setIJK(pos.i + di, pos.j + dj, pos.k + dk);
                    this->distributeInnerParticleChargeToXsurface(pos, id, q);
                } else {
                    const int di = (p.vx > 0.0) ? -1 : 1;
                    const int dj = (p.vy > 0.0) ? -1 : 1;
                    const int dk = (p.vz > 0.0) ? 0 : 1;
                    pos.setIJK(pos.i + di, pos.j + dj, pos.k + dk);
                    this->distributeInnerParticleChargeToZsurface(pos, id, q);
                }
            } else {
                if (mvy >= mvz) {
                    const int di = (p.vx > 0.0) ? -1 : 1;
                    const int dj = (p.vy > 0.0) ? 0 : 1;
                    const int dk = (p.vz > 0.0) ? -1 : 1;
                    pos.setIJK(pos.i + di, pos.j + dj, pos.k + dk);
                    this->distributeInnerParticleChargeToYsurface(pos, id, q);
                } else {
                    const int di = (p.vx > 0.0) ? -1 : 1;
                    const int dj = (p.vy > 0.0) ? -1 : 1;
                    const int dk = (p.vz > 0.0) ? 0 : 1;
                    pos.setIJK(pos.i + di, pos.j + dj, pos.k + dk);
                    this->distributeInnerParticleChargeToZsurface(pos, id, q);
                }
            }
        } else if (move_along_x && move_along_y) {
            if (fabs(p.vx) >= fabs(p.vy)) {
                const int di = (p.vx > 0.0) ? 0 : 1;
                const int dj = (p.vy > 0.0) ? -1 : 1;
                pos.setIJK(pos.i + di, pos.j + dj, pos.k);
                this->distributeInnerParticleChargeToXsurface(pos, id, q);
            } else {
                const int di = (p.vx > 0.0) ? -1 : 1;
                const int dj = (p.vy > 0.0) ? 0 : 1;
                pos.setIJK(pos.i + di, pos.j + dj, pos.k);
                this->distributeInnerParticleChargeToYsurface(pos, id, q);
            }
        } else if (move_along_x && move_along_z) {
            if (fabs(p.vx) >= fabs(p.vz)) {
                const int di = (p.vx > 0.0) ? 0 : 1;
                const int dk = (p.vz > 0.0) ? -1 : 1;
                pos.setIJK(pos.i + di, pos.j, pos.k + dk);
                this->distributeInnerParticleChargeToXsurface(pos, id, q);
            } else {
                const int di = (p.vx > 0.0) ? -1 : 1;
                const int dk = (p.vz > 0.0) ? 0 : 1;
                pos.setIJK(pos.i + di, pos.j, pos.k + dk);
                this->distributeInnerParticleChargeToZsurface(pos, id, q);
            }
        } else if (move_along_y && move_along_z) {
            if (fabs(p.vy) >= fabs(p.vz)) {
                const int dj = (p.vy > 0.0) ? 0 : 1;
                const int dk = (p.vz > 0.0) ? -1 : 1;
                pos.setIJK(pos.i, pos.j + dj, pos.k + dk);
                this->distributeInnerParticleChargeToYsurface(pos, id, q);
            } else {
                const int dj = (p.vy > 0.0) ? -1 : 1;
                const int dk = (p.vz > 0.0) ? 0 : 1;
                pos.setIJK(pos.i, pos.j + dj, pos.k + dk);
                this->distributeInnerParticleChargeToZsurface(pos, id, q);
            }
        } else if (move_along_x) {
            const int di = (p.vx > 0.0) ? 0 : 1;
            pos.setIJK(pos.i + di, pos.j, pos.k);
            this->distributeInnerParticleChargeToXsurface(pos, id, q);
        } else if (move_along_y) {
            const int dj = (p.vy > 0.0) ? 0 : 1;
            pos.setIJK(pos.i, pos.j + dj, pos.k);
            this->distributeInnerParticleChargeToYsurface(pos, id, q);
        } else if (move_along_z) {
            const int dk = (p.vz > 0.0) ? 0 : 1;
            pos.setIJK(pos.i, pos.j, pos.k + dk);
            this->distributeInnerParticleChargeToZsurface(pos, id, q);
        } else {
            throw std::logic_error("[ERROR] An Invalid Logic on Particle Collection!!!!!");
        }

        current[id] += q;
        p.makeInvalid();
    }
}

inline void Spacecraft::distributeInnerParticleChargeToXsurface(const Position& pos, const int id, const double charge) {
    charge_map[id][pos.i][pos.j    ][pos.k    ] += charge * pos.dy2 * pos.dz2;
    charge_map[id][pos.i][pos.j + 1][pos.k    ] += charge * pos.dy1 * pos.dz2;
    charge_map[id][pos.i][pos.j    ][pos.k + 1] += charge * pos.dy2 * pos.dz1;
    charge_map[id][pos.i][pos.j + 1][pos.k + 1] += charge * pos.dy1 * pos.dz1;
}

inline void Spacecraft::distributeInnerParticleChargeToYsurface(const Position& pos, const int id, const double charge) {
    charge_map[id][pos.i    ][pos.j][pos.k    ] += charge * pos.dx2 * pos.dz2;
    charge_map[id][pos.i + 1][pos.j][pos.k    ] += charge * pos.dx1 * pos.dz2;
    charge_map[id][pos.i    ][pos.j][pos.k + 1] += charge * pos.dx2 * pos.dz1;
    charge_map[id][pos.i + 1][pos.j][pos.k + 1] += charge * pos.dx1 * pos.dz1;
}

inline void Spacecraft::distributeInnerParticleChargeToZsurface(const Position& pos, const int id, const double charge) {
    charge_map[id][pos.i    ][pos.j    ][pos.k] += charge * pos.dy2 * pos.dx2;
    charge_map[id][pos.i    ][pos.j + 1][pos.k] += charge * pos.dy1 * pos.dx2;
    charge_map[id][pos.i + 1][pos.j    ][pos.k] += charge * pos.dy2 * pos.dx1;
    charge_map[id][pos.i + 1][pos.j + 1][pos.k] += charge * pos.dy1 * pos.dx1;
}

void Spacecraft::applyCharge(RhoArray& rho) const {
    //! 電荷を場に印加する
    for(const auto& node : capacity_matrix_relation) {
        const auto& pos = node.second;
        for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
            rho[pid + 1][pos.i][pos.j][pos.k] += charge_map[pid][pos.i][pos.j][pos.k];
        }
    }

    //! Glueノード上の電荷
    for(const auto& node : glue_capacity_matrix_relation) {
        const auto& pos = node.second;
        for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
            rho[pid + 1][pos.i][pos.j][pos.k] += charge_map[pid][pos.i][pos.j][pos.k];
        }
    }
}

void Spacecraft::redistributeCharge(RhoArray& rho, const tdArray& phi) {
    auto q = getTotalCharge(rho);

    if (MPIw::Environment::isRootNode(name)) {
        cout << format("[INFO] [%s] charge before redist: %16.7e") % name % q << endl;
    }

    double capacity_times_phi = 0.0;
    //! relationの中には元々内部ノードのPositionしか保存されていないので、
    //! 毎回判定しなくてよい

    if (this->isDielectricSurface()) {
        //! 誘電体の場合
        for(const auto& one_node : capacity_matrix_relation) {
            const auto j = one_node.first;
            const auto& pos = one_node.second;
            const auto capacitance = capacitance_map[j];

            if (capacitance > 0.0) {
                for(size_t i = 0; i < num_cmat; ++i) {
                    double node_charge = 0.0;
                    for (int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
                        node_charge += charge_map[pid][pos.i][pos.j][pos.k];
                    }
                    capacity_times_phi += capacity_matrix(i, j) * (phi[pos.i][pos.j][pos.k] - node_charge / capacitance);
                }
            } else {
                for(size_t i = 0; i < num_cmat; ++i) {
                    capacity_times_phi += capacity_matrix(i, j) * phi[pos.i][pos.j][pos.k];
                }
            }
        }
    } else {
        //! 完全導体の場合
        for(const auto& one_node : capacity_matrix_relation) {
            const auto j = one_node.first;
            const auto& pos = one_node.second;
            for(size_t i = 0; i < num_cmat; ++i) {
                capacity_times_phi += capacity_matrix(i, j) * phi[pos.i][pos.j][pos.k];
            }
        }
    }

    capacity_times_phi = MPIw::Environment::Comms[name].sum(capacity_times_phi);

    if (MPIw::Environment::isRootNode(name)) {
        cout << format("[INFO] [%s] Capacity x Potential: %16.7e") % name % capacity_times_phi << endl;
    }

    if (potential_fix != 0.0) {
        potential = potential_fix;
    } else {
        potential = (capacity_times_phi + q) / total_cmat_value;
    }

    if (MPIw::Environment::isRootNode(name)) {
        cout << format("[INFO] [%s] potential = %s V") % name % Normalizer::unnormalizePotential(potential) << endl;
    }

    if (this->isDielectricSurface()) {
        //! 誘電体の場合
        double sum_delta_rho = 0.0;
        for(unsigned int i = 0; i < num_cmat; ++i) {
            double delta_rho = 0.0;

            for(unsigned int j = 0; j < num_cmat; ++j) {
                if (isMyCmat(j)) {
                    const auto capacitance = capacitance_map[j];
                    const auto& pos = capacity_matrix_relation.at(j);

                    if (capacitance > 0.0) {
                        double node_charge = 0.0;
                        for (int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
                            node_charge += charge_map[pid][pos.i][pos.j][pos.k];
                        }
                        delta_rho += capacity_matrix(i, j) * (potential - phi[pos.i][pos.j][pos.k] + node_charge / capacitance_map[j]);
                    } else {
                        delta_rho += capacity_matrix(i, j) * (potential - phi[pos.i][pos.j][pos.k]);
                    }
                }
            }
            delta_rho = MPIw::Environment::Comms[name].sum(delta_rho);
            sum_delta_rho += delta_rho;

            if (isMyCmat(i)) {
                const auto& target_pos = capacity_matrix_relation.at(i);
                rho[0][target_pos.i][target_pos.j][target_pos.k] += delta_rho;
            }
        }

        if (MPIw::Environment::isRootNode(name)) {
            cout << format("[INFO] [%s] sum of redist charge: %16.7e") % name % sum_delta_rho << endl;
        }
    } else {
        //! 完全導体の場合
        double sum_delta_rho = 0.0;
        for(unsigned int i = 0; i < num_cmat; ++i) {
            double delta_rho = 0.0;

            for(unsigned int j = 0; j < num_cmat; ++j) {
                if (isMyCmat(j)) {
                    const auto& pos = capacity_matrix_relation.at(j);
                    delta_rho += capacity_matrix(i, j) * (potential - phi[pos.i][pos.j][pos.k]);
                }
            }
            delta_rho = MPIw::Environment::Comms[name].sum(delta_rho);
            sum_delta_rho += delta_rho;

            if (isMyCmat(i)) {
                const auto& target_pos = capacity_matrix_relation.at(i);
                rho[0][target_pos.i][target_pos.j][target_pos.k] += delta_rho;
            }
        }

        if (MPIw::Environment::isRootNode(name)) {
            cout << format("[INFO] [%s] sum of redist charge: %16.7e") % name % sum_delta_rho << endl;
        }
    }

    total_charge = getTotalCharge(rho);
    // update total current for output
    sumCurrent();
}

//! 粒子放出関連
void Spacecraft::emitParticles(ParticleArray& parray) {
    for(const auto& pinfo : emit_particle_info) {
        const auto id = pinfo.first;
        const auto& info = pinfo.second;

        if (Environment::isValidPosition(info.relative_emission_position)) {
            const auto emit_ptype_ptr = Environment::getEmissionParticleType(id);
            const auto max_amount = emit_ptype_ptr->getEmissionAmount();

            //! 放出時に呼び出すメンバ関数ポインタ
            std::function<void(Position&, const int, const double)> emission_func;

            //! 放出時の座標修正子
            if (fabs(info.emission_vector[0]) == 1.0) {
                emission_func = [this, &rel_pos = info.relative_emission_position](Position& pos, const int id, const double charge) {
                    this->subtractChargeOfParticleFromXsurface(rel_pos, pos, id, charge);
                };
            } else if (fabs(info.emission_vector[1]) == 1.0) {
                emission_func = [this, &rel_pos = info.relative_emission_position](Position& pos, const int id, const double charge) {
                    this->subtractChargeOfParticleFromYsurface(rel_pos, pos, id, charge);
                };
            } else if (fabs(info.emission_vector[2]) == 1.0) {
                emission_func = [this, &rel_pos = info.relative_emission_position](Position& pos, const int id, const double charge) {
                    this->subtractChargeOfParticleFromZsurface(rel_pos, pos, id, charge);
                };
            } else {
                std::string error_message = (format("[ERROR] At Spacecraft::emitParticles: Now 'emission_vector' must be [+-1,0,0] or [0,+-1,0] or [0,0,+-1].")).str();
                throw std::invalid_argument(error_message);
            }

            if (emit_ptype_ptr->getType() == "beam") {
                const auto beam_ptype_ptr = Environment::getBeamParticleType(id);
                const auto charge = beam_ptype_ptr->getChargeOfSuperParticle();

                for(int i = 0; i < max_amount; ++i) {
                    Particle p = beam_ptype_ptr->generateNewParticle(info.relative_emission_position, info.emission_vector);

                    if (this->isValidEmission(p)) {
                        auto pos = p.getPosition();
                        emission_func(pos, id, charge);
                        parray[id].push_back( std::move(p) );
                    }
                }
            }
        }
    }
}

inline void Spacecraft::subtractChargeOfParticleFromXsurface(const Position& emission_pos, const Position& pos, const int id, const double charge) {
    charge_map[id][emission_pos.i][emission_pos.j    ][emission_pos.k    ] -= charge * pos.dy2 * pos.dz2;
    charge_map[id][emission_pos.i][emission_pos.j + 1][emission_pos.k    ] -= charge * pos.dy1 * pos.dz2;
    charge_map[id][emission_pos.i][emission_pos.j    ][emission_pos.k + 1] -= charge * pos.dy2 * pos.dz1;
    charge_map[id][emission_pos.i][emission_pos.j + 1][emission_pos.k + 1] -= charge * pos.dy1 * pos.dz1;

    current[id] -= charge;
}

inline void Spacecraft::subtractChargeOfParticleFromYsurface(const Position& emission_pos, const Position& pos, const int id, const double charge) {
    charge_map[id][emission_pos.i    ][emission_pos.j    ][emission_pos.k    ] -= charge * pos.dx2 * pos.dz2;
    charge_map[id][emission_pos.i + 1][emission_pos.j    ][emission_pos.k    ] -= charge * pos.dx1 * pos.dz2;
    charge_map[id][emission_pos.i    ][emission_pos.j    ][emission_pos.k + 1] -= charge * pos.dx2 * pos.dz1;
    charge_map[id][emission_pos.i + 1][emission_pos.j    ][emission_pos.k + 1] -= charge * pos.dx1 * pos.dz1;

    current[id] -= charge;
}

inline void Spacecraft::subtractChargeOfParticleFromZsurface(const Position& emission_pos, const Position& pos, const int id, const double charge) {
    charge_map[id][emission_pos.i    ][emission_pos.j    ][emission_pos.k    ] -= charge * pos.dx2 * pos.dy2;
    charge_map[id][emission_pos.i + 1][emission_pos.j    ][emission_pos.k    ] -= charge * pos.dx1 * pos.dy2;
    charge_map[id][emission_pos.i    ][emission_pos.j + 1][emission_pos.k    ] -= charge * pos.dx2 * pos.dy1;
    charge_map[id][emission_pos.i + 1][emission_pos.j + 1][emission_pos.k    ] -= charge * pos.dx1 * pos.dy1;

    current[id] -= charge;
}

inline bool Spacecraft::isValidEmission(Particle& p) const {
    //! 放出前は内部、放出後は外部に入れば valid と見なす
    return ( (isContaining(p)) && !isContaining(p.getNewPosition()) );
}

//! I/O関連
std::string Spacecraft::getLogHeader() const {
    std::string format_string = "%16s %16s";
    for(int i = 0; i < Environment::num_of_particle_types; ++i) {
        format_string += " %16s";
    }

    auto format_base = format(std::move(format_string));
    format_base = format_base % "Potential [V]";
    format_base = format_base % "Charge [C]";
    for(int i = 0; i < Environment::num_of_particle_types; ++i) {
        format_base = format_base % (Environment::getParticleType(i)->getName() + " [A]");
    }

    std::string header = format_base.str();
    return header;
}

std::string Spacecraft::getLogEntry() const {
    std::string format_string = "%16.7e %16.7e";
    for(int i = 0; i < Environment::num_of_particle_types; ++i) {
        format_string += " %16.7e";
    }

    auto format_base = format(std::move(format_string));
    format_base = format_base % Normalizer::unnormalizePotential(potential);
    format_base = format_base % Normalizer::unnormalizeCharge(total_charge);
    for(int i = 0; i < Environment::num_of_particle_types; ++i) {
        format_base = format_base % Normalizer::unnormalizeCurrent(current[i]);
    }

    std::string entry = format_base.str();
    return entry;
}

std::ostream& operator<<(std::ostream& ost, const Spacecraft& spc) {
    ost << "  name: " << spc.name << endl;
    ost << "  object file name: " << spc.file_name << endl;
    ost << "  cmat: " << spc.num_cmat << endl;
    ost << "  surface_type: " << spc.surface_type << endl;

    if (spc.isDielectricSurface()) {
        for(const auto& mt : spc.material_names) {
            ost << "    index " << mt.first << ": " << mt.second << endl;
        }
    }

    ost << "  emit particles: \n";
    for(const auto& pair : spc.emit_particle_info) {
        ost << "    " << Environment::getParticleType(pair.first)->getName() << ": \n";

        const auto& local_pinfo = pair.second;
        ost << "      relative_emission_position: " <<
            format("%s %s %s\n") %
            local_pinfo.relative_emission_position.i %
            local_pinfo.relative_emission_position.j %
            local_pinfo.relative_emission_position.k;
        ost << "      emission_vector: " <<
            format("%s %s %s") %
            local_pinfo.emission_vector[0] %
            local_pinfo.emission_vector[1] %
            local_pinfo.emission_vector[2] << endl;
    }

    return ost;
}

void Spacecraft::insertConnectivity(SimpleVTK& gen) const {
    for(const auto& vect : connected_list) {
        gen.addVector(vect);
    }
}

void Spacecraft::insertPoints(SimpleVTK& gen) const {
    const auto unnorm = Normalizer::unnormalizeLength(1.0);
    for(size_t cmat_itr = 0; cmat_itr < num_cmat; ++cmat_itr) {
        const auto& pos = whole_capacity_matrix_relation.at(cmat_itr);
        gen.addItem(unnorm * pos.i, unnorm * pos.j, unnorm * pos.k);
    }
}

void Spacecraft::insertOffsets(SimpleVTK& gen) const {
    std::vector<int> offsets(connected_list.size());

    for(int i = 0; i < connected_list.size(); ++i) {
        offsets[i] = 4 * (i + 1);
    }

    gen.addVector(offsets);
}

void Spacecraft::insertTypes(SimpleVTK& gen) const {
    std::vector<int> types(connected_list.size(), 8);

    gen.addVector(types);
}

template<typename T>
std::vector<T> Spacecraft::getWholePotentialMap(const tdArray& phi) const {
    std::vector<T> potential_values(num_cmat);
    const auto unnorm = Normalizer::unnormalizePotential(1.0);

    for(unsigned int cmat_itr = 0; cmat_itr < num_cmat; ++cmat_itr) {
        if (isMyCmat(cmat_itr)) {
            const auto& pos = capacity_matrix_relation.at(cmat_itr);
            potential_values[cmat_itr] = static_cast<T>(unnorm * phi[pos.i][pos.j][pos.k]);
        }
    }

    auto result = MPIw::Environment::Comms[name].sum(potential_values);

    return result;
}

void Spacecraft::plotPotentialMapping(const int timestep, const tdArray& phi) const {
    if (!isDefined()) return;

    //! 全体の電位マップを集める
    std::vector<float> potential_values = this->getWholePotentialMap<float>(phi);

    if (MPIw::Environment::isRootNode(name)) {
        SimpleVTK gen;
        gen.beginVTK("UnstructuredGrid");
            gen.beginContent();
                gen.beginPiece();
                gen.setNumberOfPoints(num_cmat);
                gen.setNumberOfCells(connected_list.size());
                    gen.beginPointData();
                    gen.setScalars("potential");
                        gen.beginDataArray("potential", "Float32", "ascii");
                            gen.addVector(potential_values);
                        gen.endDataArray();
                    gen.endPointData();

                    gen.beginPoints();
                        gen.beginDataArray("Points", "Float32", "ascii");
                        gen.setNumberOfComponents("3");
                            insertPoints(gen);
                        gen.endDataArray();
                    gen.endPoints();

                    gen.beginCells();
                        gen.beginDataArray("connectivity", "Int32", "ascii");
                            this->insertConnectivity(gen);
                        gen.endDataArray();
                        gen.beginDataArray("offsets", "Int32", "ascii");
                            this->insertOffsets(gen);
                        gen.endDataArray();
                        gen.beginDataArray("types", "UInt8", "ascii");
                            this->insertTypes(gen);
                        gen.endDataArray();
                    gen.endCells();
                gen.endPiece();
            gen.endContent();
        gen.endVTK();

        constexpr char* filepath_header = "data/";
        gen.generate(filepath_header + name + "_potential_mapping_" + std::to_string(timestep));
    }
}

// Utility Functions for Objects
namespace ObjectUtils {
    ObjectDataFromFile getObjectNodesFromObjFile(const std::string& obj_file_name) {
        ObjectDataFromFile obj_data;

        if (Utils::isExistingFile(obj_file_name)) {
            const auto nx = Environment::nx;
            const auto ny = Environment::ny;
            const auto nz = Environment::nz;
            ObjectDefinedMapInt object_node_map(boost::extents[nx][ny][nz]);

            using ObjectFaceMap = boost::multi_array<bool, 3>;
            ObjectFaceMap object_xface_map(boost::extents[nx][ny - 1][nz - 1]);
            ObjectFaceMap object_yface_map(boost::extents[nx - 1][ny][nz - 1]);
            ObjectFaceMap object_zface_map(boost::extents[nx - 1][ny - 1][nz]);

            //! 初期化
            for(int i = 0; i < nx; ++i) {
                for (int j = 0; j < ny; ++j) {
                    for (int k = 0; k < nz; ++k) {
                        object_node_map[i][j][k] = -1;

                        if (j != ny - 1 && k != nz - 1) object_xface_map[i][j][k] = false;
                        if (i != nx - 1 && k != nz - 1) object_yface_map[i][j][k] = false;
                        if (i != nx - 1 && j != ny - 1) object_zface_map[i][j][k] = false;
                    }
                }
            }

            std::ifstream file_input(obj_file_name, std::ios::in);
            std::string buffer;

            //! vとf 始まりの文字列にマッチするパターン
            std::regex re_vertex("v (.*)");
            std::regex re_face("f (.*)");

            //! vertexとfaceを一時格納しておくarray
            using Vertex = std::array<int, 3>;
            using Face = std::array<int, 5>;
            std::vector<Vertex> vertices;
            std::vector<Face> faces;
            
            std::map<int, int> texture_counts{};

            while (!file_input.eof()) {
                std::getline(file_input, buffer);
                if (std::regex_match(buffer, re_vertex)) {
                    const auto& verts = Utils::split(buffer, ' ');
                    Vertex vertex;

                    vertex[0] = static_cast<int>(std::stoi(verts[1]) + Environment::nx / 2);
                    vertex[1] = static_cast<int>(std::stoi(verts[2]) + Environment::ny / 2);
                    vertex[2] = static_cast<int>(std::stoi(verts[3]) + Environment::nz / 2);

                    vertices.push_back( std::move(vertex) );
                } else if (std::regex_match(buffer, re_face)) {
                    const auto& fs = Utils::split(buffer, ' ');
                    Face face_numbers;

                    for(int i = 1; i < 5; ++i) {
                        const auto& splitted_face = Utils::split(fs[i], '/');
                        face_numbers[i - 1] = static_cast<int>(std::stoi(splitted_face[0]));
                    }

                    //! 5番目の要素は texture_index
                    {
                        const auto& splitted_face = Utils::split(fs[1], '/');
                        face_numbers[4] = static_cast<int>(std::stoi(splitted_face[1]));

                        //! テクスチャが何回出てきたかカウントする
                        texture_counts[ face_numbers[4] ] += 1;

                    }

                    faces.push_back( std::move(face_numbers) );
                }
            }

            unsigned int num_cmat = 0;
            for(const auto& face : faces) {
                const auto& vert1 = vertices[ face[0] - 1 ];
                const auto& vert2 = vertices[ face[1] - 1 ];
                const auto& vert3 = vertices[ face[2] - 1 ];
                const auto& vert4 = vertices[ face[3] - 1 ];

                //! @note: faceの分割の第3番号(1/1/3 2/1/3 3/1/3 4/1/3 の3)が法線ベクトルに対応するので
                //! そっちで判定した方が良いか
                if (vert1[0] == vert2[0] && vert1[0] == vert3[0] && vert1[0] == vert4[0]) {
                    //! X面
                    const int i = vert1[0];
                    const int miny = std::min({vert1[1], vert2[1], vert3[1], vert4[1]});
                    const int maxy = std::max({vert1[1], vert2[1], vert3[1], vert4[1]});
                    const int minz = std::min({vert1[2], vert2[2], vert3[2], vert4[2]});
                    const int maxz = std::max({vert1[2], vert2[2], vert3[2], vert4[2]});

                    for(int j = miny; j < maxy + 1; ++j) {
                        for(int k = minz; k < maxz + 1; ++k) {
                            if (object_node_map[i][j][k] < 0) {
                                obj_data.nodes[num_cmat] = {{i, j, k}};
                                obj_data.textures[num_cmat] = {face[4]};
                                object_node_map[i][j][k] = num_cmat;
                                ++num_cmat;
                            } else {
                                obj_data.textures[ object_node_map[i][j][k] ].push_back(face[4]);
                            }

                            if (j != maxy && k != maxz) {
                                object_xface_map[i][j][k] = true;
                            }
                        }
                    }
                } else if (vert1[1] == vert2[1] && vert1[1] == vert3[1] && vert1[1] == vert4[1]) {
                    //! Y面
                    const int j = vert1[1];
                    const int minx = std::min({vert1[0], vert2[0], vert3[0], vert4[0]});
                    const int maxx = std::max({vert1[0], vert2[0], vert3[0], vert4[0]});
                    const int minz = std::min({vert1[2], vert2[2], vert3[2], vert4[2]});
                    const int maxz = std::max({vert1[2], vert2[2], vert3[2], vert4[2]});

                    for(int i = minx; i < maxx + 1; ++i) {
                        for(int k = minz; k < maxz + 1; ++k) {
                            if (object_node_map[i][j][k] < 0) {
                                obj_data.nodes[num_cmat] = {{i, j, k}};
                                obj_data.textures[num_cmat] = {face[4]};
                                object_node_map[i][j][k] = num_cmat;
                                ++num_cmat;
                            } else {
                                obj_data.textures[ object_node_map[i][j][k] ].push_back(face[4]);
                            }

                            if (i != maxx && k != maxz) {
                                object_yface_map[i][j][k] = true;
                            }
                        }
                    }
                } else if (vert1[2] == vert2[2] && vert1[2] == vert3[2] && vert1[2] == vert4[2]) {
                    //! Z面
                    const int k = vert1[2];
                    const int minx = std::min({vert1[0], vert2[0], vert3[0], vert4[0]});
                    const int maxx = std::max({vert1[0], vert2[0], vert3[0], vert4[0]});
                    const int miny = std::min({vert1[1], vert2[1], vert3[1], vert4[1]});
                    const int maxy = std::max({vert1[1], vert2[1], vert3[1], vert4[1]});

                    for(int i = minx; i < maxx + 1; ++i) {
                        for(int j = miny; j < maxy + 1; ++j) {
                            if (object_node_map[i][j][k] < 0) {
                                obj_data.nodes[num_cmat] = {{i, j, k}};
                                obj_data.textures[num_cmat] = {face[4]};
                                object_node_map[i][j][k] = num_cmat;
                                ++num_cmat;
                            } else {
                                obj_data.textures[ object_node_map[i][j][k] ].push_back(face[4]);
                            }

                            if (i != maxx && j != maxy) {
                                object_zface_map[i][j][k] = true;
                            }
                        }
                    }
                } else {
                    throw std::logic_error("Face type cannot be determinted by vertices position.");
                }
            }

            //! テクスチャカウントを表示
            if (Environment::isRootNode) {
                for(const auto& pair : texture_counts) {
                    cout << format("  [OBJECT DEFINE INFO] texture index %s: %s faces") % pair.first % pair.second << endl;
                }
            }

            //! セルマップ定義のためにFaceを使ってカウントする
            using CellCountMap = boost::multi_array<int, 3>;
            CellCountMap object_cell_count_map(boost::extents[nx - 1][ny - 1][nz - 1]);

            for (int j = 0; j < ny - 1; ++j) {
                for (int k = 0; k < nz - 1; ++k) {
                    for(int lower_index = 0; lower_index < nx - 1; ++lower_index) {
                        //! 下端を見つける
                        if (object_xface_map[lower_index][j][k]) {

                            //! 下端が見つかった場合は上端を探す
                            for(int upper_index = lower_index + 1; upper_index < nx - 1; ++upper_index) {
                                if (object_xface_map[upper_index][j][k]) {
                                    //! 上端が見つかったら間を塗る
                                    for(int i = lower_index; i < upper_index; ++i) {
                                        object_cell_count_map[i][j][k] += 1;
                                    }

                                    //! 下端のindexを進めてループを抜ける
                                    lower_index = upper_index;
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            for(int i = 0; i < nx - 1; ++i) {
                for (int k = 0; k < nz - 1; ++k) {
                    for(int lower_index = 0; lower_index < ny - 1; ++lower_index) {
                        //! 下端を見つける
                        if (object_yface_map[i][lower_index][k]) {
                            //! 下端が見つかった場合は上端を探す
                            for(int upper_index = lower_index + 1; upper_index < ny - 1; ++upper_index) {
                                if (object_yface_map[i][upper_index][k]) {
                                    //! 上端が見つかったら間を塗る
                                    for(int j = lower_index; j < upper_index; ++j) {
                                        object_cell_count_map[i][j][k] += 1;
                                    }

                                    //! 下端のindexを進めてループを抜ける
                                    lower_index = upper_index;
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            for(int i = 0; i < nx - 1; ++i) {
                for (int j = 0; j < ny - 1; ++j) {
                    for(int lower_index = 0; lower_index < nz - 1; ++lower_index) {
                        //! 下端を見つける
                        if (object_zface_map[i][j][lower_index]) {
                            //! 下端が見つかった場合は上端を探す
                            for(int upper_index = lower_index + 1; upper_index < nz - 1; ++upper_index) {
                                if (object_zface_map[i][j][upper_index]) {
                                    //! 上端が見つかったら間を塗る
                                    for(int k = lower_index; k < upper_index; ++k) {
                                        object_cell_count_map[i][j][k] += 1;
                                    }

                                    //! 下端のindexを進めてループを抜ける
                                    lower_index = upper_index;
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            for(int i = 0; i < nx - 1; ++i) {
                for (int j = 0; j < ny - 1; ++j) {
                    for (int k = 0; k < nz - 1; ++k) {
                        //! countが3なら内部と定義
                        if (object_cell_count_map[i][j][k] == 3) {
                            obj_data.cells.push_back({{i, j, k}});
                        }
                    }
                }
            }

            //! Connectivity List の 計算
            for(unsigned int cmat_itr = 0; cmat_itr < num_cmat; ++cmat_itr) {
                const auto& pos = obj_data.nodes[cmat_itr];
                const auto i = pos[0];
                const auto j = pos[1];
                const auto k = pos[2];

                // x-face check
                if (object_cell_count_map[i][j][k] == 3 && object_cell_count_map[i - 1][j][k] != 3) {
                    std::vector<unsigned int> clist = {
                        cmat_itr,
                        static_cast<unsigned int>(object_node_map[i][j + 1][k   ]),
                        static_cast<unsigned int>(object_node_map[i][j    ][k + 1]),
                        static_cast<unsigned int>(object_node_map[i][j + 1][k + 1])
                    };
                    obj_data.connected_list.push_back(clist);
                } else if (object_cell_count_map[i][j][k] != 3 && object_cell_count_map[i - 1][j][k] == 3) {
                    std::vector<unsigned int> clist = {
                        cmat_itr,
                        static_cast<unsigned int>(object_node_map[i][j + 1][k    ]),
                        static_cast<unsigned int>(object_node_map[i][j    ][k + 1]),
                        static_cast<unsigned int>(object_node_map[i][j + 1][k + 1])
                    };
                    obj_data.connected_list.push_back(clist);
                }

                // y-face check
                if (object_cell_count_map[i][j][k] == 3 && object_cell_count_map[i][j - 1][k] != 3) {
                    std::vector<unsigned int> clist = {
                        cmat_itr,
                        static_cast<unsigned int>(object_node_map[i + 1][j][k    ]),
                        static_cast<unsigned int>(object_node_map[i    ][j][k + 1]),
                        static_cast<unsigned int>(object_node_map[i + 1][j][k + 1])
                    };
                    obj_data.connected_list.push_back(clist);
                } else if (object_cell_count_map[i][j][k] != 3 && object_cell_count_map[i][j - 1][k] == 3) {
                    std::vector<unsigned int> clist = {
                        cmat_itr,
                        static_cast<unsigned int>(object_node_map[i + 1][j][k    ]),
                        static_cast<unsigned int>(object_node_map[i    ][j][k + 1]),
                        static_cast<unsigned int>(object_node_map[i + 1][j][k + 1])
                    };
                    obj_data.connected_list.push_back(clist);
                }

                // z-face check
                if (object_cell_count_map[i][j][k] == 3 && object_cell_count_map[i][j][k - 1] != 3) {
                    std::vector<unsigned int> clist = {
                        cmat_itr,
                        static_cast<unsigned int>(object_node_map[i + 1][j    ][k]),
                        static_cast<unsigned int>(object_node_map[i    ][j + 1][k]),
                        static_cast<unsigned int>(object_node_map[i + 1][j + 1][k])
                    };
                    obj_data.connected_list.push_back(clist);
                } else if (object_cell_count_map[i][j][k] != 3 && object_cell_count_map[i][j][k - 1] == 3) {
                    std::vector<unsigned int> clist = {
                        cmat_itr,
                        static_cast<unsigned int>(object_node_map[i + 1][j    ][k]),
                        static_cast<unsigned int>(object_node_map[i    ][j + 1][k]),
                        static_cast<unsigned int>(object_node_map[i + 1][j + 1][k])
                    };
                    obj_data.connected_list.push_back(clist);
                }
            }
        } else {
            std::string error_message = (format("File %s does not exist.") % obj_file_name).str();
            throw std::invalid_argument(error_message);
        }

        return obj_data;
    }
}
