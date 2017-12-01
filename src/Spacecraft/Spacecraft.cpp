#include "spacecraft.hpp"
#include "particle.hpp"
#include "normalizer.hpp"
#include "grid.hpp"
#include "utils.hpp"
#include "mpiw.hpp"
#include <fstream>
#include <stdexcept>
#include <cassert>

#define USE_BOOST
#include <simple_vtk.hpp>

//! static 変数の実体
unsigned int Spacecraft::num_of_spacecraft = 0;
void Spacecraft::construct(
    const size_t nx, const size_t ny, const size_t nz, const ObjectInfo_t& obj_info,
    const ObjectNodes& nodes, const ObjectNodes& glue_nodes, const ObjectNodes& whole_nodes,
    const ObjectCells& cells, const ObjectNodeTextures& textures, const ObjectConnectivityList& clist) {

    //! このオブジェクトがプロセス内で有効かどうかを保存しておく
    is_defined_in_this_process = (nodes.size() > 0);

    ++num_of_spacecraft;
    potential = 0.0;
    total_charge = 0.0;
    plot_potential_mapping_width = obj_info.plot_potential_mapping_width;
    file_name = Utils::extractFileName(obj_info.file_name);

    is_potential_fixed = obj_info.is_potential_fixed;
    fixed_potential = Normalizer::normalizePotential(obj_info.fixed_potential);
    initial_potential_offset = Normalizer::normalizePotential(obj_info.initial_potential_offset);

    force_computation = obj_info.force_computation;

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

        //! whole_nodesも保持
        this->saveWholeNodePositions(whole_nodes);

        //! キャパシタンス行列のサイズを物体サイズに変更
        capacity_matrix.resize(num_cmat);

        tdArray::extent_gen tdExtents;
        for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {

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
        ObjectChargeMap::extent_gen cm_extent;
        //! @note: 物体の電荷配列マップは総和用の要素(index = 0)を持たない
        charge_map.resize(cm_extent[Environment::num_of_particle_types][static_cast<int>(num_cmat)]);
        temporary_charge_map.resize(cm_extent[Environment::num_of_particle_types][static_cast<int>(num_cmat)]);

        //! セルベースの物体定義 (物体内部判定用)
        for(const auto& cell_pos : cells) {
            object_cell_map[cell_pos[0]][cell_pos[1]][cell_pos[2]] = true;
        }

        //! 誘電体計算用の静電容量値をノードごとに計算
        if ( this->isDielectricSurface() ) {
            for (const auto& one_node : capacity_matrix_relation) {
                const auto& cmat_number = one_node.first;
                const auto& texture = textures.at(cmat_number);

                //! ノード周りのFaceに割り当てられたTextureの平均値から計算
                double capacitance = 0.0;

                for(const auto& texture_indicies : texture) {
                    if (material_capacitances.count( texture_indicies ) == 0) {
                        //! マッピングに使う物理定数の値を property_list から取り出しておく
                        if (material_property_list.count( material_names[ texture_indicies ] ) > 0 ) {
                            //! 各ノードの capacitance = 真空の誘電率 * 比誘電率でいい?
                            material_capacitances[ texture_indicies ] = Normalizer::normalizeCapacitance(
                                4.0 * M_PI * eps0 * material_property_list.at( material_names[ texture_indicies ] ).at("RelativePermittivity")
                            );

                            if (Environment::isRootNode) cout << format("  [INFO] Capacitance[%d] = %s F") % texture_indicies % Normalizer::unnormalizeCapacitance(material_capacitances[ texture_indicies ]) << endl;
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

//! Comm作成後に必要な初期化
void Spacecraft::initializeAfterMakeComm() {
    //! 粒子を放出するプロセス数を計算
    this->initializeEmissionParticleInfo();
}

void Spacecraft::initializeChargeMapOffset(const tdArray& phi) {
    if (initial_potential_offset != 0.0) {
        potential = initial_potential_offset;

        if (MPIw::Environment::isRootNode(name)) {
            cout << format("[INFO] [%s] offset_potential = %s V") % name % Normalizer::unnormalizePotential(potential) << endl;
        }

        for(unsigned int i = 0; i < num_cmat; ++i) {
            double delta_rho = 0.0;

            for(unsigned int j = 0; j < num_cmat; ++j) {
                if (isMyCmat(j)) {
                    const auto& pos = capacity_matrix_relation.at(j);
                    delta_rho += capacity_matrix(i, j) * potential;
                }
            }
            delta_rho = MPIw::Environment::Comms[name].sum(delta_rho);

            //! charge_mapは全プロセス同じ値にしてよい
            charge_map[0][i] += delta_rho;
        }
    }
}

void Spacecraft::saveWholeNodePositions(const ObjectNodes& whole_nodes) {
    //! 全体のノード位置を出力用、放出用に保持する
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
        throw std::invalid_argument("[ERROR] Invalid Cmat number passed to Spacecraft::getCmatPos().");
    }
}

inline unsigned int Spacecraft::getCmatNumber(const Position& pos) const {
    return this->getCmatNumber(pos.i, pos.j, pos.k);
}

unsigned int Spacecraft::getCmatNumber(const int i, const int j, const int k) const {
    //! Innerノードを探索
    for(const auto& node : capacity_matrix_relation) {
        const auto& pos = node.second;

        if (pos.i == i && pos.j == j && pos.k == k) return node.first;
    }

    //! Glueノードを探索
    for(const auto& node : glue_capacity_matrix_relation) {
        const auto& pos = node.second;

        if (pos.i == i && pos.j == j && pos.k == k) return node.first;
    }

    //! それでも見つからなければ全体ノード上を探す
    const auto abs_pos = Environment::getAbsolutePosition(i, j, k);
    for(const auto& node : whole_capacity_matrix_relation) {
        const auto& pos = node.second;

        if (pos.i == abs_pos[0] && pos.j == abs_pos[1] && pos.k == abs_pos[2]) return node.first;
    }

    std::string error_message = (format("[ERROR] %sInvalid Cmat Position (%d, %d, %d) passed to Spacecraft::getCmatNumber().") % Environment::rankStr() % i % j % k).str();
    throw std::invalid_argument(error_message);
}

bool Spacecraft::isCmatNode(const int i, const int j, const int k) const {
    try {
        Spacecraft::getCmatNumber(i, j, k);
    } catch (std::invalid_argument& e) {
        cout << Environment::rankStr() << format("It is not a cmat node: %d, %d, %d") % i % j % k << endl;
        return false;
    }

    return true;
}

auto Spacecraft::getTotalCharge() const {
    double q = 0.0;

    for(size_t cmat_itr = 0; cmat_itr < num_cmat; ++cmat_itr) {
        for(size_t pid = 0; pid < Environment::num_of_particle_types; ++pid) {
            //! charge_map は既に同期していると考えてよいので通信しなくてよい
            q += charge_map[pid][cmat_itr];
        }
    }

    return q;
}

double Spacecraft::getNodeCharge(const unsigned int cmat_itr) const {
    double node_charge = 0.0;
    for (int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
        node_charge += charge_map[pid][cmat_itr];
    }
    return node_charge;
}

double Spacecraft::getMaxPotential() const {
    double phi_max = std::numeric_limits<double>::lowest();
    const auto& phi = parent_grid->getPhi();

    for(const auto& node : capacity_matrix_relation) {
        const auto j = node.first;
        const auto& pos = node.second;
        const auto t_phi = phi[pos.i][pos.j][pos.k];
        phi_max = std::max(t_phi, phi_max);
    }

    return phi_max;
}

double Spacecraft::getMinPotential() const {
    double phi_min = std::numeric_limits<double>::max();
    const auto& phi = parent_grid->getPhi();

    for(const auto& node : capacity_matrix_relation) {
        const auto j = node.first;
        const auto& pos = node.second;
        const auto t_phi = phi[pos.i][pos.j][pos.k];
        phi_min = std::min(t_phi, phi_min);
    }

    return phi_min;
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
    //! C_ij の sum を計算して保存しておく
    total_cmat_value = capacity_matrix.total();
}

void Spacecraft::makeCmatrixInvert(void) {
    //! B行列 -> C行列に変換
    capacity_matrix = capacity_matrix.inverse();
    this->updateTotalCmatValue();
}

inline bool Spacecraft::isContaining(const Particle& p) const {
    return this->isContaining(p.getPosition());
}

inline bool Spacecraft::isContaining(const Position& pos) const {
    if (!is_defined_in_this_process) return false;

    if (Environment::isValidCellPosition(pos)) {
        return object_cell_map[pos.i][pos.j][pos.k];
    } else {
        return false;
    }
}

inline bool Spacecraft::isXsurfaceMinus(const int i, const int j, const int k) const {
    return (object_cell_map[i][j][k]) && (!object_cell_map[i - 1][j][k]);
}

inline bool Spacecraft::isXsurfacePlus(const int i, const int j, const int k) const {
    return (!object_cell_map[i][j][k]) && (object_cell_map[i - 1][j][k]);
}

inline bool Spacecraft::isYsurfaceMinus(const int i, const int j, const int k) const {
    return (object_cell_map[i][j][k]) && (!object_cell_map[i][j - 1][k]);
}

inline bool Spacecraft::isYsurfacePlus(const int i, const int j, const int k) const {
    return (!object_cell_map[i][j][k]) && (object_cell_map[i][j - 1][k]);
}

inline bool Spacecraft::isZsurfaceMinus(const int i, const int j, const int k) const {
    return (object_cell_map[i][j][k]) && (!object_cell_map[i][j][k - 1]);
}

inline bool Spacecraft::isZsurfacePlus(const int i, const int j, const int k) const {
    return (!object_cell_map[i][j][k]) && (object_cell_map[i][j][k - 1]);
}

//! Cmat NodeがXSurfaceかどうかを判定する
//! @note:Inner Cmat Nodeに対して呼ぶ == CellPositionはValidである
inline bool Spacecraft::isXsurfaceCmatNode(const Position& pos, const int sign) const {
    if (sign > 0) {
        return (
            isXsurfaceMinus(pos.i, pos.j, pos.k) ||
            isXsurfaceMinus(pos.i, pos.j - 1, pos.k) ||
            isXsurfaceMinus(pos.i, pos.j, pos.k - 1) ||
            isXsurfaceMinus(pos.i, pos.j - 1, pos.k - 1)
        );
    } else {
        return (
            isXsurfacePlus(pos.i, pos.j, pos.k) ||
            isXsurfacePlus(pos.i, pos.j - 1, pos.k) ||
            isXsurfacePlus(pos.i, pos.j, pos.k - 1) ||
            isXsurfacePlus(pos.i, pos.j - 1, pos.k - 1)
        );
    }
    return false;
}

inline bool Spacecraft::isXsurfaceCmatNode(const Position& pos) const {
    return isXsurfaceCmatNode(pos, 1) || isXsurfaceCmatNode(pos, -1);
}

//! 中間点がXSurface上の点かどうかを判定する
inline bool Spacecraft::isXsurfacePoint(const Position& pos, const int sign) const {
    constexpr double possible_error = 1e-10;
    if (pos.dx1 >= possible_error) return false;

    if (Environment::isValidCellPosition(pos) && Environment::isValidCellPosition(pos.i - 1, pos.j, pos.k)) {
        if (sign > 0) {
            return isXsurfaceMinus(pos.i, pos.j, pos.k);
        } else {
            return isXsurfacePlus(pos.i, pos.j, pos.k);
        }
    } else {
        //! Cell Positionで判定できない場合はCmatが存在するかどうかで判定する
        return isCmatNode(pos.i, pos.j, pos.k) && isCmatNode(pos.i, pos.j + 1, pos.k) && isCmatNode(pos.i, pos.j, pos.k + 1) && isCmatNode(pos.i, pos.j + 1, pos.k + 1);
    }
}

inline bool Spacecraft::isYsurfacePoint(const Position& pos, const int sign) const {
    constexpr double possible_error = 1e-10;
    if (pos.dy1 >= possible_error) return false;

    if (Environment::isValidCellPosition(pos) && Environment::isValidCellPosition(pos.i, pos.j - 1, pos.k)) {
        if (sign > 0) {
            return isYsurfaceMinus(pos.i, pos.j, pos.k);
        } else {
            return isYsurfacePlus(pos.i, pos.j, pos.k);
        }
    } else {
        return isCmatNode(pos.i, pos.j, pos.k) && isCmatNode(pos.i + 1, pos.j, pos.k) && isCmatNode(pos.i, pos.j, pos.k + 1) && isCmatNode(pos.i + 1, pos.j, pos.k + 1);
    }
}

inline bool Spacecraft::isZsurfacePoint(const Position& pos, const int sign) const {
    constexpr double possible_error = 1e-10;
    if (pos.dz1 >= possible_error) return false;

    if (Environment::isValidCellPosition(pos) && Environment::isValidCellPosition(pos.i, pos.j, pos.k - 1)) {
        if (sign > 0) {
            return isZsurfaceMinus(pos.i, pos.j, pos.k);
        } else {
            return isZsurfacePlus(pos.i, pos.j, pos.k);
        }
    } else {
        return isCmatNode(pos.i, pos.j, pos.k) && isCmatNode(pos.i + 1, pos.j, pos.k) && isCmatNode(pos.i, pos.j + 1, pos.k) && isCmatNode(pos.i + 1, pos.j + 1, pos.k);
    }
}

void Spacecraft::removeInnerParticle(Particle& p) const {
    if (isContaining(p)) { p.makeInvalid(); }
}

void Spacecraft::distributeInnerParticleCharge(Particle& p) {
    if (isContaining(p)) {
        const auto id = p.typeId;
        const auto q = p.getChargeOfSuperParticle();

        const int isign = (p.vx > 0.0) ? 1 : -1;
        const int jsign = (p.vy > 0.0) ? 1 : -1;
        const int ksign = (p.vz > 0.0) ? 1 : -1;
        bool surface_is_found = false;

        if (this->hasSecondaryParticles()) {
            //! 二次電子がある場合はIncidentの情報保存を行う必要がある
            auto cross_points = p.computeCrossPoints();
            for(auto& cross_point : cross_points) {
                if (isXsurfacePoint(cross_point, isign)) {
                    distributeInnerParticleChargeToXsurface(cross_point, id, q);
                    //! move_x_ratioが前のグリッドからの移動割合
                    this->addIncidentEvent(p.getParticleTypePtr(), p.getXMoveRatio(), cross_point, p.getVelocity(), AXIS::x);
                    surface_is_found = true;
                    break;
                } else if (isYsurfacePoint(cross_point, jsign)) {
                    distributeInnerParticleChargeToYsurface(cross_point, id, q);
                    this->addIncidentEvent(p.getParticleTypePtr(), p.getYMoveRatio(), cross_point, p.getVelocity(), AXIS::y);
                    surface_is_found = true;
                    break;
                } else if (isZsurfacePoint(cross_point, ksign)) {
                    distributeInnerParticleChargeToZsurface(cross_point, id, q);
                    this->addIncidentEvent(p.getParticleTypePtr(), p.getZMoveRatio(), cross_point, p.getVelocity(), AXIS::z);
                    surface_is_found = true;
                    break;
                }
            }

            if (!surface_is_found) {
                cout << "[ERROR] " << Environment::rankStr() <<
                    "Surface cannot detect on Spacecraft::distributeInnerParticleCharge():\n";
                cout << p << '\n';
                cout << "Cross Points:\n";
                for(auto& cross_point : cross_points) {
                    cout << cross_point << "\n";
                }
                cout << endl;
            }
        } else {
            auto cross_points = p.computeCrossPoints();
            for(auto& cross_point : cross_points) {
                if (isXsurfacePoint(cross_point, isign)) {
                    distributeInnerParticleChargeToXsurface(cross_point, id, q);
                    surface_is_found = true;
                    break;
                } else if (isYsurfacePoint(cross_point, jsign)) {
                    distributeInnerParticleChargeToYsurface(cross_point, id, q);
                    surface_is_found = true;
                    break;
                } else if (isZsurfacePoint(cross_point, ksign)) {
                    distributeInnerParticleChargeToZsurface(cross_point, id, q);
                    surface_is_found = true;
                    break;
                }
            }

            if (!surface_is_found) {
                cout << "[ERROR] " << Environment::rankStr() <<
                    "Surface cannot detect on Spacecraft::distributeInnerParticleCharge():\n";
                cout << p << '\n';
                cout << "Cross Points:\n";
                for(auto& cross_point : cross_points) {
                    cout << cross_point << "\n";
                }
                cout << endl;
            }
        }

        current[id] += q;
        p.makeInvalid();
    }
}

inline void Spacecraft::distributeInnerParticleChargeToXsurface(const Position& pos, const int id, const double charge) {
    temporary_charge_map[id][this->getCmatNumber(pos.i, pos.j    , pos.k    )] += charge * pos.dy2 * pos.dz2;
    temporary_charge_map[id][this->getCmatNumber(pos.i, pos.j + 1, pos.k    )] += charge * pos.dy1 * pos.dz2;
    temporary_charge_map[id][this->getCmatNumber(pos.i, pos.j    , pos.k + 1)] += charge * pos.dy2 * pos.dz1;
    temporary_charge_map[id][this->getCmatNumber(pos.i, pos.j + 1, pos.k + 1)] += charge * pos.dy1 * pos.dz1;
}

inline void Spacecraft::distributeInnerParticleChargeToYsurface(const Position& pos, const int id, const double charge) {
    temporary_charge_map[id][this->getCmatNumber(pos.i    , pos.j, pos.k    )] += charge * pos.dx2 * pos.dz2;
    temporary_charge_map[id][this->getCmatNumber(pos.i + 1, pos.j, pos.k    )] += charge * pos.dx1 * pos.dz2;
    temporary_charge_map[id][this->getCmatNumber(pos.i    , pos.j, pos.k + 1)] += charge * pos.dx2 * pos.dz1;
    temporary_charge_map[id][this->getCmatNumber(pos.i + 1, pos.j, pos.k + 1)] += charge * pos.dx1 * pos.dz1;
}

inline void Spacecraft::distributeInnerParticleChargeToZsurface(const Position& pos, const int id, const double charge) {
    temporary_charge_map[id][this->getCmatNumber(pos.i    , pos.j    , pos.k)] += charge * pos.dx2 * pos.dy2;
    temporary_charge_map[id][this->getCmatNumber(pos.i + 1, pos.j    , pos.k)] += charge * pos.dx1 * pos.dy2;
    temporary_charge_map[id][this->getCmatNumber(pos.i    , pos.j + 1, pos.k)] += charge * pos.dx2 * pos.dy1;
    temporary_charge_map[id][this->getCmatNumber(pos.i + 1, pos.j + 1, pos.k)] += charge * pos.dx1 * pos.dy1;
}

void Spacecraft::sumWholeCharge() {
    //! 多分ムーブ演算なのでそんなに高コストではないはず
    temporary_charge_map = MPIw::Environment::Comms[name].sum(temporary_charge_map);

    #pragma omp parallel for
    for (unsigned int cmat_itr = 0; cmat_itr < num_cmat; ++cmat_itr) {
        for (unsigned int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
            charge_map[pid][cmat_itr] += temporary_charge_map[pid][cmat_itr];
            temporary_charge_map[pid][cmat_itr] = 0.0;
        }
    }
}

void Spacecraft::redistributeCharge(RhoArray& rho, const tdArray& phi) {
    total_charge = getTotalCharge();

    if (this->isDielectricSurface()) {
        this->redistributeChargeForDielectric(rho, phi);
    } else {
        this->redistributeChargeForPerfectConductor(rho, phi);
    }

    // update total current for output
    sumCurrent();
}

void Spacecraft::redistributeChargeForDielectric(RhoArray& rho, const tdArray& phi) {
    double capacity_times_phi = 0.0;

    //! 誘電体の場合
    for(const auto& node : capacity_matrix_relation) {
        const auto j = node.first;
        const auto& pos = node.second;
        const auto capacitance = capacitance_map[j];

        if (capacitance > 0.0) {
            double node_charge = 0.0;
            for (int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
                node_charge += charge_map[pid][j];
            }
            for(size_t i = 0; i < num_cmat; ++i) {
                capacity_times_phi += capacity_matrix(i, j) * (phi[pos.i][pos.j][pos.k] - node_charge / capacitance);
            }
        } else {
            for(size_t i = 0; i < num_cmat; ++i) {
                capacity_times_phi += capacity_matrix(i, j) * phi[pos.i][pos.j][pos.k];
            }
        }
    }

    capacity_times_phi = MPIw::Environment::Comms[name].sum(capacity_times_phi);

    if (is_potential_fixed) {
        potential = fixed_potential;
    } else {
        potential = (capacity_times_phi + total_charge) / total_cmat_value;
    }

    if (MPIw::Environment::isRootNode(name)) {
        cout << format("[INFO] [%s] potential = %s V") % name % Normalizer::unnormalizePotential(potential) << endl;
    }

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
                        node_charge += charge_map[pid][j];
                    }
                    delta_rho += capacity_matrix(i, j) * (potential - phi[pos.i][pos.j][pos.k] + (node_charge / capacitance));
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
}

void Spacecraft::redistributeChargeForPerfectConductor(RhoArray& rho, const tdArray& phi) {
    if (MPIw::Environment::isRootNode(name)) {
        cout << format("[INFO] [%s] charge before redist: %16.7e") % name % total_charge << endl;
    }

    double capacity_times_phi = 0.0;
    //! 完全導体の場合
    for(const auto& one_node : capacity_matrix_relation) {
        const auto j = one_node.first;
        const auto& pos = one_node.second;
        for(size_t i = 0; i < num_cmat; ++i) {
            capacity_times_phi += capacity_matrix(i, j) * phi[pos.i][pos.j][pos.k];
        }
    }

    capacity_times_phi = MPIw::Environment::Comms[name].sum(capacity_times_phi);

    if (is_potential_fixed) {
        potential = fixed_potential;
    } else {
        potential = (capacity_times_phi + total_charge) / total_cmat_value;
    }

    if (MPIw::Environment::isRootNode(name)) {
        cout << format("[INFO] [%s] potential = %s V") % name % Normalizer::unnormalizePotential(potential) << endl;
    }

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

void Spacecraft::initializeEmissionParticleInfo() {
    for(auto& pinfo : emit_particle_info) {
        auto& info = pinfo.second;

        //! 粒子放出に参加するプロセス数を計算する
        double is_valid_node_position = 0.0;
        if (Environment::isValidNodePosition(info.relative_emission_position)) is_valid_node_position = 1.0;
        info.emission_process_number = MPIw::Environment::Comms[name].sum(is_valid_node_position);
    }
}

//! 粒子放出
void Spacecraft::emitParticles(ParticleArray& parray, const double dx) {
    for(const auto& pinfo : emit_particle_info) {
        const auto id = pinfo.first;
        const auto& info = pinfo.second;
        const auto emit_ptype_ptr = Environment::getEmissionParticleType(id);

        if (emit_ptype_ptr->getType() == "beam") {
            //! ビーム粒子
            if (Environment::isValidNodePosition(info.relative_emission_position)) {
                //! 放出時に呼び出すメンバ関数ポインタ
                std::function<void(Position&, const int, const double)> emission_func;
                std::function<Position(Particle&)> get_next_position_func;

                //! 放出時の座標修正子
                if (fabs(info.emission_vector[0]) == 1.0) {

                    emission_func = [this](Position& pos, const int id, const double charge) {
                        this->subtractChargeOfParticleFromXsurface(pos, id, charge);
                    };
                    get_next_position_func = [](const Particle& p) { return p.getNextXCrossPoint(); };

                } else if (fabs(info.emission_vector[1]) == 1.0) {

                    emission_func = [this](Position& pos, const int id, const double charge) {
                        this->subtractChargeOfParticleFromYsurface(pos, id, charge);
                    };
                    get_next_position_func = [](const Particle& p) { return p.getNextYCrossPoint(); };

                } else if (fabs(info.emission_vector[2]) == 1.0) {

                    emission_func = [this](Position& pos, const int id, const double charge) {
                        this->subtractChargeOfParticleFromZsurface(pos, id, charge);
                    };
                    get_next_position_func = [](const Particle& p) { return p.getNextZCrossPoint(); };

                } else {
                    std::string error_message = (format("[ERROR] At Spacecraft::emitParticles: Now 'emission_vector' must be [+-1,0,0] or [0,+-1,0] or [0,0,+-1].")).str();
                    throw std::invalid_argument(error_message);
                }

                const auto beam_ptype_ptr = Environment::getBeamParticleType(id);
                const auto max_amount = beam_ptype_ptr->getEmissionAmount() / info.emission_process_number;
                const auto charge = beam_ptype_ptr->getChargeOfSuperParticle();

                for(int i = 0; i < max_amount; ++i) {
                    Particle p = beam_ptype_ptr->generateNewParticle(info.relative_emission_position, info.emission_vector);

                    auto pos = get_next_position_func(p);
                    emission_func(pos, id, charge);
                    parray[id].push_back( std::move(p) );
                }
            }
        } else if (emit_ptype_ptr->getType() == "photoelectron") {
            //! 光電子
            std::function<void(Position&, const int, const double)> emission_func;
            std::function<Position(Particle&)> get_next_position_func;
            std::function<bool(const Position&)> is_surface_func;
            const auto shine_vector = Environment::getStaticField().getShineVector();

            //! 放出時の座標修正子
            if (fabs(shine_vector[0]) == 1.0) {

                emission_func = [this](Position& pos, const int id, const double charge) {
                    this->subtractChargeOfParticleFromXsurface(pos, id, charge);
                };
                get_next_position_func = [](const Particle& p) { return p.getNextXCrossPoint(); };

                const int surface_sign = (shine_vector[0] > 0.0) ? 1 : -1;
                is_surface_func = [this, sign = surface_sign](const Position& pos) { return this->isXsurfacePoint(pos, sign); };

            } else if (fabs(shine_vector[1]) == 1.0) {

                emission_func = [this](Position& pos, const int id, const double charge) {
                    this->subtractChargeOfParticleFromYsurface(pos, id, charge);
                };
                get_next_position_func = [](const Particle& p) { return p.getNextYCrossPoint(); };

                const int surface_sign = (shine_vector[1] > 0.0) ? 1 : -1;
                is_surface_func = [this, sign = surface_sign](const Position& pos) { return this->isYsurfacePoint(pos, sign); };

            } else if (fabs(shine_vector[2]) == 1.0) {

                emission_func = [this](Position& pos, const int id, const double charge) {
                    this->subtractChargeOfParticleFromZsurface(pos, id, charge);
                };
                get_next_position_func = [](const Particle& p) { return p.getNextZCrossPoint(); };

                const int surface_sign = (shine_vector[2] > 0.0) ? 1 : -1;
                is_surface_func = [this, sign = surface_sign](const Position& pos) { return this->isZsurfacePoint(pos, sign); };

            } else {

                std::string error_message = (format("[ERROR] At Spacecraft::emitParticles: Now 'shine_vector' must be [+-1,0,0] or [0,+-1,0] or [0,0,+-1].")).str();
                throw std::invalid_argument(error_message);

            }

            const auto pe_ptype_ptr = Environment::getPhotoElectronParticleType(id);
            const auto charge = pe_ptype_ptr->getChargeOfSuperParticle();
            const auto area = dx * dx;
            const auto emission_amount_per_area = area * pe_ptype_ptr->getEmissionAmount();

            //! emission_vectorはshine_vectorの逆
            const std::array<double, 3> emission_vector{ -shine_vector[0], -shine_vector[1], -shine_vector[2] };

            for(const auto& node : capacity_matrix_relation) {
                const auto j = node.first;
                const auto& pos = node.second;

                if (is_surface_func(pos)) {
                    for(int i = 0; i < emission_amount_per_area; ++i) {
                        Particle p = pe_ptype_ptr->generateNewParticle(pos, emission_vector);

                        auto emiss_pos = get_next_position_func(p);
                        emission_func(emiss_pos, id, charge);
                        parray[id].push_back( std::move(p) );
                    }
                }
            }
        } else if (emit_ptype_ptr->getType() == "secondary") {
            const auto secondary_ptype_ptr = Environment::getSecondaryParticleType(id);
            const auto charge = secondary_ptype_ptr->getChargeOfSuperParticle();

            MaterialInfo_t test_material{"test"};
            test_material.fermi_energy = 9.1;
            test_material.delta_max = 3.0;
            test_material.epsi_max = 420.0;
            test_material.atomic_number = 13.0;

            for(auto& incident : incident_events) {
                auto generated_parray = secondary_ptype_ptr->generateNewParticles(incident, test_material);
                for(auto& new_particle : generated_parray) {
                    if (isContaining(new_particle)) {
                        //! 放出後再度内部に入ってしまった場合を取り扱う
                        this->distributeInnerParticleChargeForSecondary(new_particle, incident.getAxis());
                    } else {
                        if (incident.isXsurfaceIncident()) {
                            auto cross_pos = new_particle.getOldPositionParticle().getNextXCrossPoint();
                            this->subtractChargeOfParticleFromXsurface(cross_pos, id, charge);
                        } else if (incident.isYsurfaceIncident()) {
                            auto cross_pos = new_particle.getOldPositionParticle().getNextYCrossPoint();
                            this->subtractChargeOfParticleFromYsurface(cross_pos, id, charge);
                        } else {
                            auto cross_pos = new_particle.getOldPositionParticle().getNextZCrossPoint();
                            this->subtractChargeOfParticleFromZsurface(cross_pos, id, charge);
                        }
                        parray[id].push_back(new_particle);
                    }
                }
            }
        }
    }

    if (this->hasSecondaryParticles()) {
        this->clearIncidentEvents();
    }
}

inline void Spacecraft::subtractChargeOfParticleFromXsurface(const Position& pos, const int id, const double charge) {
    temporary_charge_map[id][this->getCmatNumber(pos.i, pos.j    , pos.k    )] -= charge * pos.dy2 * pos.dz2;
    temporary_charge_map[id][this->getCmatNumber(pos.i, pos.j + 1, pos.k    )] -= charge * pos.dy1 * pos.dz2;
    temporary_charge_map[id][this->getCmatNumber(pos.i, pos.j    , pos.k + 1)] -= charge * pos.dy2 * pos.dz1;
    temporary_charge_map[id][this->getCmatNumber(pos.i, pos.j + 1, pos.k + 1)] -= charge * pos.dy1 * pos.dz1;

    current[id] -= charge;
}

inline void Spacecraft::subtractChargeOfParticleFromYsurface(const Position& pos, const int id, const double charge) {
    temporary_charge_map[id][this->getCmatNumber(pos.i    , pos.j, pos.k    )]-= charge * pos.dx2 * pos.dz2;
    temporary_charge_map[id][this->getCmatNumber(pos.i + 1, pos.j, pos.k    )]-= charge * pos.dx1 * pos.dz2;
    temporary_charge_map[id][this->getCmatNumber(pos.i    , pos.j, pos.k + 1)]-= charge * pos.dx2 * pos.dz1;
    temporary_charge_map[id][this->getCmatNumber(pos.i + 1, pos.j, pos.k + 1)]-= charge * pos.dx1 * pos.dz1;

    current[id] -= charge;
}

inline void Spacecraft::subtractChargeOfParticleFromZsurface(const Position& pos, const int id, const double charge) {
    temporary_charge_map[id][this->getCmatNumber(pos.i    , pos.j    , pos.k)] -= charge * pos.dx2 * pos.dy2;
    temporary_charge_map[id][this->getCmatNumber(pos.i + 1, pos.j    , pos.k)] -= charge * pos.dx1 * pos.dy2;
    temporary_charge_map[id][this->getCmatNumber(pos.i    , pos.j + 1, pos.k)] -= charge * pos.dx2 * pos.dy1;
    temporary_charge_map[id][this->getCmatNumber(pos.i + 1, pos.j + 1, pos.k)] -= charge * pos.dx1 * pos.dy1;

    current[id] -= charge;
}

inline bool Spacecraft::isValidEmission(Particle& p) const {
    //! 放出前は内部、放出後は外部に入れば valid と見なす
    return ( (isContaining(p)) && !isContaining(p.getNextPosition()) );
}

//! I/O関連
std::string Spacecraft::getLogHeader() const {
    std::string format_string =
        (this->isDielectricSurface()) ? "%16s %16s %16s" : "%16s %16s";
    for(int i = 0; i < Environment::num_of_particle_types; ++i) {
        format_string += " %16s";
    }

    if (this->forceComputationEnabled()) {
        format_string += " %16.7e %16.7e %16.7e";
    }

    auto format_base = format(std::move(format_string));

    if (this->isDielectricSurface()) {
        format_base = format_base % "PhiMax [V]";
        format_base = format_base % "PhiMin [V]";
    } else {
        format_base = format_base % "Potential [V]";
    }

    format_base = format_base % "Charge [C]";
    for(int i = 0; i < Environment::num_of_particle_types; ++i) {
        format_base = format_base % (Environment::getParticleType(i)->getName() + " [A]");
    }

    if (this->forceComputationEnabled()) {
        format_base = format_base % "Fx [N]" % "Fy [N]" % "Fz [N]";
    }

    std::string header = format_base.str();
    return header;
}

std::string Spacecraft::getLogEntry() const {
    std::string format_string =
        (this->isDielectricSurface()) ? "%16.7e %16.7e %16.7e" : "%16.7e %16.7e";

    for(int i = 0; i < Environment::num_of_particle_types; ++i) {
        format_string += " %16.7e";
    }

    if (this->forceComputationEnabled()) {
        format_string += " %16.7e %16.7e %16.7e";
    }

    auto format_base = format(std::move(format_string));

    if (this->isDielectricSurface()) {
        format_base = format_base % Normalizer::unnormalizePotential(this->getMaxPotential());
        format_base = format_base % Normalizer::unnormalizePotential(this->getMinPotential());
    } else {
        format_base = format_base % Normalizer::unnormalizePotential(potential);
    }

    format_base = format_base % Normalizer::unnormalizeCharge(total_charge);
    for(int i = 0; i < Environment::num_of_particle_types; ++i) {
        format_base = format_base % Normalizer::unnormalizeCurrent(current[i]);
    }

    if (this->forceComputationEnabled()) {
        auto norm = Normalizer::unnormalizeForce(1.0);
        format_base = format_base % (norm * force.fx) % (norm * force.fy) % (norm * force.fz);
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
