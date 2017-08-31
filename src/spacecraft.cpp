#include "spacecraft.hpp"
#include "particle.hpp"
#include "normalizer.hpp"
#include "utils.hpp"
#include <fstream>
#include <stdexcept>
#include <regex>

//! static 変数の実体
unsigned int Spacecraft::num_of_spacecraft = 0;

void Spacecraft::construct(const size_t nx, const size_t ny, const size_t nz, const ObjectInfo_t& obj_info, const ObjectNodes& nodes, const ObjectNodes& glue_nodes, const ObjectCells& cells) {
    //! このオブジェクトがプロセス内で有効かどうかを保存しておく
    is_defined_in_this_process = (nodes.size() > 0);
    ++num_of_spacecraft;
    potential = 0.0;
    total_charge = 0.0;
    potential_fix = Normalizer::normalizePotential(obj_info.potential_fix);

    if (is_defined_in_this_process) {
        // Node ベース, Glueセルも必要
        ObjectDefinedMapBool::extent_gen objectBoolExtents;
        object_node_map.resize(objectBoolExtents[nx + 2][ny + 2][nz + 2]);

        //! セルマップは int で texture_index を持つ
        ObjectDefinedMapInt::extent_gen objectIntExtents;
        object_cell_map.resize(objectIntExtents[nx + 1][ny + 1][nz + 1]);

        // 物体定義マップを初期化
        for(int i = 0; i < nx + 2; ++i) {
            for (int j = 0; j < ny + 2; ++j) {
                for (int k = 0; k < nz + 2; ++k) {
                    object_node_map[i][j][k] = false;

                    if (i != nx + 1 && j != ny + 1 && k != nz + 1) object_cell_map[i][j][k] = 0;
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

            object_node_map[i][j][k] = true;
            capacity_matrix_relation.emplace(std::piecewise_construct, std::make_tuple(cmat_itr), std::make_tuple(i, j, k));
        }

        //! Glueノードはobject_node_map側を更新するだけでよい
        for(const auto& node_pair : glue_nodes) {
            const auto& node_pos = node_pair.second;
            const auto i = node_pos[0];
            const auto j = node_pos[1];
            const auto k = node_pos[2];

            object_node_map[i][j][k] = true;
        }

        for(const auto& cell_pos : cells) {
            object_cell_map[cell_pos[0]][cell_pos[1]][cell_pos[2]] = cell_pos[3];
            cout << Environment::rankStr() << format("cell[%d][%d][%d] is texture: %d") % cell_pos[0] % cell_pos[1] % cell_pos[2] % cell_pos[3] << endl;
        }

        //! キャパシタンス行列のサイズを物体サイズに変更
        capacity_matrix.resize(num_cmat, num_cmat);
        
        tdArray::extent_gen tdExtents;
        for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
            //! Node ベース, glue cell ありの電荷密度マップを生成
            charge_map.emplace_back(tdExtents[nx + 2][ny + 2][nz + 2], boost::fortran_storage_order());

            //! 電流 in/out を記録する要素を初期化
            current.push_back(0.0);

            //! 放出粒子のID追加判定
            auto itr = std::find(obj_info.emit_particle_names.begin(), obj_info.emit_particle_names.end(), Environment::getParticleType(pid)->getName());
            if (itr != obj_info.emit_particle_names.end()) {
                emit_particle_ids.push_back(pid);
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
    Utils::makeInvert(capacity_matrix);
    this->updateTotalCmatValue();
}

bool Spacecraft::isContaining(const Particle& p) const {
    if (!is_defined_in_this_process) return false;

    return this->isContaining(p.getPosition());
}

bool Spacecraft::isContaining(const Position& pos) const {
    if (!is_defined_in_this_process) return false;

    if (Environment::isValidPosition(pos)) {
        return (object_cell_map[pos.i][pos.j][pos.k] > 0);
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
        const auto pos = p.getPosition();

        charge_map[id][pos.i    ][pos.j    ][pos.k    ] += q * pos.dx2 * pos.dy2 * pos.dz2;
        charge_map[id][pos.i + 1][pos.j    ][pos.k    ] += q * pos.dx1 * pos.dy2 * pos.dz2;
        charge_map[id][pos.i    ][pos.j + 1][pos.k    ] += q * pos.dx2 * pos.dy1 * pos.dz2;
        charge_map[id][pos.i + 1][pos.j + 1][pos.k    ] += q * pos.dx1 * pos.dy1 * pos.dz2;
        charge_map[id][pos.i    ][pos.j    ][pos.k + 1] += q * pos.dx2 * pos.dy2 * pos.dz1;
        charge_map[id][pos.i + 1][pos.j    ][pos.k + 1] += q * pos.dx1 * pos.dy2 * pos.dz1;
        charge_map[id][pos.i    ][pos.j + 1][pos.k + 1] += q * pos.dx2 * pos.dy1 * pos.dz1;
        charge_map[id][pos.i + 1][pos.j + 1][pos.k + 1] += q * pos.dx1 * pos.dy1 * pos.dz1;

        current[id] += q;

        p.makeInvalid();
    }
}

void Spacecraft::subtractChargeOfParticle(const Particle& p) {
    const auto id = p.typeId;
    const auto q = p.getChargeOfSuperParticle();
    const auto pos = p.getOldPosition();

    charge_map[id][pos.i    ][pos.j    ][pos.k    ] -= q * pos.dx2 * pos.dy2 * pos.dz2;
    charge_map[id][pos.i + 1][pos.j    ][pos.k    ] -= q * pos.dx1 * pos.dy2 * pos.dz2;
    charge_map[id][pos.i    ][pos.j + 1][pos.k    ] -= q * pos.dx2 * pos.dy1 * pos.dz2;
    charge_map[id][pos.i + 1][pos.j + 1][pos.k    ] -= q * pos.dx1 * pos.dy1 * pos.dz2;
    charge_map[id][pos.i    ][pos.j    ][pos.k + 1] -= q * pos.dx2 * pos.dy2 * pos.dz1;
    charge_map[id][pos.i + 1][pos.j    ][pos.k + 1] -= q * pos.dx1 * pos.dy2 * pos.dz1;
    charge_map[id][pos.i    ][pos.j + 1][pos.k + 1] -= q * pos.dx2 * pos.dy1 * pos.dz1;
    charge_map[id][pos.i + 1][pos.j + 1][pos.k + 1] -= q * pos.dx1 * pos.dy1 * pos.dz1;

    current[id] -= q;
}

void Spacecraft::applyCharge(RhoArray& rho) const {
    //! 電荷分布を場に印加する
    for(const auto& one_node : capacity_matrix_relation) {
        const auto& pos = one_node.second;

        for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
            rho[pid + 1][pos.i][pos.j][pos.k] += charge_map[pid][pos.i][pos.j][pos.k];
        }
    }
}

void Spacecraft::redistributeCharge(RhoArray& rho, const tdArray& phi) {
    auto q = getTotalCharge(rho);

    if (MPIw::Environment::isRootNode(name)) {
        cout << format("%s: %16.7e") % "charge before redist" % q << endl;
    }

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

    if (MPIw::Environment::isRootNode(name)) {
        cout << "[" << name << "] potential = " << Normalizer::unnormalizePotential(potential) << " V. " << endl;
    }

    for(unsigned int i = 0; i < num_cmat; ++i) {
        double delta_rho = 0.0;

        for(unsigned int j = 0; j < num_cmat; ++j) {
            if (isMyCmat(j)) {
                const auto& pos = capacity_matrix_relation.at(j);
                delta_rho += capacity_matrix(i, j) * (potential - phi[pos.i][pos.j][pos.k]);
            }
        }
        delta_rho = MPIw::Environment::Comms[name].sum(delta_rho);

        if (isMyCmat(i)) {
            const auto& target_pos = capacity_matrix_relation.at(i);
            rho[0][target_pos.i][target_pos.j][target_pos.k] += delta_rho;
        }
    }

    q = getTotalCharge(rho);
    // update total charge for output
    total_charge = q;

    if (MPIw::Environment::isRootNode(name)) {
        cout << format("%s: %16.7e") % "charge after redist " % q << endl;
    }
}

//! 粒子放出関連
void Spacecraft::emitParticles(ParticleArray& parray) {
    for(const auto& id : emit_particle_ids) {
        const auto emit_ptype_ptr = Environment::getEmissionParticleType(id);
        const auto max_amount = emit_ptype_ptr->getEmissionAmount();

        for(int i = 0; i < max_amount; ++i) {
            Particle p = emit_ptype_ptr->generateNewParticle();

            if (this->isValidEmission(p)) {
                this->subtractChargeOfParticle(p);
                parray[id].push_back( std::move(p) );
            }
        }
    }
}

bool Spacecraft::isValidEmission(Particle& p) const {
    //! 放出前は内部、放出後は外部に入れば valid と見なす
    return ( (!isContaining(p)) && isContaining(p.getOldPosition()) );
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
    ost << "              name: " << spc.name << endl;
    ost << "              cmat: " << spc.num_cmat << endl;
    ost << "      surface_type: " << spc.surface_type << endl;

    if (spc.isDielectricSurface()) {
        for(const auto& mt : spc.materials) {
            ost << "            index " << mt.first << ": " << mt.second << endl;
        }
    }

    ost << "    emit particles: ";
    for(const auto& id : spc.emit_particle_ids) {
        ost << Environment::getParticleType(id)->getName() << ", " << endl;
    }

    return ost;
}


// Utility Functions for Objects
namespace ObjectUtils {
    ObjectDataFromFile getObjectNodesFromObjFile(const std::string& obj_file_name) {
        ObjectDataFromFile obj_data;

        if (Utils::isExistingFile(obj_file_name)) {
            const auto nx = Environment::nx;
            const auto ny = Environment::ny;
            const auto nz = Environment::nz;
            ObjectDefinedMapBool object_node_map(boost::extents[nx][ny][nz]);

            using ObjectFaceMap = boost::multi_array<int, 4>;
            ObjectFaceMap object_xface_map(boost::extents[nx][ny - 1][nz - 1][2]);
            ObjectFaceMap object_yface_map(boost::extents[nx - 1][ny][nz - 1][2]);
            ObjectFaceMap object_zface_map(boost::extents[nx - 1][ny - 1][nz][2]);

            //! 初期化
            for(int i = 0; i < nx; ++i) {
                for (int j = 0; j < ny; ++j) {
                    for (int k = 0; k < nz; ++k) {
                        object_node_map[i][j][k] = false;
                    }
                }
            }

            std::ifstream file_input(obj_file_name, std::ios::in);
            std::string buffer;

            //! vとf 始まりの文字列にマッチするパターン
            std::regex re_vertex(R"(^v (.*)$)");
            std::regex re_face(R"(^f (.*)$)");

            //! vertexとfaceを一時格納しておくarray
            using Vertex = std::array<int, 3>;
            using Face = std::array<int, 5>;
            std::vector<Vertex> vertices;
            std::vector<Face> faces;

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
                            if (!object_node_map[i][j][k]) {
                                obj_data.nodes[num_cmat] = {{i, j, k}};
                                ++num_cmat;
                                object_node_map[i][j][k] = true;
                            }
                            if (j != maxy && k != maxz) {
                                object_xface_map[i][j][k][0] = 1;
                                object_xface_map[i][j][k][1] = face[4];
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
                            if (!object_node_map[i][j][k]) {
                                obj_data.nodes[num_cmat] = {{i, j, k}};
                                ++num_cmat;
                                object_node_map[i][j][k] = true;
                            }
                            if (i != maxx && k != maxz) {
                                object_yface_map[i][j][k][0] = 1;
                                object_yface_map[i][j][k][1] = face[4];
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
                            if (!object_node_map[i][j][k]) {
                                obj_data.nodes[num_cmat] = {{i, j, k}};
                                ++num_cmat;
                                object_node_map[i][j][k] = true;
                            }
                            if (i != maxx && j != maxy) {
                                object_zface_map[i][j][k][0] = 1;
                                object_zface_map[i][j][k][1] = face[4];
                            }
                        }
                    }
                } else {
                    throw std::logic_error("Face type cannot be determinted by vertices position.");
                }
            }

            //! セルマップ定義のためにFaceを使ってカウントする
            using CellCountMap = boost::multi_array<int, 4>;
            using CellCountMapSlice = CellCountMap::array_view<1>::type;
            using CellCountMapRange = CellCountMap::index_range;
            CellCountMap object_cell_count_map(boost::extents[nx - 1][ny - 1][nz - 1][7]);
            for (int j = 0; j < ny - 1; ++j) {
                for (int k = 0; k < nz - 1; ++k) {
                    for(int lower_index = 0; lower_index < nx - 1; ++lower_index) {
                        //! 下端を見つける
                        if (object_xface_map[lower_index][j][k][0] == 1) {
                            //! 下端が見つかった場合は上端を探す
                            for(int upper_index = lower_index + 1; upper_index < nx - 1; ++upper_index) {
                                if (object_xface_map[upper_index][j][k][0] == 1) {
                                    //! 上端が見つかったら間を塗る
                                    for(int i = lower_index; i < upper_index; ++i) {
                                        object_cell_count_map[i][j][k][0] += 1;

                                        //! 下端の面の色
                                        object_cell_count_map[i][j][k][1] = object_xface_map[lower_index][j][k][1];
                                        //! 上端の面の色
                                        object_cell_count_map[i][j][k][2] = object_xface_map[upper_index][j][k][1];
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
                        if (object_yface_map[i][lower_index][k][0] == 1) {
                            //! 下端が見つかった場合は上端を探す
                            for(int upper_index = lower_index + 1; upper_index < ny - 1; ++upper_index) {
                                if (object_yface_map[i][upper_index][k][0] == 1) {
                                    //! 上端が見つかったら間を塗る
                                    for(int j = lower_index; j < upper_index; ++j) {
                                        object_cell_count_map[i][j][k][0] += 1;

                                        //! 下端の面の色
                                        object_cell_count_map[i][j][k][3] = object_yface_map[i][lower_index][k][1];
                                        //! 上端の面の色
                                        object_cell_count_map[i][j][k][4] = object_yface_map[i][upper_index][k][1];
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
                        if (object_zface_map[i][j][lower_index][0] == 1) {
                            //! 下端が見つかった場合は上端を探す
                            for(int upper_index = lower_index + 1; upper_index < nz - 1; ++upper_index) {
                                if (object_zface_map[i][j][upper_index][0] == 1) {
                                    //! 上端が見つかったら間を塗る
                                    for(int k = lower_index; k < upper_index; ++k) {
                                        object_cell_count_map[i][j][k][0] += 1;

                                        //! 下端の面の色
                                        object_cell_count_map[i][j][k][5] = object_zface_map[i][j][lower_index][1];
                                        //! 上端の面の色
                                        object_cell_count_map[i][j][k][6] = object_zface_map[i][j][upper_index][1];
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
                        if (object_cell_count_map[i][j][k][0] == 3) {
                            int texture_index = getTextureIndex<CellCountMapSlice>( object_cell_count_map[ boost::indices[i][j][k][ CellCountMapRange(1, 7) ] ] );
                            obj_data.cells.push_back({i, j, k, texture_index});
                        }
                    }
                }
            }
        } else {
            std::string error_message = (format("File %s does not exist.") % obj_file_name).str();
            throw std::invalid_argument(error_message);
        }

        return obj_data;
    }

    template<typename T>
    int getTextureIndex(const T&& counts) {
        //! 6次元配列であることを期待する
        assert(6, counts.size());

        if ((counts[0] == counts[1]) && (counts[2] == counts[3]) && (counts[4] == counts[5])) {
            //! すべての face pair の面テクスチャが一致している場合、
            //! 数の少ない要素を優先する
            std::map<int, int> index_counts;

            for (int idx = 0; idx < 6; ++idx) {
                index_counts[counts[idx]] += 1;
            }

            int minimum_index = -1;
            int minimum_count = 100000;
            for (auto& pair : index_counts) {
                if (pair.second < minimum_count) minimum_index = pair.first;
            }

            return minimum_index;
        } else {
            // 4つのペアから1つの値を決定するラムダ

            //! 不一致のペアがあるなら、それ以外を優先
            if (counts[0] != counts[1]) {
                return determineOneValueFromFourElems(counts[2], counts[3], counts[4], counts[5]);
            } else if (counts[2] != counts[3]) {
                return determineOneValueFromFourElems(counts[0], counts[1], counts[4], counts[5]);
            } else {
                return determineOneValueFromFourElems(counts[0], counts[1], counts[2], counts[3]);
            }
        }

        throw std::logic_error("something wrong on getTextureIndex()...");
    }

    template<typename T>
    T determineOneValueFromFourElems(T v1, T v2, T v3, T v4) {
        if (v1 == v2 && v2 == v3 && v3 == v4 ) {
            //! 全部一致(基本的にこのケースが多いはず)
            return v1;
        } else if (v1 == v2 && v3 == v4) {
            return v1; // どっちも別の値でペアになっている時は早いやつを優先
        } else if (v3 != v4) {
            return v1; // この場合、v1 == v2 または v1 != v2 だが、どちらの場合も v1 を返すので判定不要
        } else {
            return v3;
        }
    }
}
