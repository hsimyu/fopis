#include "spacecraft.hpp"
#include "particle.hpp"
#include "normalizer.hpp"
#include "utils.hpp"
#include <fstream>
#include <stdexcept>
#include <regex>
#include <boost/filesystem.hpp>

//! static 変数の実体
unsigned int Spacecraft::num_of_spacecraft = 0;

void Spacecraft::construct(const size_t nx, const size_t ny, const size_t nz, const ObjectInfo_t& obj_info, const ObjectNodes& nodes, const ObjectNodes& glue_nodes, const ObjectFaces& faces) {
    //! このオブジェクトがプロセス内で有効かどうかを保存しておく
    is_defined_in_this_process = (nodes.size() > 0);
    ++num_of_spacecraft;
    potential = 0.0;
    total_charge = 0.0;
    potential_fix = Normalizer::normalizePotential(obj_info.potential_fix);

    if (is_defined_in_this_process) {
        // Node ベース, Glueセルも必要
        ObjectDefinedMap::extent_gen objectExtents;
        object_node_map.resize(objectExtents[nx + 2][ny + 2][nz + 2]);

        object_xface_map.resize(objectExtents[nx + 2][ny + 1][nz + 1]);
        object_yface_map.resize(objectExtents[nx + 1][ny + 2][nz + 1]);
        object_zface_map.resize(objectExtents[nx + 1][ny + 1][nz + 2]);

        // 物体定義マップを初期化
        for(int i = 0; i < nx + 2; ++i) {
            for (int j = 0; j < ny + 2; ++j) {
                for (int k = 0; k < nz + 2; ++k) {
                    object_node_map[i][j][k] = false;

                    if (j != ny + 1 && k != nz + 1) object_xface_map[i][j][k] = false;
                    if (i != nx + 1 && k != nz + 1) object_yface_map[i][j][k] = false;
                    if (i != nx + 1 && j != ny + 1) object_zface_map[i][j][k] = false;
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

        for(const auto& v : faces) {
            const auto face_type = v[0];
            switch(face_type) {
                case 0:
                    object_xface_map[ v[1] ][ v[2] ][ v[3] ] = true;
                    break;
                case 1:
                    object_yface_map[ v[1] ][ v[2] ][ v[3] ] = true;
                    break;
                case 2:
                    object_zface_map[ v[1] ][ v[2] ][ v[3] ] = true;
                    break;
                default:
                    break;
            }
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
    const auto i = pos.i;
    const auto j = pos.j;
    const auto k = pos.k;

    return  object_node_map[i    ][j    ][k    ] &&
            object_node_map[i + 1][j    ][k    ] &&
            object_node_map[i    ][j + 1][k    ] &&
            object_node_map[i + 1][j + 1][k    ] &&
            object_node_map[i    ][j    ][k + 1] &&
            object_node_map[i + 1][j    ][k + 1] &&
            object_node_map[i    ][j + 1][k + 1] &&
            object_node_map[i + 1][j + 1][k + 1];
}

/*
bool Spacecraft::isIncludedByFace(const Particle& p) const {
    if (!is_defined_in_this_process) return false;

    const auto pos = p.getPosition();
    const auto i = pos.i;
    const auto j = pos.j;
    const auto k = pos.k;

    return  object_xface_map[i    ][j][k] &&
            object_xface_map[i + 1][j][k] &&
            object_yface_map[i][j    ][k] &&
            object_yface_map[i][j + 1][k] &&
            object_zface_map[i][j][k    ] &&
            object_zface_map[i][j][k + 1];
}
*/

void Spacecraft::removeInnerParticle(Particle& p) const {
    if (isIncluded(p)) { p.makeInvalid(); }
}

void Spacecraft::distributeInnerParticleCharge(Particle& p) {
    if (isIncluded(p)) { 
        const auto id = p.typeId;
        const auto q = p.getCharge();
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
    cout << "charge before redist: " << q << endl;

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
        delta_rho = MPIw::Environment::Comms[name].sum(delta_rho);

        if (isMyCmat(i)) {
            const auto& target_pos = capacity_matrix_relation.at(i);
            rho[0][target_pos.i][target_pos.j][target_pos.k] += delta_rho;
        }
    }

    q = getTotalCharge(rho);
    cout << "charge after redist: " << q << endl;
    // update total charge for output
    total_charge = q;
}

//! 粒子放出関連
void Spacecraft::emitParticles(ParticleArray& parray) {
    cout << Environment::rankStr() << "called emitParticles." << endl;

    for(const auto& id : emit_particle_ids) {
        const auto emit_ptype_ptr = Environment::getEmissionParticleType(id);
        const auto max_amount = emit_ptype_ptr->getEmissionAmount();

        for(int i = 0; i < max_amount; ++i) {
            Particle p = emit_ptype_ptr->generateNewParticle();

            cout << p << endl;
            // parray.push_back()
        }
    }
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
    ost << "    name          : " << spc.name << endl;
    ost << "    cmat          : " << spc.num_cmat << endl;
    ost << "    emit particles: ";

    for(const auto& id : spc.emit_particle_ids) {
        ost << Environment::getParticleType(id)->getName() << ", " << endl;
    }

    return ost;
}


// Utility Functions for Objects
namespace ObjectUtils {
    ObjectDataFromFile getObjectNodesFromObjFile(const std::string& obj_file_name) {
        boost::filesystem::path p(obj_file_name);

        ObjectDataFromFile obj_data;

        if (boost::filesystem::exists(p)) {
            ObjectDefinedMap object_node_map(boost::extents[ Environment::nx ][ Environment::ny ][ Environment::nz ]);

            //! 初期化
            for(int i = 0; i < Environment::nx; ++i) {
                for (int j = 0; j < Environment::ny; ++j) {
                    for (int k = 0; k < Environment::nz; ++k) {
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
            using Face = std::array<int, 4>;
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
                            if (j != maxy && k != maxz) obj_data.faces.push_back({{0, i, j, k}});
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
                            if (i != maxx && k != maxz) obj_data.faces.push_back({{1, i, j, k}});
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
                            if (i != maxx && j != maxy) obj_data.faces.push_back({{2, i, j, k}});
                        }
                    }
                } else {
                    throw std::logic_error("Face type cannot be determinted by vertices position.");
                }
            }
        } else {
            std::string error_message = (format("File %s does not exist.") % obj_file_name).str();
            throw std::invalid_argument(error_message);
        }

        return obj_data;
    }
}
