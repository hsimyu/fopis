#include <regex>
#include "global.hpp"
#include "environment.hpp"
#include "spacecraft.hpp"
#include "utils.hpp"

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
