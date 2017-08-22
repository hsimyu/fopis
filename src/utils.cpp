#include "utils.hpp"
#include "mpiw.hpp"
#include "environment.hpp"
#include <fstream>
#include <regex>
#include <stdexcept>
#include <boost/filesystem.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace Utils {

    //! キャパシティ行列用だけに使うため、受け取った参照先を直接置き換える実装で良い
    void makeInvert(dMatrix& lhs) {
        namespace ublas = boost::numeric::ublas;

        dMatrix lhs_copy(lhs);
        dMatrix inv( ublas::identity_matrix<double>( lhs.size1() ) );
        ublas::permutation_matrix<> pm( lhs.size1() );

        ublas::lu_factorize(lhs, pm);
        ublas::lu_substitute(lhs, pm, inv);

        lhs.assign_temporary(inv);
    }

    void initializeTdarray(tdArray& x) {
        for(int i = 0; i < x.shape()[0]; ++i) {
            for(int j = 0; j < x.shape()[1]; ++j) {
                for(int k = 0; k < x.shape()[2]; ++k) {
                    x[i][j][k] = 0.0;
                }
            }
        }
    }

    void clearBoundaryValues(tdArray& x, const int nx, const int ny, const int nz) {
        for(int i = 0; i < nx; i += nx - 1) {
            for(int j = 0; j < ny; ++j){
                for(int k = 0; k < nz; ++k){
                    x[i][j][k] = 0.0;
                }
            }
        }

        for(int j = 0; j < ny; j += ny - 1){
            for(int i = 1; i < nx - 1; ++i) {
                for(int k = 0; k < nz; ++k){
                    x[i][j][k] = 0.0;
                }
            }
        }

        for(int k = 0; k < nz; k += nz - 1){
            for(int j = 1; j < ny - 1; ++j){
                for(int i = 1; i < nx - 1; ++i) {
                    x[i][j][k] = 0.0;
                }
            }
        }
    }

    int getAxisIndex(const AXIS axis) {
        if (axis == AXIS::x) {
            return 0;
        } else if (axis == AXIS::y) {
            return 1;
        } else if (axis == AXIS::z) {
            return 2;
        } else {
            throw std::invalid_argument( "Unknown axis type was passed." );
        }
    }

    int getLowOrUpIndex(const AXIS_SIDE low_or_up) {
        if (low_or_up == AXIS_SIDE::low) {
            return 0;
        } else if (low_or_up == AXIS_SIDE::up) {
            return 1;
        } else {
            throw std::invalid_argument( "Unknown low_or_up type was passed." );
        }
    }

    float* getTrueEdges2(tdArray const& xvalue, tdArray const& yvalue, tdArray const& zvalue){
        int nx = yvalue.shape()[0];
        int ny = xvalue.shape()[1];
        int nz = xvalue.shape()[2];

        const int c_indexing = 0;
        const int fortran_indexing = 1;
        const int indexing = c_indexing;

        const int length = (nx-2)*(ny-2)*(nz-2);
        const int yoffset = length;
        const int zoffset = 2 * length;
        float* x1D = new float[length];

        if(indexing == fortran_indexing) {
            for(int i = 1; i < nx - 1; ++i){
                for(int j = 1; j < ny - 1; ++j){
                    for(int k = 1; k < nz - 1; ++k){
                        if(i != (nx - 2)) {
                            x1D[(k-1) + (j-1)*(nz-2) + (i-1)*(nz-2)*(ny-2)] = static_cast<float>(xvalue[i][j][k]);
                        } else {
                            x1D[(k-1) + (j-1)*(nz-2) + (i-1)*(nz-2)*(ny-2)] = 0.0f;
                        }

                        if(j != (ny - 2)) {
                            x1D[yoffset + (k-1) + (j-1)*(nz-2) + (i-1)*(nz-2)*(ny-2)] = static_cast<float>(yvalue[i][j][k]);
                        } else {
                            x1D[yoffset + (k-1) + (j-1)*(nz-2) + (i-1)*(nz-2)*(ny-2)] = 0.0f;
                        }

                        if(k != (nz - 2)) {
                            x1D[zoffset + (k-1) + (j-1)*(nz-2) + (i-1)*(nz-2)*(ny-2)] = static_cast<float>(zvalue[i][j][k]);
                        } else {
                            x1D[zoffset + (k-1) + (j-1)*(nz-2) + (i-1)*(nz-2)*(ny-2)] = 0.0f;
                        }
                    }
                }
            }
        } else {
            for(int k = 1; k < nz - 1; ++k){
                for(int j = 1; j < ny - 1; ++j){
                    for(int i = 1; i < nx - 1; ++i){
                        if(i != (nx - 2)) {
                            x1D[(i-1) + (j-1)*(nx-2) + (k-1)*(nx-2)*(ny-2)] = static_cast<float>(xvalue[i][j][k]);
                        } else {
                            x1D[(i-1) + (j-1)*(nx-2) + (k-1)*(nx-2)*(ny-2)] = 0.0f;
                        }
                        if(j != (ny - 2)) {
                            x1D[yoffset + (i-1) + (j-1)*(nx-2) + (k-1)*(nx-2)*(ny-2)] = static_cast<float>(yvalue[i][j][k]);
                        } else {
                            x1D[yoffset + (i-1) + (j-1)*(nx-2) + (k-1)*(nx-2)*(ny-2)] = 0.0f;
                        }
                        if(k != (nz - 2)) {
                            x1D[zoffset + (i-1) + (j-1)*(nx-2) + (k-1)*(nx-2)*(ny-2)] = static_cast<float>(zvalue[i][j][k]);
                        } else {
                            x1D[zoffset + (i-1) + (j-1)*(nx-2) + (k-1)*(nx-2)*(ny-2)] = 0.0f;
                        }
                    }
                }
            }
        }

        return x1D;
    }

    float* getTrueEdges(const tdArray& x3D, const int axis){
        int nx = x3D.shape()[0];
        int ny = x3D.shape()[1];
        int nz = x3D.shape()[2];

        const int c_indexing = 0;
        const int fortran_indexing = 1;
        const int indexing = c_indexing;

        // Add extra slot to axis
        // switch(axis){
        //     case 0:
        //         nx += 1;
        //         break;
        //     case 1:
        //         ny += 1;
        //         break;
        //     case 2:
        //         nz += 1;
        //         break;
        //     default:
        //         throw std::invalid_argument("[ERROR] Unknown edge axis was passed to getTrueEdges.");
        //         break;
        // }
        //
        float* x1D = new float[(nx-2)*(ny-2)*(nz-2)];

        if(indexing == fortran_indexing) {
            for(int i = 1; i < nx - 1; ++i){
                for(int j = 1; j < ny - 1; ++j){
                    for(int k = 1; k < nz - 1; ++k){
                        switch(axis){
                            case 0:
                                if(i != (nx - 2)) {
                                    x1D[(k-1) + (j-1)*(nz-2) + (i-1)*(nz-2)*(ny-2)] = static_cast<float>(x3D[i][j][k]);
                                } else {
                                    x1D[(k-1) + (j-1)*(nz-2) + (i-1)*(nz-2)*(ny-2)] = 0.0f;
                                }
                                break;
                            case 1:
                                if(j != (ny - 2)) {
                                    x1D[(k-1) + (j-1)*(nz-2) + (i-1)*(nz-2)*(ny-2)] = static_cast<float>(x3D[i][j][k]);
                                } else {
                                    x1D[(k-1) + (j-1)*(nz-2) + (i-1)*(nz-2)*(ny-2)] = 0.0f;
                                }
                                break;
                            case 2:
                                if(k != (nz - 2)) {
                                    x1D[(k-1) + (j-1)*(nz-2) + (i-1)*(nz-2)*(ny-2)] = static_cast<float>(x3D[i][j][k]);
                                } else {
                                    x1D[(k-1) + (j-1)*(nz-2) + (i-1)*(nz-2)*(ny-2)] = 0.0f;
                                }
                                break;
                            default:
                                throw std::invalid_argument("[ERROR] Unknown edge axis was passed to getTrueEdges.");
                                break;
                        }
                    }
                }
            }
        } else {
            for(int k = 1; k < nz - 1; ++k){
                for(int j = 1; j < ny - 1; ++j){
                    for(int i = 1; i < nx - 1; ++i){
                        switch(axis){
                            case 0:
                                if(i != (nx - 2)) {
                                    x1D[(i-1) + (j-1)*(nx-2) + (k-1)*(nx-2)*(ny-2)] = static_cast<float>(x3D[i][j][k]);
                                } else {
                                    x1D[(i-1) + (j-1)*(nx-2) + (k-1)*(nx-2)*(ny-2)] = 0.0f;
                                }
                                break;
                            case 1:
                                if(j != (ny - 2)) {
                                    x1D[(i-1) + (j-1)*(nx-2) + (k-1)*(nx-2)*(ny-2)] = static_cast<float>(x3D[i][j][k]);
                                } else {
                                    x1D[(i-1) + (j-1)*(nx-2) + (k-1)*(nx-2)*(ny-2)] = 0.0f;
                                }
                                break;
                            case 2:
                                if(k != (nz - 2)) {
                                    x1D[(i-1) + (j-1)*(nx-2) + (k-1)*(nx-2)*(ny-2)] = static_cast<float>(x3D[i][j][k]);
                                } else {
                                    x1D[(i-1) + (j-1)*(nx-2) + (k-1)*(nx-2)*(ny-2)] = 0.0f;
                                }
                                break;
                            default:
                                throw std::invalid_argument("[ERROR] Unknown edge axis was passed to getTrueEdges.");
                                break;
                        }
                    }
                }
            }
        }

        return x1D;
    }

    float* getTrueFaces(const tdArray& x3D, const int axis){
        int nx = x3D.shape()[0];
        int ny = x3D.shape()[1];
        int nz = x3D.shape()[2];

        enum class INDEXING { C, FORTRAN };

        const auto indexing = INDEXING::C;

        // Add extra slot to axis
        switch(axis){
            case 0:
                ny += 1; nz += 1;
                break;
            case 1:
                nx += 1; nz += 1;
                break;
            case 2:
                ny += 1; nz += 1;
                break;
            default:
                throw std::invalid_argument("[ERROR] Unknown edge axis was passed to getTrueEdges.");
                break;
        }

        float* x1D = new float[(nx-2)*(ny-2)*(nz-2)];

        if(indexing == INDEXING::FORTRAN) {
            for(int i = 1; i < nx - 1; ++i){
                for(int j = 1; j < ny - 1; ++j){
                    for(int k = 1; k < nz - 1; ++k){
                        switch(axis){
                            case 0:
                                if(j != (ny - 2) && k != (nz - 2)) {
                                    x1D[(k-1) + (j-1)*(nz-2) + (i-1)*(nz-2)*(ny-2)] = static_cast<float>(x3D[i][j][k]);
                                }
                                break;
                            case 1:
                                if(i != (nx - 2) && k != (nz - 2)) {
                                    x1D[(k-1) + (j-1)*(nz-2) + (i-1)*(nz-2)*(ny-2)] = static_cast<float>(x3D[i][j][k]);
                                }
                                break;
                            case 2:
                                if(i != (nx - 2) && j != (ny - 2)) {
                                    x1D[(k-1) + (j-1)*(nz-2) + (i-1)*(nz-2)*(ny-2)] = static_cast<float>(x3D[i][j][k]);
                                }
                                break;
                            default:
                                throw std::invalid_argument("[ERROR] Unknown edge axis was passed to getTrueEdges.");
                                break;
                        }
                    }
                }
            }
        } else {
            for(int k = 1; k < nz - 1; ++k){
                for(int j = 1; j < ny - 1; ++j){
                    for(int i = 1; i < nx - 1; ++i){
                        switch(axis){
                            case 0:
                                if(j != (ny - 2) && k != (nz - 2)) {
                                    x1D[(i-1) + (j-1)*(nx-2) + (k-1)*(nx-2)*(ny-2)] = static_cast<float>(x3D[i][j][k]);
                                }
                                break;
                            case 1:
                                if(i != (nx - 2) && k != (nz - 2)) {
                                    x1D[(i-1) + (j-1)*(nx-2) + (k-1)*(nx-2)*(ny-2)] = static_cast<float>(x3D[i][j][k]);
                                }
                                break;
                            case 2:
                                if(i != (nx - 2) && j != (ny - 2)) {
                                    x1D[(i-1) + (j-1)*(nx-2) + (k-1)*(nx-2)*(ny-2)] = static_cast<float>(x3D[i][j][k]);
                                }
                                break;
                            default:
                                throw std::invalid_argument("[ERROR] Unknown edge axis was passed to getTrueEdges.");
                                break;
                        }
                    }
                }
            }
        }

        return x1D;
    }

    void convert1Dto3Darray(double* x1D, const int nx, const int ny, const int nz, tdArray& x3D){
        // convert to 3D array
        for(int i = 0; i < nx; ++i){
            for(int j = 0; j < ny; ++j){
                for(int k = 0; k < nz; ++k){
                    x3D[i][j][k] = x1D[i + j*nx + k*nx*ny];
                }
            }
        }
    }

    std::string prettyMemoryString(double mem){
        std::string suffix = "B";
        if(mem > 1048.0){
            mem /= 1048.0;
            suffix = "kB";
        }
        if(mem > 1048.0){
            mem /= 1048.0;
            suffix = "MB";
        }
        if(mem > 1048.0){
            mem /= 1048.0;
            suffix = "GB";
        }

        return (format("%-6.2f") % mem).str() + suffix;
    }

    std::string readFile(const std::string& file_name){
        boost::filesystem::path p(file_name);

        if(boost::filesystem::exists(p)) {
            std::ifstream ifs;
            std::string res;
            std::string ifs_buffer;

            ifs.open(file_name, std::ios::in);

            while(!ifs.eof()){
                std::getline(ifs, ifs_buffer);
                res += ifs_buffer;
            }
            return res;
        } else {
            std::string error_message = (format("File %s does not exist.") % file_name).str();
            throw std::invalid_argument(error_message);
            MPIw::Environment::abort(1);
        }
    }

    picojson::value::object readJSONFile(const std::string& filename){
        std::string json = readFile(filename);

        picojson::value v;
        std::string error = picojson::parse(v, json);
        picojson::value::object& o = v.get<picojson::object>();

        return o;
    }

    void createDir(std::string dirname) {
        boost::filesystem::path dir(dirname);
        if(!boost::filesystem::is_directory(dir)) {
            boost::filesystem::create_directory(dir);
        }
    }

    std::vector<std::string> split(const std::string& target, char delim) {
        std::stringstream ss(target);
        std::string buffer;
        std::vector<std::string> res;

        while( std::getline(ss, buffer, delim) ) {
            res.push_back(buffer);
        }

        return res;
    }

    ObjectNodes getObjectNodesFromObjFile(const std::string& obj_file_name) {
        boost::filesystem::path p(obj_file_name);
        ObjectNodes temp_obj_node_array;

        if (boost::filesystem::exists(p)) {
            ObjectDefinedMap temp_object_map(boost::extents[ Environment::nx ][ Environment::ny ][ Environment::nz ]);
            //! 初期化
            for(int i = 0; i < Environment::nx; ++i) {
                for (int j = 0; j < Environment::ny; ++j) {
                    for (int k = 0; k < Environment::nz; ++k) {
                        temp_object_map[i][j][k] = false;
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
                    const auto& verts = split(buffer, ' ');

                    Vertex vertex;

                    vertex[0] = static_cast<int>(std::stoi(verts[1]) + Environment::nx / 2);
                    vertex[1] = static_cast<int>(std::stoi(verts[2]) + Environment::ny / 2);
                    vertex[2] = static_cast<int>(std::stoi(verts[3]) + Environment::nz / 2);

                    vertices.emplace_back(vertex);
                } else if (std::regex_match(buffer, re_face)) {
                    const auto& fs = split(buffer, ' ');
                    Face face_numbers;
                    for(int i = 1; i < 5; ++i) {
                        const auto& splitted_face = split(fs[i], '/');
                        face_numbers[i - 1] = static_cast<int>(std::stoi(splitted_face[0]));
                    }
                    faces.emplace_back(face_numbers);
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
                            if (!temp_object_map[i][j][k]) {
                                temp_obj_node_array[num_cmat] = {{i, j, k}};
                                ++num_cmat;
                                temp_object_map[i][j][k] = true;
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
                            if (!temp_object_map[i][j][k]) {
                                temp_obj_node_array[num_cmat] = {{i, j, k}};
                                ++num_cmat;
                                temp_object_map[i][j][k] = true;
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
                            if (!temp_object_map[i][j][k]) {
                                temp_obj_node_array[num_cmat] = {{i, j, k}};
                                ++num_cmat;
                                temp_object_map[i][j][k] = true;
                            }
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

        return temp_obj_node_array;
    }
}
