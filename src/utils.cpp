#include "utils.hpp"
#include "mpiw.hpp"
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


    std::string readFile(const std::string& filename){
        std::ifstream ifs;
        std::string res;
        std::string ifs_buffer;

        boost::filesystem::path p(filename);

        if(boost::filesystem::exists(p)) {
            ifs.open(filename, std::ios::in);

            while(!ifs.eof()){
                std::getline(ifs, ifs_buffer);
                res += ifs_buffer;
            }
        } else {
            throw std::invalid_argument("[ERROR] input.json does not exist.");
            MPIw::Environment::exitWithFinalize(1);
        }
        return res;
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
}
