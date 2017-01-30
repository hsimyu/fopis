#include <tdpic.h>
#include <fstream>
#include <boost/filesystem.hpp>

namespace Utils {
    // Normalizedのstatic変数の実体
    //! @note: Normalizer は基本的にRoot Gridのグリッド幅しか保持していない
    double Normalizer::x_unit = 1.0;
    double Normalizer::t_unit = 1.0;
    double Normalizer::e_unit = 1.0;
    double Normalizer::m_unit = 1.0;

    // Normalize Utilities
    double Normalizer::normalizeLength(double raw_x) {
        return raw_x / x_unit;
    }

    double Normalizer::unnormalizeLength(double normalized_x) {
        return normalized_x * x_unit;
    }

    double Normalizer::normalizeVelocity(double raw_v) {
        return t_unit * raw_v / x_unit;
    }

    double Normalizer::unnormalizeVelocity(double normalized_v) {
        return normalized_v * x_unit / t_unit;
    }

    double Normalizer::normalizeTime(double raw_t) {
        return raw_t / t_unit;
    }

    double Normalizer::unnormalizeTime(double normalized_t) {
        return normalized_t * t_unit;
    }

    double Normalizer::normalizeCharge(double raw_e) {
        return raw_e / e_unit;
    }

    double Normalizer::unnormalizeCharge(double normalized_e) {
        return normalized_e * e_unit;
    }

    double Normalizer::normalizeMass(double raw_mass) {
        return raw_mass / m_unit;
    }

    double Normalizer::unnormalizeMass(double normalized_mass) {
        return normalized_mass * m_unit;
    }

    //! kg*m^2/s^2 -> 1;
    double Normalizer::normalizeEnergy(double raw_energy) {
        return raw_energy * pow(t_unit, 2) / (m_unit * pow(x_unit, 2));
    }

    //! 1 -> kg*m^2/s^2;
    double Normalizer::unnormalizeEnergy(double normalized_energy) {
        return normalized_energy * m_unit * pow(x_unit, 2) / pow(t_unit, 2); //kg * m^2/s^2;
    }

    //! s^4*A^2/kg*m^3 == C^2*s^2/kg*m^3 -> 1;
    double Normalizer::normalizeEpsilon(double raw_epsilon) {
        return raw_epsilon * (pow(x_unit, 3) * m_unit) / (pow(e_unit, 2) * pow(t_unit, 2));
    }

    //! 1 -> C^2*s^2/kg*m^3;
    double Normalizer::unnormalizeEpsilon(double normalized_eps) {
        return normalized_eps * (pow(e_unit, 2) * pow(t_unit, 2)) / (pow(x_unit, 3) * m_unit);
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

    //! for DATA IO
    float* getTrueCells(const tdArray& x3D){
        int nx = x3D.shape()[0];
        int ny = x3D.shape()[1];
        int nz = x3D.shape()[2];
        float* x1D = new float[(nx-2)*(ny-2)*(nz-2)];

        //! Fortran-based indicing
        //! without glue cells
        for(int k = 1; k < nz - 1; ++k){
            for(int j = 1; j < ny - 1; ++j){
                for(int i = 1; i < nx - 1; ++i){
                    x1D[(i-1) + (j-1)*(nx-2) + (k-1)*(nx-2)*(ny-2)] = static_cast<float>(x3D[i][j][k]);
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
        const int fortran_indexing = 0;
        const int indexing = c_indexing;

        // Add extra slot to axis
        switch(axis){
            case 0:
                nx += 1;
                break;
            case 1:
                ny += 1;
                break;
            case 2:
                nz += 1;
                break;
            default:
                throw std::invalid_argument("[ERROR] Unknown edge axis was passed to getTrueEdges.");
                break;
        }

        float* x1D = new float[(nx-2)*(ny-2)*(nz-2)];

        if(indexing == c_indexing) {
            for(int i = 1; i < nx - 1; ++i){
                for(int j = 1; j < ny - 1; ++j){
                    for(int k = 1; k < nz - 1; ++k){
                        switch(axis){
                            case 0:
                                if(i != (nx - 2)) {
                                    x1D[(k-1) + (j-1)*(nz-2) + (i-1)*(nz-2)*(ny-2)] = static_cast<float>(x3D[i][j][k]);
                                }
                                break;
                            case 1:
                                if(j != (ny - 2)) {
                                    x1D[(k-1) + (j-1)*(nz-2) + (i-1)*(nz-2)*(ny-2)] = static_cast<float>(x3D[i][j][k]);
                                }
                                break;
                            case 2:
                                if(k != (nz - 2)) {
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
                                if(i != (nx - 2)) {
                                    x1D[(i-1) + (j-1)*(nx-2) + (k-1)*(nx-2)*(ny-2)] = static_cast<float>(x3D[i][j][k]);
                                }
                                break;
                            case 1:
                                if(j != (ny - 2)) {
                                    x1D[(i-1) + (j-1)*(nx-2) + (k-1)*(nx-2)*(ny-2)] = static_cast<float>(x3D[i][j][k]);
                                }
                                break;
                            case 2:
                                if(k != (nz - 2)) {
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

    float* getTrueFaces(const tdArray& x3D, const int axis){
        int nx = x3D.shape()[0];
        int ny = x3D.shape()[1];
        int nz = x3D.shape()[2];

        const int c_indexing = 0;
        const int fortran_indexing = 0;
        const int indexing = c_indexing;

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

        if(indexing == c_indexing) {
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

    static std::string computeMemory(double mem){
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

        return (format("%6.2f") % mem).str() + suffix;
    }

    void printTotalMemory(const ParticleType& pinfo){
        cout << "Memory Info:" << endl;
        printParticleMemory(pinfo);
    }

    void printParticleMemory(const ParticleType& pinfo){
        static const double memory_per_particle = 50.0;

        double pmem = static_cast<double>(pinfo.getTotalNumber()) * memory_per_particle;

        cout << "    Particle[" << pinfo.getTotalNumber() << "]:";
        cout << computeMemory(pmem) << endl;
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
