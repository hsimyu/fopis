#include <tdpic.h>
#include <fstream>

using std::cout;
using std::endl;
using boost::format;

namespace Utils {
    // Normalizedのstatic変数の実体
    double Normalizer::x_unit = 1.0;
    double Normalizer::t_unit = 1.0;
    double Normalizer::e_unit = 1.0;

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

    threeDArray* create3DArray(const int nx, const int ny, const int nz) {
        threeDArray* x = new threeDArray;

        // resize to [nx][ny][nz]
        threeDArray::extent_gen extents;
        x->resize(extents[nx][ny][nz]);

        return x;
    }

    void delete3DArray(threeDArray* x) {
        delete x;
    }

    void clearBoundaryValues(threeDArray* x, const int nx, const int ny, const int nz) {
        for(int i = 0; i < nx; i += nx - 1) {
            for(int j = 0; j < ny; ++j){
                for(int k = 0; k < nz; ++k){
                    (*x)[i][j][k] = 0.0;
                }
            }
        }

        for(int j = 0; j < ny; j += ny - 1){
            for(int i = 1; i < nx - 1; ++i) {
                for(int k = 0; k < nz; ++k){
                    (*x)[i][j][k] = 0.0;
                }
            }
        }

        for(int k = 0; k < nz; k += nz - 1){
            for(int j = 1; j < ny - 1; ++j){
                for(int i = 1; i < nx - 1; ++i) {
                    (*x)[i][j][k] = 0.0;
                }
            }
        }
    }

    void convert3Dto1Darray(threeDArray* x3D, const int nx, const int ny, const int nz, double* x1D){
        // convert to 1D array
        // Fortran-based indicing
        for(int k = 0; k < nz; ++k){
            for(int j = 0; j < ny; ++j){
                for(int i = 0; i < nx; ++i){
                    x1D[i + j*nx + k*nx*ny] = (*x3D)[i][j][k];
                }
            }
        }
    }

    void convert1Dto3Darray(double* x1D, const int nx, const int ny, const int nz, threeDArray* x3D){
        // convert to 3D array
        for(int i = 0; i < nx; ++i){
            for(int j = 0; j < ny; ++j){
                for(int k = 0; k < nz; ++k){
                    (*x3D)[i][j][k] = x1D[i + j*nx + k*nx*ny];
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

        ifs.open(filename, std::ios::in);

        while(!ifs.eof()){
            std::getline(ifs, ifs_buffer);
            res += ifs_buffer;
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
}
