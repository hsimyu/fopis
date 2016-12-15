#include <tdpic.h>
#include <fstream>

using std::cout;
using std::endl;
using boost::format;

namespace Utils {
    double*** create3DArray(const int nx, const int ny, const int nz) {
        double*** x;

        x = new double**[nx];
        x[0] = new double*[nx * ny];
        x[0][0] = new double[nx * ny * nz];

        for(int i = 0; i < nx; ++i) {
            x[i] = x[0] + i * ny;
            for(int j = 0; j  < ny; ++j){
                x[i][j] = x[0][0] + i * ny * nz + j * nz;
            }
        }

        for(int i = 0; i < nx; ++i) {
            for(int j = 0; j  < ny; ++j){
                for(int k = 0; k  < nz; ++k){
                    x[i][j][k] = 0.0;
                }
            }
        }

        return x;
    }

    void delete3DArray(double*** x) {
        delete [] x[0][0];
        delete [] x[0];
        delete [] x;
    }

    void clearBoundaryValues(double*** x, const int nx, const int ny, const int nz) {
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

    void convert3Dto1Darray(double*** x3D, const int nx, const int ny, const int nz, double* x1D){
        // convert to 1D array
        // Fortran-based indicing
        for(int k = 0; k < nz; ++k){
            for(int j = 0; j < ny; ++j){
                for(int i = 0; i < nx; ++i){
                    x1D[i + j*nx + k*nx*ny] = x3D[i][j][k];
                }
            }
        }
    }

    void convert1Dto3Darray(double* x1D, const int nx, const int ny, const int nz, double*** x3D){
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
