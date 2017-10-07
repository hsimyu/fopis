#include "utils.hpp"
#include "mpiw.hpp"
#include "environment.hpp"
#include <fstream>
#include <stdexcept>
#include <boost/filesystem.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <Eigen/Core>

namespace Utils {
    TimeCounter* TimeCounter::instance = nullptr;

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

    void initialize3DArray(tdArray& x) {
        #pragma omp parallel for shared(x)
        for(int i = 0; i < x.shape()[0]; ++i) {
            for(int j = 0; j < x.shape()[1]; ++j) {
                for(int k = 0; k < x.shape()[2]; ++k) {
                    x[i][j][k] = 0.0;
                }
            }
        }
    }

    void initializeRhoArray(std::vector<tdArray>& x) {
        #pragma omp parallel shared(x)
        for(int id = 0; id < x.size(); ++id) {
            #pragma omp for
            for(int i = 0; i < x[0].shape()[0]; ++i) {
                for(int j = 0; j < x[0].shape()[1]; ++j) {
                    for(int k = 0; k < x[0].shape()[2]; ++k) {
                        x[id][i][j][k] = 0.0;
                    }
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

    bool isExistingFile(const std::string& file_name) {
        boost::filesystem::path p(file_name);
        return boost::filesystem::exists(p);
    }

    bool isExistingDirectory(const std::string& dir_name) {
        boost::filesystem::path p(dir_name);
        return boost::filesystem::is_directory(p);
    }

    void createDir(std::string dir_name) {
        if(!isExistingDirectory(dir_name)) {
            boost::filesystem::path dir(dir_name);
            boost::filesystem::create_directory(dir);
        }
    }

    std::string extractFileName(const std::string& target) {
        std::vector<std::string> strings = split(target, '/');
        return strings[ strings.size() - 1 ];
    }

    std::string readFile(const std::string& file_name){
        if(isExistingFile(file_name)) {
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

    //! convert picojson::array to std::vector
    std::vector<double> convertPicoJSONArrayToVectorDouble(const picojson::array& pico_array) {
        std::vector<double> vect;
        for(const auto& v : pico_array) {
            vect.push_back(v.get<double>());
        }
        return vect;
    }

    std::vector<std::string> convertPicoJSONArrayToVectorString(const picojson::array& pico_array) {
        std::vector<std::string> vect;
        for(const auto& v : pico_array) {
            vect.push_back(v.to_str());
        }
        return vect;
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
}
