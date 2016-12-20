#include <tdpic.h>
#include <fstream>
#include <H5Cpp.h> // C++ HDF5 Library
#include <write_hdf5.hpp>

using std::cout;
using std::cin;
using std::endl;
using boost::format;

namespace IO {
    const int MAX_NAME_LENGTH = 100;

    typedef struct {
        int age;
        char sex;
        char name[MAX_NAME_LENGTH];
        float height;
    } PersonalInformation;

    void writeFieldData(Grid* g, const std::string filename) {
        auto shape = g->getField()->getPhi().shape();
        threeDArray temp(boost::extents[ shape[0] ][ shape[1] ][ shape[2] ]);

        //! copy to temp array
        //! To make C-like ordering array
        temp = g->getField()->getPhi();

        H5::H5File file(filename, H5F_ACC_TRUNC);
        write_hdf5(file, "potential", temp);
    }

    void print3DArray(const threeDArray& data, const int nx, const int ny, const int nz){
        for (int k = 0; k < nz; ++k ) {
            cout << "[z:" << k << "] " << endl;

            for ( int i = 0 ; i < nx; ++i ) {
                    if(i == 0) {
                        cout << "     [x/y]";
                        for ( int j = 0 ; j < ny; ++j ) {
                            cout << "[" << j << "]";
                        }
                        cout << endl;
                    }
                    cout << "     [" << i << "]  ";
                for ( int j = 0 ; j < ny; ++j ) {
                    cout << " " << data[i][j][k] << " ";
                }
                cout << endl;
            }
            cout << endl;
        }
    }

    void outputParticlePositions(const Environment* env, const ParticleArray& parray, std::string filename){
        std::ofstream ofs(filename);

        for(int id = 0; id < env->num_of_particle_types; ++id){

            ofs << "## " << env->ptype[id].getName() << endl;

            for(int i = 0; i < parray[id].size(); ++i){
                ofs << format("%9.4f %9.4f %9.4f") % parray[id][i].getX() % parray[id][i].getY() % parray[id][i].getZ();
                ofs << format("%13.4e %13.4e %13.4e") % parray[id][i].getVX() % parray[id][i].getVY() % parray[id][i].getVZ();
                ofs << endl;
            }

            ofs << endl << endl;
        }
    }
}
