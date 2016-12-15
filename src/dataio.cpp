#include <tdpic.h>
#include <fstream>

using std::cout;
using std::cin;
using std::endl;
using boost::format;

namespace IO {
    void print3DArray(double*** const data, const int nx, const int ny, const int nz){
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
                ofs << format("%8.3f %8.3f %8.3f") % parray[id][i].getX() % parray[id][i].getY() % parray[id][i].getZ();
                ofs << format("%8.3f %8.3f %8.3f") % parray[id][i].getVX() % parray[id][i].getVY() % parray[id][i].getVZ();
                ofs << endl;
            }

            ofs << endl << endl;
        }
    }
}
