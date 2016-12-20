#include <tdpic.h>
#include <fstream>
#include <H5Cpp.h> // C++ HDF5 Library

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

    void testHDF5Write(){
        const std::string filename = "test.h5";
        const std::string dataset_name = "Matrix in file";

        PersonalInformation person_list[] = {
            { 18, 'M', "Mary",  152.0   },
            { 32, 'F', "Tom",   178.6   },
            { 29, 'M', "Tarou", 166.6   }
        };

        int length = sizeof(person_list) / sizeof(PersonalInformation); // => 3

        hsize_t dim[1];
        dim[0] = length;

        int rank = sizeof(dim) / sizeof(hsize_t); // => 1

        H5::CompType mtype(sizeof(PersonalInformation));
        mtype.insertMember("Age", HOFFSET(PersonalInformation, age), H5::PredType::NATIVE_INT);
        mtype.insertMember("Sex", HOFFSET(PersonalInformation, sex), H5::PredType::C_S1);
        mtype.insertMember("Name", HOFFSET(PersonalInformation, name), H5::StrType(H5::PredType::C_S1, MAX_NAME_LENGTH));
        mtype.insertMember("Height", HOFFSET(PersonalInformation, height), H5::PredType::NATIVE_FLOAT);

        // data space
        H5::DataSpace space(rank, dim);

        // file pointer
        H5::H5File* file = new H5::H5File(filename, H5F_ACC_TRUNC);

        // data set
        H5::DataSet* dataset = new H5::DataSet(file->createDataSet(dataset_name, mtype, space));

        dataset->write(person_list, mtype);

        delete dataset;
        delete file;
    }

    void print3DArray(const threeDArray* data, const int nx, const int ny, const int nz){
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
                    cout << " " << (*data)[i][j][k] << " ";
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
