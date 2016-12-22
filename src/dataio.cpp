#include <tdpic.h>
#include <fstream>
#include <silo.h>

using std::cout;
using std::cin;
using std::endl;
using boost::format;

namespace IO {
    void writeMultimesh(DBfile* file, int total_blocknum, char** meshnames, char** varnames) {
        // make options list
        DBoptlist* optList = DBMakeOptlist(1);
        int meshtype = DB_QUAD_RECT;
        DBAddOption(optList, DBOPT_MB_BLOCK_TYPE, &meshtype);
        DBPutMultimesh(file, "multimesh", total_blocknum, meshnames, NULL, optList);

        int vartype = DB_QUADVAR;
        DBAddOption(optList, DBOPT_MB_BLOCK_TYPE, &vartype);
        DBPutMultivar(file, "multivar", total_blocknum, varnames, NULL, optList);
        DBFreeOptlist(optList);
    }

    void writeBlock(DBfile* file, int blocknum, Grid* g, char* meshname, char* varname){
        // dimension
        const int dim = 3;

        // Create new directory for a block
        std::string dirname = (format("block%d") % blocknum).str();
        DBMkDir(file, dirname.c_str());
        DBSetDir(file, dirname.c_str());

        threeDArray& tdArray = g->getField()->getPhi();
        int nx = tdArray.shape()[0];
        int ny = tdArray.shape()[1];
        int nz = tdArray.shape()[2];

        // dimensions
        int dimensions[dim];
        dimensions[0] = nx;
        dimensions[1] = ny;
        dimensions[2] = nz;

        // names of the coordinates
        char* coordnames[dim];
        coordnames[0] = const_cast<char*>("x");
        coordnames[1] = const_cast<char*>("y");
        coordnames[2] = const_cast<char*>("z");

        // names of the variables
        char* varnames[1];
        varnames[0] = "potential";

        // generate the coordinate arrays
        double dx = g->getDX();
        float node_x[nx];
        float node_y[ny];
        float node_z[nz];

        double base = g->getBaseX();
        for(int i = 0; i < nx; ++i) {
            node_x[i] = base + dx * i;
        }
        base = g->getBaseY();
        for(int i = 0; i < ny; ++i) {
            node_y[i] = base + dx * i;
        }
        base = g->getBaseZ();
        for(int i = 0; i < nz; ++i) {
            node_z[i] = base + dx * i;
        }

        // the array of coordinate arrays
        float* coordinates[dim];
        coordinates[0] = node_x;
        coordinates[1] = node_y;
        coordinates[2] = node_z;

        // make options list
        DBoptlist* optList = DBMakeOptlist(2);
        char* unit = "V";
        DBAddOption(optList, DBOPT_UNITS, unit);
        int major_order = 1;
        DBAddOption(optList, DBOPT_MAJORORDER, &major_order); // column-major (Fortran) order

        double* vars[] = {tdArray.data()};

        DBPutQuadmesh(file, meshname, coordnames, coordinates, dimensions, dim, DB_FLOAT, DB_COLLINEAR, NULL);
        DBPutQuadvar(file, varname, meshname, 1, varnames, vars, dimensions, dim, NULL, 0, DB_DOUBLE, DB_NODECENT, optList);

        // Free optList
        DBFreeOptlist(optList);

        DBSetDir(file, "..");

    }

    void writeData(Grid* g, int timestep) {
        std::string filename = (format("test%d.silo") % timestep).str();
        DBfile* file = DBCreate(filename.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);

        int total_blocknum = 1;
        char* meshnames[total_blocknum];
        char* varnames[total_blocknum];
        meshnames[0] = "mesh1";
        varnames[0] = "var1";

        writeMultimesh(file, total_blocknum, meshnames, varnames);

        for(int i = 0; i < total_blocknum; ++i){
            writeBlock(file, i, g, meshnames[i], varnames[i]);
        }

        DBClose(file);
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
