#include <tdpic.h>
#include <fstream>
#include <silo.h>

#include <mpi.h>
#include <pmpio.h>

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
        // std::string dirname = (format("block%d") % blocknum).str();
        // DBMkDir(file, dirname.c_str());
        // DBSetDir(file, dirname.c_str());

        tdArray& tdArray = g->getField()->getPhi();

        // dimensions
        // glue cellも出力
        int dimensions[3];
        dimensions[0] = g->getNX() + 2;
        dimensions[1] = g->getNY() + 2;
        dimensions[2] = g->getNZ() + 2;

        // names of the coordinates
        char* coordnames[3];
        coordnames[0] = const_cast<char*>("x");
        coordnames[1] = const_cast<char*>("y");
        coordnames[2] = const_cast<char*>("z");

        // names of the variables
        char* varnames[1];
        varnames[0] = "potential";

        // the array of coordinate arrays
        float** coordinates = g->getMeshNodes(dim);

        // make options list
        DBoptlist* optList = DBMakeOptlist(2);
        char* unit = "V";
        DBAddOption(optList, DBOPT_UNITS, unit);
        int major_order = 1;
        DBAddOption(optList, DBOPT_MAJORORDER, &major_order); // column-major (Fortran) order

        double* vars[] = {tdArray.data()};

        // directoryが1階層下がった分を考慮
        // std::string mesh_path = "../" + meshname;
        // std::string var_path = "../" + meshname;

        DBPutQuadmesh(file, meshname, coordnames, coordinates, dimensions, dim, DB_FLOAT, DB_COLLINEAR, NULL);
        DBPutQuadvar(file, varname, meshname, 1, varnames, vars, dimensions, dim, NULL, 0, DB_DOUBLE, DB_NODECENT, optList);

        // Free optList
        DBFreeOptlist(optList);

        // free mesh nodes
        delete [] coordinates[0];
        delete [] coordinates[1];
        delete [] coordinates[2];
        delete [] coordinates;

        // DBSetDir(file, "..");

    }

    // void writeData(Grid* g, int timestep) {
    //     std::string filename = (format("potential%d_%04d.silo") % timestep % MPI::Environment::rank).str();
    //     DBfile* file = DBCreate(filename.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);
    //
    //     int total_blocknum = MPI::Environment::numprocs;
    //     char* meshnames[total_blocknum];
    //     char* varnames[total_blocknum];
    //
    //     for(int i = 0; i < total_blocknum; ++i) {
    //         meshnames[i] = ((format("mesh%d") % i).str()).c_str();
    //         varnames[i] = ((format("var%d") % i).str()).c_str();
    //     }
    //
    //     writeMultimesh(file, total_blocknum, meshnames, varnames);
    //
    //     for(int i = 0; i < total_blocknum; ++i){
    //         writeBlock(file, i, g, meshnames[i], varnames[i]);
    //     }
    //
    //     DBClose(file);
    // }

    void writeDataInParallel(Grid* g, int timestep, std::string dataTypeName) {
        const int maxIOUnit = 2;
        int numfiles = (MPI::Environment::numprocs <= maxIOUnit) ? MPI::Environment::numprocs : maxIOUnit;

        PMPIO_baton_t* bat = PMPIO_Init(numfiles, PMPIO_WRITE, MPI_COMM_WORLD, 1000, PMPIO_DefaultCreate, PMPIO_DefaultOpen, PMPIO_DefaultClose, NULL);
        int groupRank = PMPIO_GroupRank(bat, MPI::Environment::rank);
        int rankInGroup = PMPIO_RankInGroup(bat, MPI::Environment::rank);
        std::string filename = (format("data/%s_%d_%04d.silo") % dataTypeName % timestep % groupRank).str();
        std::string blockname = (format("block%04d") % rankInGroup).str();
        DBfile* file = static_cast<DBfile*>(PMPIO_WaitForBaton(bat, filename.c_str(), blockname.c_str()));

        PMPIO_HandOffBaton(bat, file);
        PMPIO_Finish(bat);
    }

    void print3DArray(const tdArray& data, const int nx, const int ny, const int nz){
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
