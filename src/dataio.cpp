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

    void writeBlock(DBfile* file, Grid* g, char* meshname, char* varname){
        // dimension
        const int dim = 3;

        double* tdArray = Utils::getTrueCells(g->getField()->getPhi());

        // dimensions
        // glue cellも出力
        int dimensions[3];
        dimensions[0] = g->getNX();
        dimensions[1] = g->getNY();
        dimensions[2] = g->getNZ();

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

        // make options list for mesh
        DBoptlist* optListMesh = DBMakeOptlist(1);

        // make options list for var
        DBoptlist* optListVar = DBMakeOptlist(2);
        char* unit = "V";
        DBAddOption(optListVar, DBOPT_UNITS, unit);
        int major_order = 1;
        DBAddOption(optListVar, DBOPT_MAJORORDER, &major_order); // column-major (Fortran) order

        double* vars[] = {tdArray};

        DBPutQuadmesh(file, meshname, coordnames, coordinates, dimensions, dim, DB_FLOAT, DB_COLLINEAR, NULL);
        DBPutQuadvar(file, varname, meshname, 1, varnames, vars, dimensions, dim, NULL, 0, DB_DOUBLE, DB_NODECENT, optListVar);

        // Free optList
        DBFreeOptlist(optListMesh);
        DBFreeOptlist(optListVar);

        // free mesh nodes
        delete [] coordinates[0];
        delete [] coordinates[1];
        delete [] coordinates[2];
        delete [] coordinates;
    }

    void writeDataInParallel(Grid* g, int timestep, std::string dataTypeName) {
        const int maxIOUnit = 8;
        int numfiles = (MPI::Environment::numprocs <= maxIOUnit) ? MPI::Environment::numprocs : maxIOUnit;

        PMPIO_baton_t* bat = PMPIO_Init(numfiles, PMPIO_WRITE, MPI_COMM_WORLD, 1000, PMPIO_DefaultCreate, PMPIO_DefaultOpen, PMPIO_DefaultClose, NULL);
        int groupRank = PMPIO_GroupRank(bat, MPI::Environment::rank);
        int rankInGroup = PMPIO_RankInGroup(bat, MPI::Environment::rank);
        std::string filename = (format("data/%s_%04d_%04d.silo") % dataTypeName % groupRank % timestep).str();
        std::string blockname = (format("block%04d") % rankInGroup).str();
        DBfile* file = static_cast<DBfile*>(PMPIO_WaitForBaton(bat, filename.c_str(), blockname.c_str()));

        if(rankInGroup == 0 && groupRank == 0) {
            int total_blocknum = MPI::Environment::numprocs;
            char** meshnames = new char*[total_blocknum];
            char** varnames = new char*[total_blocknum];

            for(int i = 0; i < total_blocknum; ++i) {
                int tmpGroupRank = PMPIO_GroupRank(bat, i);
                int tmpRankInGroup = PMPIO_RankInGroup(bat, i);
                std::string tmpfilename = (format("%s_%04d_%04d.silo") % dataTypeName % tmpGroupRank % timestep).str();

                std::string tmpstring;
                if(filename == tmpstring) {
                    tmpstring = (format("/block%04d/mesh") % tmpRankInGroup).str();
                } else {
                    tmpstring = (format("%s:/block%04d/mesh") % tmpfilename % tmpRankInGroup).str();
                }
                meshnames[i] = new char[tmpstring.size() + 1];
                std::strcpy(meshnames[i], tmpstring.c_str());

                if(filename == tmpstring) {
                    tmpstring = (format("/block%04d/%s") % tmpRankInGroup % dataTypeName).str();
                } else {
                    tmpstring = (format("%s:/block%04d/%s") % tmpfilename % tmpRankInGroup % dataTypeName).str();
                }
                varnames[i] = new char[tmpstring.size() + 1];
                std::strcpy(varnames[i], tmpstring.c_str());
            }

            writeMultimesh(file, total_blocknum, meshnames, varnames);

            for(int i = 0; i < total_blocknum; ++i) {
                delete [] meshnames[i];
                delete [] varnames[i];
            }
            delete [] meshnames;
            delete [] varnames;
        }

        char* meshname = "mesh";
        char* varname = const_cast<char*>(dataTypeName.c_str());
        writeBlock(file, g, meshname, varname);

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
