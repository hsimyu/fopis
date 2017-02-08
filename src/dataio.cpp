#include "global.hpp"
#include "environment.hpp"
#include "dataio.hpp"
#include "grid.hpp"
#include "mpiw.hpp"
#include "utils.hpp"
#include <numeric>
#include <silo.h>
#include <mpi.h>
#include <pmpio.h>

namespace IO {
    void writeGroupelMap(
            DBfile* file, Grid* root_grid, const int maxAMRLevel, const int numOfPatches,
            int* numOfPatchesOnLevel, std::vector< std::vector<int> >& idMap, std::map<int, std::vector<int> >& childMap)
    {
        // ARMの最大レベルと、パッチの総数、各レベルでのパッチ数は事前にMPI通信で確認しておく
        // Mrgvarは各ファイルに分割して記入する
        // -- make level map --
        int* levelSegTypes = new int[maxAMRLevel];
        for(int i = 0; i < maxAMRLevel; i++){
            levelSegTypes[i] = DB_BLOCKCENT;
        }

        int** segData = new int*[maxAMRLevel];
        // numPatchesOnLevelもいらないというかidMapから生成できる
        for(int i = 0; i < maxAMRLevel; ++i){
            segData[i] = new int[ numOfPatchesOnLevel[i] ];

            for(int j = 0; j < numOfPatchesOnLevel[i]; ++j){
                segData[i][j] = idMap[i][j];
            }
        }

        DBPutGroupelmap(file, "levelMap", maxAMRLevel, levelSegTypes, numOfPatchesOnLevel, 0, segData, 0, 0, 0);

        // -- make child map --
        int* patchSegTypes = new int[numOfPatches];
        for(int i = 0; i < numOfPatches; i++){
            patchSegTypes[i] = DB_BLOCKCENT;
        }

        // 1D child number array
        int* childOfPatches = new int[numOfPatches];
        // segData must have child IDs of each patch
        segData = new int*[numOfPatches];

        for(int i = 0; i < numOfPatches; i++){
            childOfPatches[i] = childMap[i].size();
            segData[i] = new int[ childOfPatches[i] ];

            for(int j = 0; j < childOfPatches[i]; ++j) {
                segData[i][j] = childMap[i][j];
            }
        }

        DBPutGroupelmap(file, "childMap", numOfPatches, patchSegTypes, childOfPatches, 0, segData, 0, 0, 0);

        // free segData
        for(int i = 0; i < maxAMRLevel; i++){
            delete [] segData[i];
        }
        delete [] segData;

        // make mrgTree
        /* Create an mrg tree for inclusion in the file */
        DBmrgtree* mrgTree = DBMakeMrgtree(DB_MULTIMESH, 0, 2, 0);
        // 2 is maximum number of children

        /* Add a region for AMR configuration */
        DBAddRegion(mrgTree, "amr_decomp", 0, 2, 0, 0, 0, 0, 0, 0);
        DBSetCwr(mrgTree, "amr_decomp");
        DBAddRegion(mrgTree, "levels", 0, maxAMRLevel, 0, 0, 0, 0, 0, 0);
        DBSetCwr(mrgTree, "levels");

        // levels
        // segTypes は使い回し
        {
            char* levelRegnNames[1];
            int* segIds = new int[maxAMRLevel];

            // これはセグメントへのID付けであって、データ自体へのID付けではない
            for (int i = 0; i < maxAMRLevel; ++i) {
                segIds[i] = i;
            }

            // printf style
            levelRegnNames[0] = const_cast<char*>("@level%d@n");
            DBAddRegionArray(mrgTree, maxAMRLevel, levelRegnNames, 0, "levelMap", 1, segIds, numOfPatchesOnLevel, levelSegTypes, 0);

            delete [] segIds;
        }
        DBSetCwr(mrgTree, "..");
        DBAddRegion(mrgTree, "patches", 0, numOfPatches, 0, 0, 0, 0, 0, 0);
        DBSetCwr(mrgTree, "patches");

        // patches
        {
            char* patchRegnNames[1];
            int* segIds = new int[numOfPatches];

            //! これはセグメントへのID付けであって、データ自体へのID付けではないが
            //! childOfPatchesがID順に格納されているとしているので、自動的に一致する
            for (int i = 0; i < numOfPatches; ++i) {
                segIds[i] = i;
            }
            patchRegnNames[0] = const_cast<char*>("@patch%d@n");
            DBAddRegionArray(mrgTree, numOfPatches, patchRegnNames, 0, "childMap", 1, segIds, childOfPatches, patchSegTypes, 0);

            delete [] segIds;
        }

        /*{
            DBoptlist* optList = DBMakeOptlist(10);
            char* mrgv_onames[5];
            mrgv_onames[0] = "lvlRatios";
            mrgv_onames[1] = "ijkExts";
            mrgv_onames[2] = "xyzExts";
            mrgv_onames[3] = "rank";
            mrgv_onames[4] = 0;

            DBAddOption(optList, DBOPT_MRGV_ONAMES, mrgv_onames);
            DBPutMrgtree(file, "mrgTree", "amr_mesh", mrgTree, optList);
            DBFreeMrgtree(mrgTree);
            DBFreeOptlist(optList);
        }*/

        const int numOfDim = 3;
        /* Output level refinement ratios as an mrg variable on the array of regions
           representing the levels */
        /*{
            char* compnames[3] = {"iRatio","jRatio","kRatio"};
            char* levelRegnNames[1];
            int* data[3];
            for (int i = 0; i < numOfDim; i++) {
                data[i] = new int[maxAMRLevel];
            }

            for (int i = 0; i < maxAMRLevel; i++)
            {
                if(i == 0) {
                    data[0][i] = 1;
                    data[1][i] = 1;
                    data[2][i] = 1;
                } else {
                    data[0][i] = 2;
                    data[1][i] = 2;
                    data[2][i] = 2;
                }
            }
            // data = { {1, 2, 2}, {1, 2, 2}, {1, 2, 2}};
            levelRegnNames[0] = const_cast<char*>("@level%d@n");
            DBPutMrgvar(file, "lvlRatios", "mrgTree", numOfDim, compnames, maxAMRLevel, levelRegnNames, DB_INT, data, 0);
            for (int i = 0; i < numOfDim; i++) {
                delete [] data[i];
            }
        }*/

        //! logical Extents
        //! Output logical extents of the patches as an mrg variable on the
        //!   array of regions representing the patches
        /*{
            char*  compnames[6] = {"iMin","iMax","jMin","jMax","kMin","kMax"};
            char* scompnames[6] = {"xMin","xMax","yMin","yMax","zMin","zMax"};
            char* patchRegnNames[1];
            int* data[6];
            float* rdata[1];
            float* sdata[6];
            patchRegnNames[0] = const_cast<char*>("@patch%d@n");

            for (int i = 0; i < 2 * numOfDim; i++) {
                data[i] = new int[numOfPatches];
                sdata[i] = new float[numOfPatches];
            }

            rdata[0] = new float[numOfPatches];

            root_grid->addExtent(data, sdata, rdata);

            DBPutMrgvar(file, "ijkExts", "mrgTree", 6, compnames, numOfPatches, patchRegnNames, DB_INT, data, 0);
            DBPutMrgvar(file, "xyzExts", "mrgTree", 6, scompnames, numOfPatches, patchRegnNames, DB_FLOAT, sdata, 0);
            DBPutMrgvar(file, "rank", "mrgTree", 1, 0, numOfPatches, patchRegnNames, DB_FLOAT, rdata, 0);

            for (int i = 0; i < 6; i++)
            {
                delete [] data[i];
                delete [] sdata[i];
            }
            delete [] rdata[0];
        }*/

        delete [] levelSegTypes;
        delete [] patchSegTypes;
        delete [] numOfPatchesOnLevel;
        delete [] childOfPatches;
    }

    void writeDomainGroupelMap(DBfile* file) {
        int numprocs = MPIw::Environment::numprocs;
        int* segTypes = new int[numprocs];
        for(int i = 0; i < numprocs; i++){
            segTypes[i] = DB_BLOCKCENT;
        }

        // とりあえず全部1メッシュとする
        int* numPatchesOnProcess = new int[numprocs];

        // 各プロセスの持つメッシュ数?
        int** segData = new int*[numprocs];
        for(int i = 0; i < numprocs; ++i){
            numPatchesOnProcess[i] = 1;
            segData[i] = new int[ numPatchesOnProcess[i] ];

            for(int j = 0; j < 1; ++j){
                segData[i][j] = 10000*MPIw::Environment::rank; // + idMap[i][j];
            }
        }

        DBPutGroupelmap(file, "domainMap", numprocs, segTypes, numPatchesOnProcess, 0, segData, 0, 0, 0);
    }

    void writeMRGTree(DBfile* file) {
        int maxGroupChildren = 1;
        int numprocs = MPIw::Environment::numprocs;
        DBmrgtree* mrgTree = DBMakeMrgtree(DB_MULTIMESH, 0, maxGroupChildren, 0);

        // Add a region for Multimesh Groupings
        int maxBlock = numprocs;
        DBAddRegion(mrgTree, "groupings", 0, maxBlock, 0, 0, 0, 0, 0, 0);
        DBSetCwr(mrgTree, "groupings");

        /*
        {
            char* regnNames[1];
            int* segIds = new int[numprocs];
            int* segTypes = new int[numprocs];

            // これはセグメントへのID付けであって、データ自体へのID付けではない
            for (int i = 0; i < numprocs; ++i) {
                segIds[i] = i;
                segTypes[i] = DB_BLOCKCENT;
            }

            // printf style
            regnNames[0] = const_cast<char*>("@grouping_block%04d@n");
            DBAddRegionArray(mrgTree, numprocs, regnNames, 0, "domainMap", 0, 0, 0, segTypes, 0);

            delete [] segIds;
        }
        */

        {
            DBoptlist* optList = DBMakeOptlist(1);
            const char* mrgv_onames[2];
            mrgv_onames[0] = "groupings";
            mrgv_onames[1] = 0;
            // mrgv_onames[1] = "ijkExts";
            // mrgv_onames[2] = "xyzExts";
            // mrgv_onames[3] = "rank";
            // mrgv_onames[4] = 0;

            DBAddOption(optList, DBOPT_MRGV_ONAMES, mrgv_onames);
            DBPutMrgtree(file, "mrgTree", "amr_mesh", mrgTree, optList);
            DBFreeMrgtree(mrgTree);
            DBFreeOptlist(optList);
        }
    }

    void writeMultimesh(DBfile* file, int total_blocknum, char** meshnames, char** varnames, std::string varlabel, float datatime) {
        // make options list
        DBoptlist* optList = DBMakeOptlist(2);
        int meshtype = DB_QUAD_RECT;
        DBAddOption(optList, DBOPT_MB_BLOCK_TYPE, &meshtype);
        DBAddOption(optList, DBOPT_TIME, &datatime);

        const char* blockname = "multimesh";
        DBPutMultimesh(file, blockname, total_blocknum, meshnames, NULL, optList);
        DBClearOptlist(optList);

        int vartype = DB_QUADVAR;
        const char* varblockname = varlabel.c_str();
        DBAddOption(optList, DBOPT_MB_BLOCK_TYPE, &vartype);
        DBAddOption(optList, DBOPT_TIME, &datatime);
        DBPutMultivar(file, varblockname, total_blocknum, varnames, NULL, optList);
        DBFreeOptlist(optList);
    }

    void writeBlock(DBfile* file, Grid* g, int rankInGroup, std::string dataTypeName){
        // dimension
        const int dim = 3;
        // names of the coordinates
        const char* coordnames[3] = {"x", "y", "z"};

        // make options list for mesh
        DBoptlist* optListMesh = DBMakeOptlist(1);
        // char* mrgTreeName = const_cast<char*>("mrgTree");
        // DBAddOption(optListMesh, DBOPT_MRGTREE_NAME, mrgTreeName);

        // make options list for var
        DBoptlist* optListVar = DBMakeOptlist(2);

        // set unit
        char* unit;
        if(dataTypeName == "potential") {
            unit = const_cast<char*>("V");
        } else if(dataTypeName == "rho") {
            unit = const_cast<char*>("/m^3");
        } else if (dataTypeName == "efield") {
            unit = const_cast<char*>("V/m");
        } else if (dataTypeName == "bfield") {
            unit = const_cast<char*>("T");
        } else {
            throw std::invalid_argument("[ERROR] Invalid dataTypeName was passed.");
        }
        DBAddOption(optListVar, DBOPT_UNITS, unit);

        int major_order = 1;
        DBAddOption(optListVar, DBOPT_MAJORORDER, &major_order); // column-major (Fortran) order

        g->putQuadMesh(file, dataTypeName, coordnames, rankInGroup, optListMesh, optListVar);

        // Free optList
        DBFreeOptlist(optListMesh);
        DBFreeOptlist(optListVar);
    }


    // Callback functions for PMPIO
    void* createFileCallback(const char* fname, const char* dname, void* udata){
        DBfile* file = DBCreate(fname, DB_CLOBBER, DB_LOCAL, NULL, DB_HDF5);
        DBMkDir(file, dname);
        DBSetDir(file, dname);
        return file;
    }

    void* openFileCallback(const char* fname, const char* dname, PMPIO_iomode_t iomode, void* udata){
        DBfile* file = DBOpen(fname, DB_UNKNOWN, DB_APPEND);
        DBMkDir(file, dname);
        DBSetDir(file, dname);
        return file;
    }

    void closeFileCallback(void* file, void* udata){
        DBClose(static_cast<DBfile*>(file));
    }

    void writeDataInParallel(Grid* g, int timestep, std::string dataTypeName) {
        //! @note: 事前に全てのプロセスの持つパッチ数などを云々しておく
        const float datatime = static_cast<float>(timestep * Environment::dt + 101);
        const int maxIOUnit = 4;
        int numfiles = (MPIw::Environment::numprocs <= maxIOUnit) ? MPIw::Environment::numprocs : maxIOUnit;
        int rank = MPIw::Environment::rank;

        PMPIO_baton_t* bat = PMPIO_Init(numfiles, PMPIO_WRITE, MPI_COMM_WORLD, 10000, createFileCallback, openFileCallback, closeFileCallback, NULL);
        int groupRank = PMPIO_GroupRank(bat, rank);
        int rankInGroup = PMPIO_RankInGroup(bat, rank);
        std::string filename = (format("%s_%04d_%04d.silo") % dataTypeName % groupRank % timestep).str();
        std::string blockname = (format("block%04d") % rankInGroup).str();
        DBfile* file = (DBfile*)(PMPIO_WaitForBaton(bat, ("data/" + filename).c_str(), blockname.c_str()));

        // 自分の持つroot_gridを再帰的にPutQuadmesh, PutQuadvarする
        writeBlock(file, g, rankInGroup, dataTypeName);

        if (rank == 0) {
            std::vector<int> patchesOnEachProcess(MPIw::Environment::numprocs);

            for(int process_num = 0; process_num < MPIw::Environment::numprocs; ++process_num) {
                // 各プロセスの持つパッチ数
                if(process_num == 0) {
                    patchesOnEachProcess[process_num] = g->getSumOfChild() + 1;
                } else {
                    patchesOnEachProcess[process_num] = 1;
                }
            }

            // 自分の持つroot_gridを統合したMultimeshを作成する
            int numAllPatches = std::accumulate(patchesOnEachProcess.begin(), patchesOnEachProcess.end(), 0);

            char** meshnames = new char*[numAllPatches];
            char** varnames = new char*[numAllPatches];

            //! 全パッチ数を保存する変数
            int index = 0;
            for(int process_num = 0; process_num < MPIw::Environment::numprocs; ++process_num) {
                for(int id = 0; id < patchesOnEachProcess[process_num]; ++id){
                    std::string tmpmeshname, tmpvarname;
                    if(process_num == 0) {
                        tmpmeshname = (format("/block%04d/mesh%04d%04d") % rankInGroup % process_num % id).str();
                        tmpvarname = (format("/block%04d/%s%04d%04d") % rankInGroup % dataTypeName % process_num % id).str();
                    } else {
                        int tmpGroupRank = PMPIO_GroupRank(bat, process_num);
                        int tmpRankInGroup = PMPIO_RankInGroup(bat, process_num);
                        std::string tmpfilename = (format("%s_%04d_%04d.silo") % dataTypeName % tmpGroupRank % timestep).str();

                        tmpmeshname = (format("%s:/block%04d/mesh%04d%04d") % tmpfilename % tmpRankInGroup % process_num % id).str();
                        tmpvarname = (format("%s:/block%04d/%s%04d%04d") % tmpfilename % tmpRankInGroup % dataTypeName % process_num % id).str();
                    }
                    meshnames[index] = new char[tmpmeshname.size() + 1];
                    std::strcpy(meshnames[index], tmpmeshname.c_str());

                    varnames[index] = new char[tmpvarname.size() + 1];
                    std::strcpy(varnames[index], tmpvarname.c_str());
                    ++index;
                }
            }

            writeMultimesh(file, numAllPatches, meshnames, varnames, dataTypeName, datatime);

            for(int i = 0; i < numAllPatches; ++i) {
                delete [] meshnames[i];
                delete [] varnames[i];
            }
            delete [] meshnames;
            delete [] varnames;
        }

        /*
        std::string mrgTreeRegionName;
        if (rank == 0) {
            // mrgTreeRegionName = (format("grouping_block%04d") % rank).str();
            mrgTreeRegionName = "groupings";
            writeDomainGroupelMap(file);
            writeMRGTree(file);
        } else {
            // mrgTreeRegionName = (format("potential_0000_0000.silo:grouping_block%04d") % rank).str();
            mrgTreeRegionName = "potential_0000_0000.silo:groupings";
        }
        */

        /*
        int maxAMRLevel = g->getMaxLevel();
        int* numPatchesOnLevel = g->getNumOfPatches();
        std::map<int, std::vector<int> > childMap = g->getChildMapOnRoot();
        std::vector< std::vector<int> > idMap = g->getIDMapOnRoot();
        writeGroupelMap(file, g, maxAMRLevel, g->getSumOfChild() + 1, numPatchesOnLevel, idMap, childMap);
        */
        PMPIO_HandOffBaton(bat, file);
        PMPIO_Finish(bat);
    }

    void plotEnergy(Grid* g, int timestep){
        const double datatime = timestep * Environment::dt;
        double particleEnergy = g->getParticleEnergy();
        double fieldEnergy = g->getFieldEnergy();
        double receivedParticleEnergy = MPIw::Environment::Comms["world"]->sum(particleEnergy, 0);
        double receivedFieldEnergy = MPIw::Environment::Comms["world"]->sum(fieldEnergy, 0);

        if(Environment::isRootNode) {
            std::string filename = "data/energy.txt";
            auto openmode = (timestep == 0) ? std::ios::out : std::ios::app;
            std::ofstream ofs(filename, openmode);

            if(timestep == 0) {
                ofs << "# " << format("%8s %10s %15s %15s %15s") % "timestep" % "time" % "Energy [J]" % "Particle [J]" % "Field [J]" << endl;
            }

            ofs << format("%10d %10.5f %15.7e %15.7e %15.7e") % timestep % datatime %
                Utils::Normalizer::unnormalizeEnergy(receivedParticleEnergy + receivedFieldEnergy) %
                Utils::Normalizer::unnormalizeEnergy(receivedParticleEnergy) %
                Utils::Normalizer::unnormalizeEnergy(receivedFieldEnergy) << endl;
        }
    }

    void print3DArray(const tdArray& data){
        const int nx = data.shape()[0];
        const int ny = data.shape()[1];
        const int nz = data.shape()[2];

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

    void outputParticlePositions(const ParticleArray& parray, std::string filename){
        std::ofstream ofs(filename);

        for(int id = 0; id < Environment::num_of_particle_types; ++id){

            ofs << "## " << Environment::ptype[id].getName() << endl;

            for(int i = 0; i < parray[id].size(); ++i){
                ofs << format("%9.4f %9.4f %9.4f") % parray[id][i].x % parray[id][i].y % parray[id][i].z;
                ofs << format("%13.4e %13.4e %13.4e") % parray[id][i].vx % parray[id][i].vy % parray[id][i].vz;
                ofs << endl;
            }

            ofs << endl << endl;
        }
    }
}
