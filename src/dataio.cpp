#include <tdpic.h>
#include <numeric>
#include <fstream>

#include <mpi.h>
#include <pmpio.h>

namespace IO {
    void writeMRGTree(DBfile* file) {
        int maxGroupChildren = 1;
        DBmrgtree* mrgTree = DBMakeMrgtree(DB_MULTIMESH, 0, maxGroupChildren, 0);

        // Add a region for Multimesh Groupings
        int maxBlock = 1;
        DBAddRegion(mrgTree, "groupings", 0, maxBlock, 0, 0, 0, 0, 0, 0);
        DBSetCwr(mrgTree, "groupings");

        int maxChildren = 1;
        DBAddRegion(mrgTree, "grouping_block0000", 0, maxChildren, 0, 0, 0, 0, 0, 0);
        DBSetCwr(mrgTree, "grouping_block0000");

        {
            DBoptlist* optList = DBMakeOptlist(10);
            // char* mrgv_onames[1];
            // mrgv_onames[0] = "grouping_block0000";
            // mrgv_onames[1] = "ijkExts";
            // mrgv_onames[2] = "xyzExts";
            // mrgv_onames[3] = "rank";
            // mrgv_onames[4] = 0;

            // DBAddOption(optList, DBOPT_MRGV_ONAMES, mrgv_onames);
            DBPutMrgtree(file, "mrgTree", "amr_mesh", mrgTree, NULL);
            DBFreeMrgtree(mrgTree);
            DBFreeOptlist(optList);
        }
    }

    void writeMultimesh(DBfile* file, int total_blocknum, char** meshnames, char** varnames, int rank) {
        // make options list
        DBoptlist* optList = DBMakeOptlist(2);
        int meshtype = DB_QUAD_RECT;
        // multimesh を region nameと紐付ける
        // char* pname = const_cast<char*>("grouping_block0000");
        DBAddOption(optList, DBOPT_MB_BLOCK_TYPE, &meshtype);
        // DBAddOption(optList, DBOPT_REGION_PNAMES, &pname);
        std::string tmpstring = (format("multiblock%04d") % rank).str();
        char* blockname = const_cast<char*>( tmpstring.c_str() );
        DBPutMultimesh(file, blockname, total_blocknum, meshnames, NULL, optList);
        DBClearOptlist(optList);

        int vartype = DB_QUADVAR;
        tmpstring = (format("var%04d") % rank).str();
        char* varblockname = const_cast<char*>( tmpstring.c_str() );
        DBAddOption(optList, DBOPT_MB_BLOCK_TYPE, &vartype);
        DBPutMultivar(file, varblockname, total_blocknum, varnames, NULL, optList);
        DBFreeOptlist(optList);
    }

    void writeBlock(DBfile* file, Grid* g, std::string dataTypeName, int rankInGroup){
        // dimension
        const int dim = 3;

        // names of the coordinates
        char* coordnames[3];
        coordnames[0] = const_cast<char*>("x");
        coordnames[1] = const_cast<char*>("y");
        coordnames[2] = const_cast<char*>("z");

        // names of the variables
        char* varnames[1];
        varnames[0] = const_cast<char*>(dataTypeName.c_str());

        // make options list for mesh
        // DBoptlist* optListMesh = DBMakeOptlist(1);

        // make options list for var
        DBoptlist* optListVar = DBMakeOptlist(2);
        char* unit = const_cast<char*>("V");
        DBAddOption(optListVar, DBOPT_UNITS, unit);
        int major_order = 1;
        DBAddOption(optListVar, DBOPT_MAJORORDER, &major_order); // column-major (Fortran) order

        g->putQuadMesh(file, coordnames, varnames, optListVar, rankInGroup);

        // Free optList
        // DBFreeOptlist(optListMesh);
        DBFreeOptlist(optListVar);
    }

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

        {
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
        }

        const int numOfDim = 3;
        /* Output level refinement ratios as an mrg variable on the array of regions
           representing the levels */
        {
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
        }

        //! logical Extents
        //! Output logical extents of the patches as an mrg variable on the
        //!   array of regions representing the patches
        {
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
        }

        delete [] levelSegTypes;
        delete [] patchSegTypes;
        delete [] numOfPatchesOnLevel;
        delete [] childOfPatches;
    }

    // Callback functions for PMPIO
    void* createFileCallback(const char* fname, const char* dname, void* udata){
        DBfile* file = DBCreate(fname, DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);
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
        DBClose((DBfile*)file);
    }

    void writeDataInParallel(Grid* g, int timestep, std::string dataTypeName) {
        //! @note: 事前に全てのプロセスの持つパッチ数などを云々しておく
        const int maxIOUnit = 4;
        int numfiles = (MPIw::Environment::numprocs <= maxIOUnit) ? MPIw::Environment::numprocs : maxIOUnit;

        PMPIO_baton_t* bat = PMPIO_Init(numfiles, PMPIO_WRITE, MPI_COMM_WORLD, 10000, createFileCallback, openFileCallback, closeFileCallback, NULL);
        int groupRank = PMPIO_GroupRank(bat, MPIw::Environment::rank);
        int rankInGroup = PMPIO_RankInGroup(bat, MPIw::Environment::rank);
        std::string filename = (format("data/%s_%04d_%04d.silo") % dataTypeName % groupRank % timestep).str();
        std::string blockname = (format("block%04d") % rankInGroup).str();
        DBfile* file = (DBfile*)(PMPIO_WaitForBaton(bat, filename.c_str(), blockname.c_str()));

        // 自分の持つroot_gridを統合したMultimeshを作成する
        // int total_blocknum = MPI::Environment::numprocs;
        // std::vector<int> patchNumOfEachProcess;
        // patchNumOfEachProcess.resize(total_blocknum);
        //
        // patchNumOfEachProcess[0] = g->getSumOfChild() + 1;
        // patchNumOfEachProcess[1] = 1;
        //
        // int numAllPatches = std::accumulate(patchNumOfEachProcess.begin(), patchNumOfEachProcess.end(), 0);
        int numAllPatches = g->getSumOfChild() + 1;
        /*char** meshnames = new char*[numAllPatches];
        char** varnames = new char*[numAllPatches];

        //! 全パッチ数を保存する変数
        for(int i = 0; i < numAllPatches; ++i) {
            std::string tmpstring;
            tmpstring = (format("/block%04d/mesh%04d%04d") % rankInGroup % MPIw::Environment::rank % i).str();

            meshnames[i] = new char[tmpstring.size() + 1];
            std::strcpy(meshnames[i], tmpstring.c_str());

            tmpstring = (format("/block%04d/%s%04d%04d") % rankInGroup % dataTypeName % MPIw::Environment::rank % i).str();
            varnames[i] = new char[tmpstring.size() + 1];
            std::strcpy(varnames[i], tmpstring.c_str());

            cout << MPIw::Environment::rankStr() << "meshnames[" << i << "] = " << meshnames[i] << endl;
            cout << MPIw::Environment::rankStr() << "varnames[" << i << "] = " << varnames[i] << endl;
        }*/

        // writeMultimesh(file, numAllPatches, meshnames, varnames, rankInGroup);
        // writeMRGTree(file);

        /*
        int maxAMRLevel = g->getMaxLevel();
        int* numPatchesOnLevel = g->getNumOfPatches();
        std::map<int, std::vector<int> > childMap = g->getChildMapOnRoot();
        std::vector< std::vector<int> > idMap = g->getIDMapOnRoot();
        writeGroupelMap(file, g, maxAMRLevel, g->getSumOfChild() + 1, numPatchesOnLevel, idMap, childMap);
        */

        writeBlock(file, g, dataTypeName, rankInGroup);

        PMPIO_HandOffBaton(bat, file);
        PMPIO_Finish(bat);

        // for(int i = 0; i < numAllPatches; ++i) {
        //     delete [] meshnames[i];
        //     delete [] varnames[i];
        // }
        // delete [] meshnames;
        // delete [] varnames;
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
