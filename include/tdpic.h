#ifndef __TDPIC_H_INCLUDED__
#define __TDPIC_H_INCLUDED__
#include <iostream>
#include <stdexcept>
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include <memory>
#include <map>
#include <mpi.h>
#include <boost/format.hpp>
#include <boost/multi_array.hpp>
#include <silo.h>
#include <picojson.h>
#include <mkl.h> // Intel Math Kernel Library

#define ARRAY_LENGTH(ARR) (sizeof(ARR) / sizeof((ARR)[0]))

// Prototype宣言のみ書く
class Particle;
class ParticleType;

using std::cout;
using std::endl;
using boost::format;

typedef boost::multi_array<double, 3> tdArray;
typedef std::vector< std::vector<Particle> > ParticleArray;

//! constants
const double e = 1.60217733e-19;
const double eps0 = 8.854187817e-12;
const double me = 9.1093818872e-31;
const double c = 2.99792458e8;

struct Poisson {
    public:
        ~Poisson();

        MKL_INT nx;
        MKL_INT ny;
        MKL_INT nz;

        MKL_INT* ipar;
        double* dpar;
        MKL_INT stat;
        DFTI_DESCRIPTOR_HANDLE* xhandle;
        DFTI_DESCRIPTOR_HANDLE* yhandle;
        double* b_lx;
        double* b_ly;
        double* b_lz;

        double* rho1D;
};

struct Environment {
    public:
        static int num_of_particle_types;
        static double dx;
        static double dt;
        static int nx, ny, nz;
        static int proc_x, proc_y, proc_z;
        static int cell_x, cell_y, cell_z;
        static int max_iteration;
        static bool isRootNode;

        static bool onLowXedge, onHighXedge;
        static bool onLowYedge, onHighYedge;
        static bool onLowZedge, onHighZedge;

        static std::string jobtype;
        static std::string solver_type;
        static std::string boundary;
        static std::string dimension;

        static ParticleType* ptype;

        //! 中身はMPI::Environment::getRankStr()と同様
        static std::string rankStr(void);
        static void printInfo(void);
};

class Field {
    private:
        tdArray rho;
        tdArray phi;

        tdArray ex;
        tdArray ey;
        tdArray ez;
        tdArray exref;
        tdArray eyref;
        tdArray ezref;

        tdArray bx;
        tdArray by;
        tdArray bz;

        Poisson* psn = nullptr;
    public:
        Field();
        ~Field();

        // scalars
        tdArray& getPhi();
        void setPhi(tdArray&);
        tdArray& getRho();
        void setRho(tdArray&);
        tdArray& getScalar(std::string);

        // vectors
        tdArray& getEx();
        void setEx(tdArray&);
        tdArray& getEy();
        void setEy(tdArray&);
        tdArray& getEz();
        void setEz(tdArray&);

        // reference efield on nodes
        tdArray& getExRef();
        tdArray& getEyRef();
        tdArray& getEzRef();
        void setExRef(tdArray&);
        void setEyRef(tdArray&);
        void setEzRef(tdArray&);

        tdArray& getBx();
        void setBx(tdArray&);
        tdArray& getBy();
        void setBy(tdArray&);
        tdArray& getBz();
        void setBz(tdArray&);

        void initializePoisson(const int, const int, const int);
        void updateEfield(const int, const int, const int);
        void updateBfield(const int, const int, const int);

        //! solver の実体
        void solvePoissonMKL(const int, const int, const int);
        void solvePoissonSOR(const int, const double);
        void solvePoisson(const int, const double);
};

class Velocity {
    private:
        double vx, vy, vz;
    public:
        Velocity();
        ~Velocity();
        void set(double, double, double);
        double getVX(void) const;
        double getVY(void) const;
        double getVZ(void) const;
};

class Position {
    private:
        void updateDelta(void);

    public:
        int i, j, k;
        double x, y, z;
        double dx1, dy1, dz1, dx2, dy2, dz2;
        Position(const double, const double, const double);
        Position(const int, const int, const int);
        Position(const Particle&);
        ~Position();
        void setXYZ(const double, const double, const double);
        void setIJK(const int, const int, const int);
};

class Particle {
    public:
        // MPIで送るためにoffsetofを使う必要があり、
        // publicにする必要がある
        int typeId;
        int isValid;
        double x, y, z;
        double vx, vy, vz;

        Particle();
        Particle(Particle const&);
        ~Particle();

        void setPosition(Position const&);
        void setPosition(const double, const double, const double);

        void setVelocity(Velocity const&);
        void setVelocity(const double, const double, const double);

        double getEnergy(void) const;
        double getSquaredMagnitudeOfVelocity(void) const;

        void updatePosition(void);
        Particle& operator=(Particle const&);
        friend std::ostream& operator<<(std::ostream&, Particle const&);
};

class ParticleType {
    private:
        int id;
        std::string name;

        // plasma info
        std::string type;
        double charge;
        double mass;
        double density;
        double temperature;
        int size;

        // for computation
        int particle_per_cell;
        int totalNumber;
    public:
        ParticleType();

        void setId(int);
        void setCharge(double);
        void setMass(double);
        void setDensity(double);
        void setTemperature(double);
        void setSize(int);
        void setName(std::string);
        void setType(std::string);
        void setTotalNumber(int);
        void setPcell(int);

        int getId() const;
        std::string getType() const;
        std::string getName() const;
        double getCharge() const;
        double getMass() const;
        double getDensity() const;
        double getTemperature() const;
        double getTrueTemperature() const;
        int getSize() const;
        int getTotalNumber() const;
        int getPcell() const;

        int calcSize(void);
        int calcTotalNumber(void);
        double calcDeviation(void) const;

        friend std::ostream& operator<<(std::ostream&, const ParticleType&);
        friend std::istream& operator<<(std::istream&, const ParticleType&);
};

//! @class Grid
class Grid {
    private:
        Grid* parent;
        std::vector<Grid*> children;

        //! 親のどの座標にくっついているか
        //! @{
        int from_ix;
        int from_iy;
        int from_iz;
        int to_ix;
        int to_iy;
        int to_iz;
        double base_x;
        double base_y;
        double base_z;
        //! @}

        unsigned int id;
        int nx, ny, nz;
        int level;
        int sumTotalNumOfChildGrids;
        double dx;

        Field* field;

        // Class Unique ID
        static unsigned int nextID;

    public:
        Grid(void);
        Grid(Grid*, const int, const int, const int, const int, const int, const int);

        ~Grid();

        unsigned int getID(void) const;
        unsigned int getNextID(void);
        // IDへのsetterは提供しない、代わりにresetterを用意
        static void resetNextID(void);

        void setFromIX(int);
        void setFromIY(int);
        void setFromIZ(int);
        int  getFromIX() const;
        int  getFromIY() const;
        int  getFromIZ() const;
        void setToIX(int);
        void setToIY(int);
        void setToIZ(int);
        int  getToIX() const;
        int  getToIY() const;
        int  getToIZ() const;

        void setBaseX(double);
        void setBaseY(double);
        void setBaseZ(double);

        double getBaseX() const;
        double getBaseY() const;
        double getBaseZ() const;

        void setNX(int);
        void setNY(int);
        void setNZ(int);

        int  getNX() const;
        int  getNY() const;
        int  getNZ() const;

        void setLevel(int);
        int getLevel() const;

        void setDX(double);
        double getDX() const;

        void setParent(Grid*);
        Grid* getParent();

        void setField(Field*);
        Field* getField();

        // Field 初期化
        void initializeField(void);

        // 親子でのScalarやりとり用
        void copyScalarToChildren(std::string);
        void copyScalarToParent(std::string);

        // 子供管理メソッド
        void makeChild(const int, const int, const int, const int, const int, const int);
        void addChild(Grid*);
        std::vector<Grid*>& getChildren();
        int getChildrenLength(void) const;
        void removeChild(const int);
        void checkGridValidness(void);

        // child gridの総数を保存するためのメソッド
        int getSumOfChild(void) const;
        void setSumOfChild(const int);
        void incrementSumOfChild(void);
        void decrementSumOfChild(void);

        //! Particleを格納したstd::vectorを、粒子種ごとに保持したstd::vector
        ParticleArray particles;

        // update fields
        void updateRho(void);
        void solvePoisson(void);
        void updateEfield(void);
        void updateBfield(void);

        // update particles
        void updateParticleVelocity(void);
        void updateParticlePosition(void);

        //! エネルギーを計算する
        double getParticleEnergy(void);
        double getFieldEnergy(void);

        // create mesh nodes array
        float** getMeshNodes(int);
        int* getChildOfPatches(void);
        int* getNumOfPatches(void);
        int* getChildIdMap(void);
        int getMaxLevel(void);
        int getMaxChildLevel(void);

        // {id: [子のIDを格納したvector]}のmapを作成する
        std::map<int, std::vector<int> > getChildMapOnRoot(void);
        void addChildrenIDToMap(std::map<int, std::vector<int> >&);
        std::vector<int> getChildrenIDs(void);

        // [level: [IDを格納したvector]]のvectorを作成する
        std::vector< std::vector<int> > getIDMapOnRoot(void);
        void addIDToVector(std::vector< std::vector<int> >&);

        // Extentを指定されたポインタに格納する (Silo MRG Tree出力用)
        void addExtent(int* logicalExtent[6], float* spatialExtent[6], float* rank[1]);

        // QuadMeshとVarをDBfileに突っ込む
        void putQuadMesh(DBfile*, std::string, const char* coordnames[3], int, DBoptlist*, DBoptlist*);

        // std out
        friend std::ostream& operator<<(std::ostream&, Grid*);

        //! inline functions
        void checkXBoundary(std::vector< std::vector<Particle> >& pbuff, Particle& p, const double slx) {
            if(p.x < 0.0) {
                if(!Environment::onLowXedge) pbuff[0].push_back(p);
                p.isValid = 0;
            } else if (p.x > slx) {
                if(!Environment::onHighXedge) pbuff[1].push_back(p);
                p.isValid = 0;
            }
        }

        void checkYBoundary(std::vector< std::vector<Particle> >& pbuff, Particle& p, const double sly) {
            if(p.y < 0.0) {
                if(!Environment::onLowYedge) pbuff[2].push_back(p);
                p.isValid = 0;
            } else if (p.y > sly) {
                if(!Environment::onHighYedge) pbuff[3].push_back(p);
                p.isValid = 0;
            }
        }

        void checkZBoundary(std::vector< std::vector<Particle> >& pbuff, Particle& p, const double slz) {
            if(p.z < 0.0) {
                if(!Environment::onLowZedge) pbuff[4].push_back(p);
                p.isValid = 0;
            } else if (p.z > slz) {
                if(!Environment::onHighZedge) pbuff[5].push_back(p);
                p.isValid = 0;
            }
        }
};

namespace Initializer {
    void initTDPIC(Grid*&);
    void setTestEnvirontment(void);
    void setMPIInfoToEnv(void);
    void loadEnvironment(picojson::object&);
    void loadParticleType(picojson::object&);
    Grid* initializeGrid(void);
}

namespace Utils {
    void printTotalMemory(const ParticleType&);
    void printParticleMemory(const ParticleType&);
    std::string readFile(const std::string&);
    picojson::value::object readJSONFile(const std::string&);

    float* getTrueCells(const tdArray&);
    float* getTrueEdges(const tdArray&, const int);
    float* getTrueFaces(const tdArray&, const int);
    void convert1Dto3Darray(double*, const int, const int, const int, tdArray&);
    void clearBoundaryValues(tdArray&, const int, const int, const int);

    void createDir(std::string);

    class Normalizer {
        protected:
            // static class
            Normalizer();
            ~Normalizer();

        public:
            // static member
            static double x_unit;
            static double t_unit;
            static double m_unit;
            static double e_unit;

            static double normalizeLength(const double);
            static double unnormalizeLength(const double);

            static double normalizeVelocity(const double);
            static double unnormalizeVelocity(const double);

            static double normalizeTime(const double);
            static double unnormalizeTime(const double);

            static double normalizeCharge(const double);
            static double unnormalizeCharge(const double);

            static double normalizeMass(const double);
            static double unnormalizeMass(const double);

            static double normalizeEnergy(const double);
            static double unnormalizeEnergy(const double);

            static double normalizeEpsilon(const double);
            static double unnormalizeEpsilon(const double);
    };
}

namespace IO {
    void plotEnergy(Grid*, int);
    void writeDataInParallel(Grid*, int, std::string);
    void print3DArray(const tdArray&);
    void outputParticlePositions(const ParticleArray&, std::string filename = "data/particlePositions.csv");
}

namespace MPIw {
    class Communicator {
        private:
            MPI_Comm comm;

        public:
            Communicator();
            Communicator(MPI_Comm);
            //デストラクタは特に必要ない

            MPI_Comm getComm(void);
            void setComm(MPI_Comm);

            // -- communicate methods --
            void barrier();

            // summation
            double sum(double, const int);
            int sum(int, const int);

            // send
            void send(Particle const&, const int);
            void recv(Particle&, const int);
            void sendVector(std::vector<Particle> const&, const int);
            void recvVector(std::vector<Particle>&, const int);
            void sendRecvVector(std::vector<Particle> const&, std::vector<Particle>&, const int, const int);
    };

    MPI_Datatype registerParticleType(void);
    void deregisterMpiType(MPI_Datatype);

    //! MPI ランク保持用のクラス
    class Environment {
        private:
            static void finalize(void);

        public:
            Environment(int, char**);
            ~Environment();

            //! MPIのランクとプロセス数はstaticに持つ
            //! int rank = MPI::Environment::rank; でアクセスする
            static int rank;
            static int numprocs;
            static int xrank, yrank, zrank;
            static int adj[6];

            //! 通信用の型
            static MPI_Datatype MPI_PARTICLE;

            //! コミュニケータのリスト
            static std::map<std::string, Communicator*> Comms;

            //! Utility member functions
            static std::string rankStr(void);
            static void exitWithFinalize(int);
            static void sendRecvParticlesX(std::vector< std::vector<Particle> > const&, std::vector< std::vector<Particle> >&);
            static void sendRecvParticlesY(std::vector< std::vector<Particle> > const&, std::vector< std::vector<Particle> >&);
            static void sendRecvParticlesZ(std::vector< std::vector<Particle> > const&, std::vector< std::vector<Particle> >&);
            static void sendRecvParticles(std::vector< std::vector<Particle> > const&, std::vector< std::vector<Particle> >&, const int, const int, std::string);

    };

};
#endif
