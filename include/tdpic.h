#ifndef __TDPIC_H_INCLUDED__
#define __TDPIC_H_INCLUDED__
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <memory>
#include <boost/format.hpp>
#include <boost/multi_array.hpp>
#include <picojson.h>
#include <mkl.h> // Intel Math Kernel Library
#include <mpi_wrapper.hpp>

#define ARRAY_LENGTH(ARR) (sizeof(ARR) / sizeof((ARR)[0]))

// proto types
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
        int num_of_particle_types;
        double dx;
        double dt;
        int nx, ny, nz;
        int proc_x, proc_y, proc_z;
        int cell_x, cell_y, cell_z;
        int max_iteration;
        bool isRootNode;

        bool onLowXedge, onHighXedge;
        bool onLowYedge, onHighYedge;
        bool onLowZedge, onHighZedge;

        // MPI info
        int numprocs;
        int rank;
        int xrank, yrank, zrank;

        std::string jobtype;
        std::string solver_type;
        std::string boundary;
        std::string dimension;

        ParticleType* ptype;

        //! 中身はMPI::Environment::getRankStr()と同様
        std::string rankStr(void) const;

        friend std::ostream& operator<<(std::ostream&, const Environment&);
        friend std::ostream& operator<<(std::ostream&, const Environment*);
        friend std::istream& operator<<(std::istream&, const Environment&);
};

class Field {
    private:
        tdArray rho;
        tdArray phi;
        tdArray ex;
        tdArray ey;
        tdArray ez;
        tdArray bx;
        tdArray by;
        tdArray bz;

        Poisson* psn = nullptr;
    public:
        Field();
        ~Field();

        tdArray& getPhi();
        void setPhi(tdArray&);

        tdArray& getRho();
        void setRho(tdArray&);

        tdArray& getEx();
        void setEx(tdArray&);
        tdArray& getEy();
        void setEy(tdArray&);
        tdArray& getEz();
        void setEz(tdArray&);

        tdArray& getBx();
        void setBx(tdArray&);
        tdArray& getBy();
        void setBy(tdArray&);
        tdArray& getBz();
        void setBz(tdArray&);

        void initializePoisson(const Environment*);
        void solvePoisson(const Environment*);
        void updateEfield(const Environment*);
        void updateBfield(const Environment*);
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
        double x, y, z;
    public:
        Position();
        ~Position();
        void set(double, double, double);
        double getX(void) const;
        double getY(void) const;
        double getZ(void) const;
};

class Particle {
    private:
        // 8byte * 6
        // + 2byte   = 50 byte
        double x, y, z;
        double vx, vy, vz;
        unsigned short int id;
    public:
        Particle();
        ~Particle();

        void setPosition(const Position&);
        void setPosition(double, double, double);

        void setVelocity(const Velocity&);
        void setVelocity(double, double, double);
        Position getPosition();
        Velocity getVelocity();

        double getX(void) const;
        double getY(void) const;
        double getZ(void) const;
        double getVX(void) const;
        double getVY(void) const;
        double getVZ(void) const;
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
        int getSize() const;
        int getTotalNumber() const;
        int getPcell() const;

        int calcSize(const Environment*);
        int calcTotalNumber(const Environment*);
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
        int base_ix;
        int base_iy;
        int base_iz;
        double base_x;
        double base_y;
        double base_z;
        //! @}

        int nx, ny, nz;
        int level;
        int sumTotalNumOfChildGrids;
        double dx;

        Field* field;
    public:
        Grid(const Environment*);
        ~Grid();

        Grid(Grid*, const int, const int, const int, const int, const int, const int);

        void setBaseIX(int);
        void setBaseIY(int);
        void setBaseIZ(int);

        int  getBaseIX() const;
        int  getBaseIY() const;
        int  getBaseIZ() const;

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
        void updateRho(const Environment*);

        // create mesh nodes array
        float** getMeshNodes(int);
        int* getChildrenNumbers();

        // std out
        friend std::ostream& operator<<(std::ostream&, Grid*);
};

namespace Initializer {
    void setMPIInfoToEnv(Environment*);
    double getSizeOfSuperParticle(int, double, double);
    Environment* loadEnvironment(picojson::object&);
    ParticleType* loadParticleType(picojson::object&, Environment*);
    void initializeRootField(const Environment*, Grid*);
    Grid* initializeGrid(const Environment*);
}

namespace Utils {
    void printTotalMemory(const ParticleType&);
    void printParticleMemory(const ParticleType&);
    std::string readFile(const std::string&);
    picojson::value::object readJSONFile(const std::string&);

    double* getTrueCells(const tdArray&);
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
            static double e_unit;

            static double normalizeLength(const double);
            static double unnormalizeLength(const double);

            static double normalizeVelocity(const double);
            static double unnormalizeVelocity(const double);

            static double normalizeTime(const double);
            static double unnormalizeTime(const double);

            static double normalizeCharge(const double);
            static double unnormalizeCharge(const double);

    };
}

namespace IO {
    void writeDataInParallel(Grid*, int, std::string);
    void print3DArray(const tdArray&, const int, const int, const int);
    void outputParticlePositions(const Environment*, const ParticleArray&, std::string filename = "data/particlePositions.csv");
}
#endif
