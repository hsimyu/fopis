#ifndef __TDPIC_H_INCLUDED__
#define __TDPIC_H_INCLUDED__
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <memory>
#include <boost/format.hpp>
#include <picojson.h>

#define ARRAY_LENGTH(ARR) (sizeof(ARR) / sizeof((ARR)[0]))

// proto types
class Particle;
class ParticleType;

typedef double*** threeDArray;
typedef std::vector< std::vector<Particle> > ParticleArray;

struct Environment {
    public:
        int num_of_particle_types;
        double dx;
        double dt;
        int nx, ny, nz;
        int proc_x, proc_y, proc_z;
        int cell_x, cell_y, cell_z;
        std::string jobtype;
        int max_iteration;

        ParticleType* ptype;

        friend std::ostream& operator<<(std::ostream&, const Environment&);
        friend std::ostream& operator<<(std::ostream&, const Environment*);
        friend std::istream& operator<<(std::istream&, const Environment&);
};

class Field {
    private:
        threeDArray pPhi;
        threeDArray pRho;
        threeDArray pEx;
        threeDArray pEy;
        threeDArray pEz;
        threeDArray pBx;
        threeDArray pBy;
        threeDArray pBz;
    public:
        ~Field();

        threeDArray getPhi();
        void setPhi(threeDArray);

        threeDArray getRho();
        void setRho(threeDArray);

        threeDArray getEx();
        void setEx(threeDArray);
        threeDArray getEy();
        void setEy(threeDArray);
        threeDArray getEz();
        void setEz(threeDArray);

        threeDArray getBx();
        void setBx(threeDArray);
        threeDArray getBy();
        void setBy(threeDArray);
        threeDArray getBz();
        void setBz(threeDArray);
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

        friend std::ostream& operator<<(std::ostream&, const ParticleType&);
        friend std::istream& operator<<(std::istream&, const ParticleType&);
};

class Grid {
    private:
        Grid* parent;
        std::vector< std::unique_ptr<Grid*> > children;

        //! 親のどの座標にくっついているか
        //! @{
        int base_x;
        int base_y;
        int base_z;
        //! @}

        int nx, ny, nz;
        int level;

        Field* field;

    public:
        Grid(const Environment*);
        ~Grid();

        void setBaseX(int);
        void setBaseY(int);
        void setBaseZ(int);

        int  getBaseX();
        int  getBaseY();
        int  getBaseZ();

        void setLevel(int);
        int getLevel();

        void setParent(Grid*);
        Grid* getParent();

        void setField(Field*);
        Field* getField();

        void addChild(Grid*);
        std::vector< std::unique_ptr<Grid*> > getChildren();
        void removeChild();

        //! Particleを格納したstd::vectorを、粒子種ごとに保持したstd::vector
        ParticleArray particles;

        // update fields
        void updateRho(const Environment*);
};

namespace Initializer {
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

    double*** create3DArray(const int nx, const int ny, const int nz);

    void delete3DArray(double*** x);
}

namespace IO {
    void print3DArray(double***, const int, const int, const int);
    void outputParticlePositions(const Environment*, const ParticleArray&, std::string filename = "particlePositions.csv");
}
#endif
