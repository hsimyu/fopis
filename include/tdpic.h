#ifndef __TDPIC_H_INCLUDED__
#define __TDPIC_H_INCLUDED__
#include <iostream>
#include <string>
#include <math.h>
#include <memory>
#include <boost/multi_array.hpp>
#include <picojson.h>
#include <vector>

#define ARRAY_LENGTH(ARR) (sizeof(ARR) / sizeof((ARR)[0]))

typedef boost::multi_array<double, 3> threeD_array;

struct Environment {
    public:
        int particle_types;
        double dx;
        double dt;
        int nx, ny, nz;
        int proc_x, proc_y, proc_z;
        int cell_x, cell_y, cell_z;
        std::string jobtype;
        int max_iteration;

        friend std::ostream& operator<<(std::ostream&, const Environment&);
        friend std::ostream& operator<<(std::ostream&, const Environment*);
        friend std::istream& operator<<(std::istream&, const Environment&);
};

class Field {
    private:
        threeD_array* pPhi;
        threeD_array* pRho;
        threeD_array* pEx;
        threeD_array* pEy;
        threeD_array* pEz;
        threeD_array* pBx;
        threeD_array* pBy;
        threeD_array* pBz;
    public:
        ~Field();

        threeD_array* getPhi();
        void setPhi(threeD_array*);

        threeD_array* getRho();
        void setRho(threeD_array*);

        threeD_array* getEx();
        void setEx(threeD_array*);
        threeD_array* getEy();
        void setEy(threeD_array*);
        threeD_array* getEz();
        void setEz(threeD_array*);

        threeD_array* getBx();
        void setBx(threeD_array*);
        threeD_array* getBy();
        void setBy(threeD_array*);
        threeD_array* getBz();
        void setBz(threeD_array*);
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

        //! 粒子種ごとのParticle配列を格納したstd::vectorへのunique_ptrを保持する
        // std::vector< Particle[] > particles;
    public:
        //! vector< unique_ptr<> > はpublicにすべき？
        std::vector< std::unique_ptr<Particle[]> > particles;
        Grid(const Environment*, const ParticleType*);
        ~Grid();

        //! 親のX座標を設定します.
        void setBaseX(int);
        //! 親のY座標を設定します.
        void setBaseY(int);
        //! 親のZ座標を設定します.
        void setBaseZ(int);

        int  getBaseX();
        int  getBaseY();
        int  getBaseZ();

        void setLevel(int);
        void setParent(Grid*);
        void addChild(Grid*);
        int getLevel();
        Grid* getParent();
        std::vector< std::unique_ptr<Grid*> > getChildren();

        // void removeChild();
        // std::unique_ptr<Particle[]> particles(new Particle[ptype.getTotalNumber()]);
};

namespace Initializer {
    double getSizeOfSuperParticle(int, double, double);
    Environment* loadEnvironment(picojson::object&);
    ParticleType* loadParticleType(picojson::object&, Environment*);
    Field* initializeField(const Environment*);
    Grid* initializeGrid(const Environment*, const ParticleType*);
}

namespace Utils {
    void print3DArray(threeD_array*);
    void printTotalMemory(const ParticleType&);
    void printParticleMemory(const ParticleType&);
    std::string readFile(const std::string&);
    picojson::value::object readJSONFile(const std::string&);
}
#endif
