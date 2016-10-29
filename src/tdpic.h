#ifndef __TDPIC_H_INCLUDED__
#define __TDPIC_H_INCLUDED__
#include <boost/multi_array.hpp>
#include <iostream>
#include <string>
#include <math.h>
#include <memory>
#include <boost/multi_array.hpp>

#define ARRAY_LENGTH(ARR) (sizeof(ARR) / sizeof((ARR)[0]))

typedef boost::multi_array<double, 3> threeD_array;

class FieldPointers {
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
        void setTotalNumber(int);
        void setPcell(int);

        int getId() const;
        std::string getName() const;
        double getCharge() const;
        double getMass() const;
        double getDensity() const;
        double getTemperature() const;
        int getSize() const;
        int getTotalNumber() const;
        int getPcell() const;

        int calcTotalNumber(int, int, int, int);

        friend std::ostream& operator<<(std::ostream&, const ParticleType&);
        friend std::istream& operator<<(std::istream&, const ParticleType&);
};

namespace Initializer {
    double getSizeOfSuperParticle(int, double, double);
    void setFieldPointers(FieldPointers*, int, int, int);
}

namespace Utils {
    void print3DArray(threeD_array*);
    void printTotalMemory(const ParticleType&);
    void printParticleMemory(const ParticleType&);
}
#endif
