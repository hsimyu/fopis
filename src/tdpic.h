#ifndef __TDPIC_H_INCLUDED__
#define __TDPIC_H_INCLUDED__
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
        // 8byte * 6 = 48 byte
        double x, y, z;
        double vx, vy, vz;
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

class ParticleInfo {
    private:
        int sizeOfSuperParticle;
    public:
        ParticleInfo();
        void setParticleSize(int);
        int getParticleSize();
};

namespace Initializer {
    double getSizeOfSuperParticle(int, double, double);
    void setParticleInfo(ParticleInfo*);
    void setFieldPointers(FieldPointers*, int, int, int);
}

namespace Utils {
    void print3DArray(threeD_array*);
}
#endif
