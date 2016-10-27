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

class ParticleInfo {
    public:
        ParticleInfo();
        int sizeOfSuperParticle;
};

namespace Initializer {
    double getSizeOfSuperParticle(int, double, double);
    void setFieldPointers(FieldPointers*, int, int, int);
}

namespace Utils {
    void print3DArray(threeD_array*);
}
#endif
