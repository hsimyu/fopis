#ifndef __2DPIC_H_INCLUDED__
#define __2DPIC_H_INCLUDED__
#include <boost/multi_array.hpp>

#define ARRAY_LENGTH(ARR) (sizeof(ARR) / sizeof((ARR)[0]))

typedef boost::multi_array<double, 3> threeD_array;

class FieldPointers {
    private:
        threeD_array* pPhi;
        threeD_array* pRho;
    public:
        threeD_array* getPhi();
        void setPhi(threeD_array*);

        threeD_array* getRho();
        void setRho(threeD_array*);
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
