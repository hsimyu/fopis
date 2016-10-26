#ifndef __2DPIC_H_INCLUDED__
#define __2DPIC_H_INCLUDED__

#define ARRAY_LENGTH(ARR) (sizeof(ARR) / sizeof((ARR)[0]))

class ParticleInfo {
    public:
        ParticleInfo();
        int sizeOfSuperParticle;
};

namespace Initializer {
    double getSizeOfSuperParticle(int, double, double);
}

namespace Utils {
    template <typename T>
    void print2DArray(T**, int, int);
}
#endif
