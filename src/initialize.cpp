#include <math.h>
#include "2dpic.h"

namespace Initializer {
    double getSizeOfSuperParticle(int nr, double density, double dx){
        return pow(dx, 2.0) * density / nr;
    }
}

