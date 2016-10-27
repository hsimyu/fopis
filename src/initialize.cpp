#include <math.h>
#include "2dpic.h"

namespace Initializer {
    double getSizeOfSuperParticle(int nr, double density, const double dx){
        return pow(dx, 2.0) * density / nr;
    }

    void setFieldPointers(FieldPointers* field, int cx, int cy, int cz){
        auto extents = boost::extents[cx][cy][cz];
        threeD_array* phi = new threeD_array(extents);
        threeD_array* rho = new threeD_array(extents);

        field->setPhi(phi);
        field->setRho(rho);
    }
}

