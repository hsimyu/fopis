#include <math.h>
#include "tdpic.h"

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

        threeD_array* ex = new threeD_array(boost::extents[cx-1][cy][cz]);
        threeD_array* ey = new threeD_array(boost::extents[cx][cy-1][cz]);
        threeD_array* ez = new threeD_array(boost::extents[cx][cy][cz-1]);

        field->setEx(ex);
        field->setEy(ey);
        field->setEz(ez);

        threeD_array* bx = new threeD_array(boost::extents[cx][cy-1][cz-1]);
        threeD_array* by = new threeD_array(boost::extents[cx-1][cy][cz-1]);
        threeD_array* bz = new threeD_array(boost::extents[cx-1][cy-1][cz]);

        field->setBx(bx);
        field->setBy(by);
        field->setBz(bz);
    }
}

