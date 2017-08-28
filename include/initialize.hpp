#ifndef __TDPIC_INITIALIZE_H_INCLUDED__
#define __TDPIC_INITIALIZE_H_INCLUDED__
#include <picojson.h>

class Grid;

namespace Initializer {
    void initTDPIC(Grid*&);
    void setMPIInfoToEnv(void);
    void loadEnvironment(picojson::object&);
    void loadParticleType(picojson::object&);
    Grid* initializeGrid(void);
}
#endif
