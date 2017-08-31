#ifndef __TDPIC_INITIALIZE_H_INCLUDED__
#define __TDPIC_INITIALIZE_H_INCLUDED__
#include <picojson.h>

class RootGrid;

namespace Initializer {
    std::shared_ptr<RootGrid> initTDPIC();
    void setMPIInfoToEnv(void);
    void loadEnvironment(picojson::object&);
    void loadParticleType(picojson::object&);
}
#endif
