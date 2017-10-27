#ifndef __TDPIC_INITIALIZE_H_INCLUDED__
#define __TDPIC_INITIALIZE_H_INCLUDED__
#include <picojson.h>

class RootGrid;

namespace Initializer {
    std::shared_ptr<RootGrid> initTDPIC(const std::string& input_filename);
    void setMPIInfoToEnv(void);
    void loadEnvironment(picojson::object&);
    void loadParticleType(picojson::object&);
    void loadStaticField(picojson::object&);
}
#endif
