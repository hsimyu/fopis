#include "spacecraft.hpp"

const std::map<std::string, PropertyPair> Spacecraft::material_property_list = {
    { "PerfectConductor", 
        PropertyPair({
            {"RelativePermittivity", 1e11},
        })
    },
    { "Aluminum10", 
        PropertyPair({
            {"RelativePermittivity", 10.0},
        })
    },
    { "Aluminum7.5", 
        PropertyPair({
            {"RelativePermittivity", 7.5},
        })
    },
    { "Aluminum5", 
        PropertyPair({
            {"RelativePermittivity", 5.0},
        })
    },
    { "Aluminum1", 
        PropertyPair({
            {"RelativePermittivity", 5.0},
        })
    },
};
