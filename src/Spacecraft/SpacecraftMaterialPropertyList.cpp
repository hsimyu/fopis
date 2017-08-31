#include "spacecraft.hpp"

const std::map<std::string, PropertyPair> Spacecraft::material_property_list = {
    { "PerfectConductor", 
        PropertyPair({
            {"Capacitance", 0.0},
            {"Resistance", 0.0},
        })
    },
    { "material_name_1", 
        PropertyPair({
            {"Capacitance", 20.0},
            {"Resistance", 10.0},
        })
    },
};
