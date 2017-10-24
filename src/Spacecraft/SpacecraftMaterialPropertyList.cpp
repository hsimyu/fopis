#include "spacecraft.hpp"

const std::map<std::string, PropertyPair> Spacecraft::material_property_list = {
    { "Silicon", 
        PropertyPair({
            {"RelativePermittivity", 12.0},
        })
    },
    { "Aluminum", 
        PropertyPair({
            {"RelativePermittivity", 9.5},
        })
    },
    { "EpoxyGlass", 
        PropertyPair({
            {"RelativePermittivity", 5.0},
        })
    },
    { "Bakelite", 
        PropertyPair({
            {"RelativePermittivity", 4.8},
        })
    },
    { "GlassReinforcedTeflon", 
        PropertyPair({
            {"RelativePermittivity", 2.32},
        })
    },
    { "Teflon", 
        PropertyPair({
            {"RelativePermittivity", 2.1},
        })
    },
};
