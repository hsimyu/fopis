#include "spacecraft.hpp"
#include "particle.hpp"

bool Spacecraft::hasSecondaryParticles() const {
    for(const auto& pinfo : emit_particle_info) {
        const auto emit_ptype_ptr = Environment::getEmissionParticleType(pinfo.first);
        if (emit_ptype_ptr->getType() == "secondary") return true;
    }

    return false;
}

void Spacecraft::addIncidentEvent(const Position& incident_pos, const Velocity& incident_vel) {
    IncidentInfo_t new_incident{incident_pos, incident_vel};
    incident_events.push_back( std::move(new_incident) );
}

void Spacecraft::clearIncidentEvents() {
    incident_events.clear();
}