#include "normalizer.hpp"
#include <Spacecraft/Incident.hpp>

IncidentInfo_t::IncidentInfo_t(const double _mass, const double _charge, const Position& _pos, const Velocity& _vel) : pos{_pos}, vel{_vel}, mass{_mass}, charge{_charge} {}

IncidentInfo_t::IncidentInfo_t(const double _mass, const double _charge, const Position& _pos, const Velocity& _vel, AXIS _axis) : pos{_pos}, vel{_vel}, axis{_axis}, mass{_mass}, charge{_charge} {}

double IncidentInfo_t::getIncidentEnergyInElectronVolt() const {
    return Normalizer::unnormalizeEnergy(0.5 * mass * vel.powered(), "eV");
}

//! 法線からの角度を返す
double IncidentInfo_t::getIncidentAngle() const {
    switch(axis) {
        case AXIS::x:
            return std::acos( std::fabs(vel.vx) / vel.abs() );
        case AXIS::y:
            return std::acos( std::fabs(vel.vy) / vel.abs() );
        case AXIS::z:
            return std::acos( std::fabs(vel.vz) / vel.abs() );
    }
    throw std::logic_error("[ERROR] Incident Axis is set to Unknown value.");
};

void IncidentInfo_t::setSurface(AXIS _axis) {
    axis = _axis;
};

bool IncidentInfo_t::isXsurfaceIncident() const {
    return (axis == AXIS::x); 
}

bool IncidentInfo_t::isYsurfaceIncident() const {
    return (axis == AXIS::y); 
}

bool IncidentInfo_t::isZsurfaceIncident() const {
    return (axis == AXIS::z); 
}

inline Position IncidentInfo_t::getPosition() const { return pos; }
inline double IncidentInfo_t::getX() const {return pos.x;}
inline double IncidentInfo_t::getY() const {return pos.y;}
inline double IncidentInfo_t::getZ() const {return pos.z;}

inline Velocity IncidentInfo_t::getVelocity()) const { return vel; }
inline double IncidentInfo_t::getVx() const {return vel.vx;}
inline double IncidentInfo_t::getVy() const {return vel.vy;}
inline double IncidentInfo_t::getVz() const {return vel.vz;}