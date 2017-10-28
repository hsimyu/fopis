#ifndef __TDPIC_SPACECRAFT_INCIDENT_H_INCLUDED__
#define __TDPIC_SPACECRAFT_INCIDENT_H_INCLUDED__

#include "position.hpp"

//! 物体表面への衝突に関する情報を格納する構造体
class IncidentInfo_t {
    private:
        Position pos;
        Velocity vel;
        std::array<bool, 3> surface_vector{false, false, false};

    public:
        IncidentInfo_t(const Position& _pos, const Velocity& _vel) : pos{_pos}, vel{_vel} {}

        double energy;
        double angle;

        bool isXsurfaceIncident() const {
            return surface_vector[0]; 
        }

        bool isYsurfaceIncident() const {
            return surface_vector[1]; 
        }

        bool isZsurfaceIncident() const {
            return surface_vector[2]; 
        }

        double getX() const {return pos.x;}
        double getY() const {return pos.y;}
        double getZ() const {return pos.z;}

        double getVx() const {return vel.vx;}
        double getVy() const {return vel.vy;}
        double getVz() const {return vel.vz;}
};

#endif