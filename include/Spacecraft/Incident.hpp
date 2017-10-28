#ifndef __TDPIC_SPACECRAFT_INCIDENT_H_INCLUDED__
#define __TDPIC_SPACECRAFT_INCIDENT_H_INCLUDED__

#include "global.hpp"
#include "position.hpp"

//! 物体表面への衝突に関する情報を格納する構造体
class IncidentInfo_t {
    private:
        Position pos;
        Velocity vel;
        AXIS axis;
        double time;
        double mass;
        double charge;

    public:
        IncidentInfo_t(const double mass, const double charge, const double time, const Position& _pos, const Velocity& _vel);
        IncidentInfo_t(const double mass, const double charge, const double time, const Position& _pos, const Velocity& _vel, AXIS axis);

        //! eV単位で返す
        double getIncidentEnergyInElectronVolt() const;
        double getIncidentEnergy() const;
        double getIncidentAngle() const;
        double getRemainingTime() const {return time;}

        void setSurface(AXIS axis);
        bool isXsurfaceIncident() const;
        bool isYsurfaceIncident() const;
        bool isZsurfaceIncident() const;

        Position getPosition() const;
        double getX() const;
        double getY() const;
        double getZ() const;

        Velocity getVelocity() const;
        double getVx() const;
        double getVy() const;
        double getVz() const;
};

#endif