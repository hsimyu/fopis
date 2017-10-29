#include "spacecraft.hpp"
#include "particle.hpp"
#include "particle_type.hpp"

bool Spacecraft::hasSecondaryParticles() const {
    for(const auto& pinfo : emit_particle_info) {
        const auto emit_ptype_ptr = Environment::getEmissionParticleType(pinfo.first);
        if (emit_ptype_ptr->getType() == "secondary") return true;
    }

    return false;
}

void Spacecraft::addIncidentEvent(
    const std::shared_ptr<ParticleType> ptype_ptr,
    const double remaining_time,
    const Position& incident_pos,
    const Velocity& incident_vel,
    AXIS axis) {
    IncidentInfo_t new_incident{
        ptype_ptr->getMass(),
        ptype_ptr->getCharge(),
        remaining_time,
        incident_pos,
        incident_vel,
        axis
    };

    incident_events.push_back( std::move(new_incident) );
}

void Spacecraft::clearIncidentEvents() {
    incident_events.clear();
}

void Spacecraft::distributeInnerParticleChargeForSecondary(Particle& p, AXIS axis) {
    const auto id = p.typeId;
    const auto q = p.getChargeOfSuperParticle();

    const int isign = (p.vx > 0.0) ? 1 : -1;
    const int jsign = (p.vy > 0.0) ? 1 : -1;
    const int ksign = (p.vz > 0.0) ? 1 : -1;
    bool surface_is_found = false;

    if (this->hasSecondaryParticles()) {
        //! 二次電子がある場合はIncidentの情報保存を行う必要がある
        auto cross_points = p.computeCrossPoints();
        for(auto& cross_point : cross_points) {
            if (axis != AXIS::x && isXsurfacePoint(cross_point, isign)) {
                distributeInnerParticleChargeToXsurface(cross_point, id, q);
                //! move_x_ratioが前のグリッドからの移動割合
                this->addIncidentEvent(p.getParticleTypePtr(), p.getXMoveRatio(), cross_point, p.getVelocity(), AXIS::x);
                surface_is_found = true;
                break;
            } else if (axis != AXIS::y && isYsurfacePoint(cross_point, jsign)) {
                distributeInnerParticleChargeToYsurface(cross_point, id, q);
                this->addIncidentEvent(p.getParticleTypePtr(), p.getYMoveRatio(), cross_point, p.getVelocity(), AXIS::y);
                surface_is_found = true;
                break;
            } else if (axis != AXIS::z && isZsurfacePoint(cross_point, ksign)) {
                distributeInnerParticleChargeToZsurface(cross_point, id, q);
                this->addIncidentEvent(p.getParticleTypePtr(), p.getZMoveRatio(), cross_point, p.getVelocity(), AXIS::z);
                surface_is_found = true;
                break;
            }
        }

        if (!surface_is_found) {
            cout << "[ERROR] " << Environment::rankStr() <<
                "Surface cannot detect on Spacecraft::distributeInnerParticleChargeForSecondary():\n";
            cout << p << '\n';
            cout << "Cross Points:\n";
            for(auto& cross_point : cross_points) {
                cout << cross_point << "\n";
            }
            cout << endl;
        }
    } else {
        auto cross_points = p.computeCrossPoints();
        for(auto& cross_point : cross_points) {
            if (axis != AXIS::x && isXsurfacePoint(cross_point, isign)) {
                distributeInnerParticleChargeToXsurface(cross_point, id, q);
                surface_is_found = true;
                break;
            } else if (axis != AXIS::y && isYsurfacePoint(cross_point, jsign)) {
                distributeInnerParticleChargeToYsurface(cross_point, id, q);
                surface_is_found = true;
                break;
            } else if (axis != AXIS::z && isZsurfacePoint(cross_point, ksign)) {
                distributeInnerParticleChargeToZsurface(cross_point, id, q);
                surface_is_found = true;
                break;
            }
        }

        if (!surface_is_found) {
            cout << "[ERROR] " << Environment::rankStr() <<
                "Surface cannot detect on Spacecraft::distributeInnerParticleChargeForSecondary():\n";
            cout << p << '\n';
            cout << "Cross Points:\n";
            for(auto& cross_point : cross_points) {
                cout << cross_point << "\n";
            }
            cout << endl;
        }
    }

    current[id] += q;
    p.makeInvalid();
}
