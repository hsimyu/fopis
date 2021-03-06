#ifndef __TDPIC_PARTICLE_TYPE_HPP_INCLUDED__
#define __TDPIC_PARTICLE_TYPE_HPP_INCLUDED__

#include "global.hpp"
#include "position.hpp"
#include "random_distribution.hpp"
#include <random>
#include <array>

class Grid;
class Particle;

//! 粒子系のベースクラス
class ParticleType {
    protected:
        int id;
        std::string name;

        // plasma info
        std::string type;
        double charge;
        double mass;
        double density;
        double temperature;
        int size;

        // 計算用
        int particle_per_cell;

        // random number generator
        std::mt19937 mt_x;
        std::mt19937 mt_y;
        std::mt19937 mt_z;
        std::mt19937 mt_vx;
        std::mt19937 mt_vy;
        std::mt19937 mt_vz;

        // 擬似乱数生成回数を保存しておく
        using GeneratedCount_t = std::vector<size_t>;
        GeneratedCount_t generated_counts{0, 0, 0, 0, 0, 0};
        void incrementGeneratedCount(const unsigned int id) {
            ++generated_counts[id];
        }

        void incrementPositionGeneratedCount() {
            ++generated_counts[0];
            ++generated_counts[1];
            ++generated_counts[2];
        }

        void incrementVelocityGeneratedCount() {
            ++generated_counts[3];
            ++generated_counts[4];
            ++generated_counts[5];
        }

    public:
        ParticleType();

        ParticleType(ParticleType const& ptype) {
            id = ptype.getId();
            type = ptype.getType();
            charge = ptype.getCharge();
            mass = ptype.getMass();
            density = ptype.getDensity();
            temperature = ptype.getTemperature();
            size = ptype.getSize();
        };

        // setters
        void setId(int);
        void setCharge(double _charge){ charge = _charge; }
        void setMass(double _mass){ mass = _mass; }
        void setDensity(double _density){ density = _density; }
        void setTemperature(double _temp){
            //! eV形式で入力されると仮定
            //! 内部的にはkB Teの値で持つ => eをかけて保存
            temperature = e * _temp;
        }
        void setName(std::string _name){ name = _name; }
        void setType(std::string _type){ type = _type; }
        void setSize(int _size){ size = _size; }
        void setPcell(int _pcell){ particle_per_cell = _pcell; }

        // getters
        int getId() const { return id; }
        std::string getType() const { return type; }
        std::string getName() const { return name; }
        double getCharge() const { return charge; }
        double getChargeOfSuperParticle() const { return charge * size; }
        double getMass() const { return mass; }
        double getDensity() const { return density; }
        double getTemperature() const { return temperature / e; }
        double getTrueTemperature() const { return temperature; }

        int getSize() const { return size; }
        int getPcell() const { return particle_per_cell; }
        size_t getTotalNumber() const;

        int updateSize(void);
        double calcDebyeLength(void) const;
        double calcThermalVelocity(void) const;
        double calcDeviation(void) const;
        double calcPlasmaFrequency(void) const;
        std::string calcMemory(void) const;

        auto getGeneratedCounts() const {
            return generated_counts;
        }

        void proceedGeneratedCounts(const GeneratedCount_t& target_counts);

        virtual void printInfo() const;
};

//! 背景プラズマ
class AmbientParticleType : public ParticleType {
    protected:
        Velocity drift_velocity{0.0, 0.0, 0.0};
        std::array<bool, 6> injection_axis{true, true, true, true, true, true};

    public:
        AmbientParticleType() : ParticleType(){}
        std::vector<double> calcFlux() const;

        template<typename T>
        void setDriftVelocity(const T& array) {
            drift_velocity.vx = array[0];
            drift_velocity.vy = array[1];
            drift_velocity.vz = array[2];
        };

        template<typename T>
        void setInjectionAxis(const T& array) {
            assert(array.size() == 6);

            for(int i = 0; i < 6; ++i) {
                injection_axis[i] = array[i];
            }
        };

        //! factory関数
        Particle generateNewParticle(const double min_x, const double max_x, const double min_y, const double max_y, const double min_z, const double max_z);
        Particle generateNewParticle(const double min_x, const double max_x, const double min_y, const double max_y, const double min_z, const double max_z, const Velocity& vel);
        Position generateNewPosition(const double min_x, const double max_x, const double min_y, const double max_y, const double min_z, const double max_z);
        Velocity generateNewVelocity(void);
        Velocity generateNewInjectionVelocity(AXIS axis, AXIS_SIDE low_or_up);

        virtual void printInfo() const override;
};

//! 放出用プラズマ粒子クラス
class EmissionParticleType : public ParticleType {
    protected:

    public:
        EmissionParticleType() : ParticleType(){}
};

//! 光電子
class PhotoElectronParticleType : public EmissionParticleType {
    protected:
        double current_density = 0.0;
        RandomDistribution::CosineEmission angle_gen;
    
    public:
        PhotoElectronParticleType();

        double getCurrentDensity() const {
            return current_density;
        }

        void setCurrentDensity(const double value) {
            current_density = value;
        }

        double getEmissionAmount() const;
        Particle generateNewParticle(const Position& relative_emission_position, const std::array<double, 3>& emission_vector);
        Position generateNewPosition(const Position& relative_emission_position, const std::array<double, 3>& emission_vector, const Velocity& vel);
        Velocity generateNewVelocity(const std::array<double, 3>& emission_vector);
        virtual void printInfo() const override;
};

//! 二次電子
class MaterialInfo_t;
class IncidentInfo_t;

class SecondaryParticleType : public EmissionParticleType {
    protected:
        std::mt19937_64 mt_ptype_gen;
        std::uniform_real_distribution<> dist_uniform{0.0, 1.0};

        std::mt19937_64 mt_inelastic_gen;
        RandomDistribution::CosineEmission inelastic_angle_gen;

        std::mt19937_64 mt_true_see_gen;
        RandomDistribution::CosineEmission truesee_angle_gen;

        double computeWhippleTrueSecondaryCoeff(const MaterialInfo_t& material, const double incident_energy, const double incident_angle) const;
        double computeWhippleElasticCoeff(const MaterialInfo_t& material, const double incident_energy, const double incident_angle) const;
        double computeWhippleInelasticCoeff(const MaterialInfo_t& material, const double incident_energy, const double incident_angle) const;

        double getTrueSecondaryCoeff(const MaterialInfo_t& material, const double incident_energy, const double incident_angle) const;
        double getElasticBackscatterCoeff(const MaterialInfo_t& material, const double incident_energy, const double incident_angle) const;
        double getInelasticBackscatterCoeff(const MaterialInfo_t& material, const double incident_energy, const double incident_angle) const;

        Particle generateElasticBackscatterParticle(const IncidentInfo_t& incident, const MaterialInfo_t& material);
        Particle generateInelasticBackscatterParticle(const IncidentInfo_t& incident, const MaterialInfo_t& material);
        Particle generateTrueSecondaryParticle(const IncidentInfo_t& incident, const MaterialInfo_t& material);

    public:
        SecondaryParticleType();

        std::vector<Particle> generateNewParticles(const IncidentInfo_t& incident, const MaterialInfo_t& material);
        // Position generateNewPosition(const Position& relative_emission_position, const std::array<double, 3>& emission_vector, const Velocity& vel);
        // Velocity generateNewVelocity(const std::array<double, 3>& emission_vector);
        virtual void printInfo() const override;
};

//! 放出されるビーム用クラス
class BeamParticleType : public EmissionParticleType {
    protected:
        std::string emission_type;
        double accel_potential;
        double beam_current;
        double beam_divergence;
        double emission_radius;

    public:
        BeamParticleType() : EmissionParticleType(), emission_radius(0.0) {}

        void setEmissionType(const std::string& emiss_type) {
            emission_type = emiss_type;
        }

        std::string getEmissionType() const {return emission_type;}

        void setAcceleratingPotential(const double value) {
            accel_potential = value;
        }

        double getAcceleratingPotential() const {return accel_potential;}
        double getAcceleration() const;

        void setBeamCurrent(const double value) { beam_current = value; }
        double getBeamCurrent() const {return beam_current;}

        void setBeamDivergence(const double value) { beam_divergence = value; }
        double getBeamDivergence() const {return beam_divergence;}

        void setEmissionRadius(const double value) { emission_radius = value; }
        double setEmissionRadius() const {return emission_radius;}

        double getEmissionAmount() const;

        Particle generateNewParticle(const Position& relative_emission_position, const std::array<double, 3>& emission_vector);
        Position generateNewPosition(const Position& relative_emission_position, const std::array<double, 3>& emission_vector, const Velocity& vel);
        Velocity generateNewVelocity(const std::array<double, 3>& emission_vector);
        virtual void printInfo() const override;
};
#endif
