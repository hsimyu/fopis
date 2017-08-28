#ifndef __TDPIC_PARTICLE_TYPE_HPP_INCLUDED__
#define __TDPIC_PARTICLE_TYPE_HPP_INCLUDED__

#include <random>
#include "position.hpp"

class Position;
class Velocity;
class Grid;

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
        int totalNumber;

        // random number generator
        std::mt19937 mt_x;
        std::mt19937 mt_y;
        std::mt19937 mt_z;
        std::mt19937 mt_vx;
        std::mt19937 mt_vy;
        std::mt19937 mt_vz;

    public:
        ParticleType(void);
        ParticleType(ParticleType const& ptype) {
            id = ptype.getId();
            type = ptype.getType();
            charge = ptype.getCharge();
            mass = ptype.getMass();
            density = ptype.getDensity();
            temperature = ptype.getTemperature();
            size = ptype.getSize();

            // particle_per_cell = ptype.getPcell();
            totalNumber = ptype.getTotalNumber();
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
        void setTotalNumber(int _num){ totalNumber = _num; }

        // getters
        int getId() const { return id; }
        std::string getType() const { return type; }
        std::string getName() const { return name; }
        double getCharge() const { return charge; }
        double getMass() const { return mass; }
        double getDensity() const { return density; }
        double getTemperature() const { return temperature / e; }
        double getTrueTemperature() const { return temperature; }

        int getSize() const { return size; }
        int getPcell() const { return particle_per_cell; }
        int getTotalNumber() const { return totalNumber; }

        int updateSize(void);
        int updateTotalNumber(void);
        double calcDebyeLength(void) const;
        double calcThermalVelocity(void) const;
        double calcDeviation(void) const;
        double calcPlasmaFrequency(void) const;
        std::vector<double> calcFlux(Grid const&) const;
        std::string calcMemory(void) const;

        Position generateNewPosition(const double min_x, const double max_x, const double min_y, const double max_y, const double min_z, const double max_z);
        Velocity generateNewVelocity(void);

        friend std::ostream& operator<<(std::ostream&, const ParticleType&);
        friend std::istream& operator<<(std::istream&, const ParticleType&);
};

//! 各種派生クラス

class AmbientParticle : public ParticleType {
    protected:

    public:
        AmbientParticle() : ParticleType(){}


};

class BeamParticle : public ParticleType {
    protected:
        Position emission_position;
        std::vector<double> emission_vector;

    public:
        BeamParticle() : ParticleType(){}

        void setEmissionPosition(const std::vector<double>& pos) {
            if (pos.size() != 3) {
                throw std::invalid_argument("The size of the argument for setEmissionPosition() must be 3.");
            }

            emission_position = Position(pos[0], pos[1], pos[2]);
        }
        Position getEmissionPosition() const {return emission_position;}

        void setEmissionVector(const std::vector<double>& pos) {
            if (pos.size() != 3) {
                throw std::invalid_argument("The size of the argument for setEmissionVector() must be 3.");
            }

            emission_vector = pos;
        }
        std::vector<double> getEmissionVector() const {return emission_vector;}

        friend std::ostream& operator<<(std::ostream&, const BeamParticle&);
};
#endif
