#ifndef __TDPIC_PARTICLE_H_INCLUDED__
#define __TDPIC_PARTICLE_H_INCLUDED__
#include <random>
#include "global.hpp"

class Position;
class Velocity;
class Grid;

class ParticleType {
    private:
        int id;
        std::string name;

        // plasma info
        std::string type;
        double charge;
        double mass;
        double density;
        double temperature;
        int size;

        // for computation
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

            particle_per_cell = ptype.getPcell();
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
        void setTotalNumber(int _num){ totalNumber = _num; }
        void setPcell(int _pcell){ particle_per_cell = _pcell; }

        // getters
        int getId() const { return id; }
        std::string getType() const { return type; }
        std::string getName() const { return name; }
        double getCharge() const { return charge; }
        double getMass() const { return mass; }
        double getDensity() const { return density; }
        double getTemperature() const { return temperature / e; }
        double getTrueTemperature() const { return temperature; }

        int getPcell() const { return particle_per_cell; }
        int getSize() const { return size; }
        int getTotalNumber() const { return totalNumber; }

        int calcSize(void);
        int calcTotalNumber(void);
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

class Particle {
    public:
        //! - MPIで送るためにoffsetofを使う必要があり、データメンバをpublicにする必要がある
        double x, y, z;
        double vx, vy, vz;
        int typeId;

        //! - MPI_BoolがCpp Bindingでしか有効でないため、int型でisValidを保存する
        //! - 0 == false, 1 == true
        int isValid;

        Particle(const int _typeId) {
            typeId = _typeId;
            isValid = 1;
        }

        Particle(void) { isValid = 0; }

        ~Particle(){};

        // Copy Constructer
        Particle(Particle const& p){
            x = p.x;
            y = p.y;
            z = p.z;
            vx = p.vx;
            vy = p.vy;
            vz = p.vz;
            typeId = p.typeId;
            isValid = p.isValid;
        }

        Particle& operator=(Particle const& rhs){
            x = rhs.x;
            y = rhs.y;
            z = rhs.z;
            vx = rhs.vx;
            vy = rhs.vy;
            vz = rhs.vz;
            typeId = rhs.typeId;
            isValid = rhs.isValid;

            return *this;
        }

        void setPosition(const double _x, const double _y, const double _z){
            x = _x;
            y = _y;
            z = _z;
        }

        void setVelocity(const double _vx, const double _vy, const double _vz){
            vx = _vx;
            vy = _vy;
            vz = _vz;
        }

        void setVelocity(Velocity const&);
        void setPosition(Position const&);

        //! 位置の更新
        void updatePosition(void) {
            x += vx;
            y += vy;
            z += vz;
        }

        Position&& getPosition(void) const;
        Position&& getOldPosition(void) const;
        Position&& getNewPosition(void) const;

        //! 電流配分
        void distributeCurrentAtOldPosition(const double, tdArray&, tdArray&, tdArray&) const;
        void distributeCurrentAtNewPosition(const double, tdArray&, tdArray&, tdArray&) const;

        //! 個別にエネルギーを計算するためのメンバ関数
        double getEnergy(void) const;

        //! まとめてエネルギーを計算する時のためのメンバ関数
        double getSquaredMagnitudeOfVelocity(void) const {
            return (vx*vx + vy*vy + vz*vz);
        }

        double getMagnitudeOfVelocity(void) const {
            return sqrt(vx*vx + vy*vy + vz*vz);
        }

        //! 粒子生成時用の位置速度生成関数
        //! - 内部的にはParticleTypeの同名メンバ関数を呼び出すだけ
        void generateNewVelocity(void);
        void generateNewPosition(const double, const double, const double, const double, const double, const double);

        friend std::ostream& operator<<(std::ostream&, Particle const&);
};

typedef std::vector< std::vector<Particle> > ParticleArray;
#endif
