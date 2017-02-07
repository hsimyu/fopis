#ifndef __TDPIC_PARTICLE_H_INCLUDED__
#define __TDPIC_PARTICLE_H_INCLUDED__
#include <vector>
#include <string>
#include "global.hpp"

class Position;

class Particle {
    public:
        // MPIで送るためにoffsetofを使う必要があり、
        // publicにする必要がある
        int typeId;
        int isValid;
        double x, y, z;
        double vx, vy, vz;

        Particle(void){
            isValid = 1;
        }

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
        //! Position class はまだ定義されていないのでinlineで書けない
        void setPosition(Position const&);

        void setVelocity(const double _vx, const double _vy, const double _vz){
            vx = _vx;
            vy = _vy;
            vz = _vz;
        }

        //! 位置の更新
        void updatePosition(void) {
            x += vx;
            y += vy;
            z += vz;
        }

        double getEnergy(void) const;
        double getSquaredMagnitudeOfVelocity(void) const;

        friend std::ostream& operator<<(std::ostream&, Particle const&);
};

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
    public:
        ParticleType(){}

        // setters
        void setId(int _id){ id = _id; }
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
        double calcDeviation(void) const;
        std::string calcMemory(void) const;

        friend std::ostream& operator<<(std::ostream&, const ParticleType&);
        friend std::istream& operator<<(std::istream&, const ParticleType&);
};

typedef std::vector< std::vector<Particle> > ParticleArray;

class Position {
    private:
        void updateDelta(void){
            // delta xyz is automatically updated
            dx1 = x - (i - 1);
            dy1 = y - (j - 1);
            dz1 = z - (k - 1);
            dx2 = 1.0 - dx1;
            dy2 = 1.0 - dy1;
            dz2 = 1.0 - dz1;
        }

    public:
        int i, j, k;
        double x, y, z;
        double dx1, dy1, dz1, dx2, dy2, dz2;

        // Position Class
        Position(const double _x, const double _y, const double _z){
            this->setXYZ(_x, _y, _z);
        }

        Position(const int _i, const int _j, const int _k){
            this->setIJK(_i, _j, _k);
        }

        Position(const Particle& p){
            this->setXYZ(p.x, p.y, p.z);
        }

        ~Position(){}

        void setXYZ(const double _x, const double _y, const double _z){
            x = _x;
            y = _y;
            z = _z;

            i = floor(x) + 1;
            j = floor(y) + 1;
            k = floor(z) + 1;

            this->updateDelta();
        }

        void setIJK(const int _i, const int _j, const int _k){
            i = _i;
            j = _j;
            k = _k;

            // x, y, zはglue cellなしの座標系にする
            x = static_cast<double>(i - 1);
            y = static_cast<double>(j - 1);
            z = static_cast<double>(k - 1);

            this->updateDelta();
        }
};

#endif
