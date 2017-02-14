#ifndef __TDPIC_POSITION_H_INCLUDED__
#define __TDPIC_POSITION_H_INCLUDED__
#include <math.h>
class Particle;

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

        //! Particle 用のコンストラクタはinlineで書けない
        Position(const Particle&);
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

class Velocity {
    private:
        void set(const double _vx, const double _vy, const double _vz){
            vx = _vx;
            vy = _vy;
            vz = _vz;
        }

    public:
        double vx, vy, vz;

        Velocity(const double _vx, const double _vy, const double _vz){
            this->set(_vx, _vy, _vz);
        }

        //! Particle用のコンストラクタはinlineで書けない
        Velocity(const Particle&);

        ~Velocity(){}

        double getMagnitude(void) {
            return vx*vx + vy*vy + vz*vz;
        }
};
#endif
