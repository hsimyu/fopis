#ifndef __TDPIC_POSITION_H_INCLUDED__
#define __TDPIC_POSITION_H_INCLUDED__

#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <algorithm>
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
        //! glue cell ありの座標系に乗る整数値
        int i, j, k;

        //! 実数座標
        double x, y, z;
        double dx1, dy1, dz1, dx2, dy2, dz2;

        // 実数座標からのコンストラクタ
        Position(const double _x, const double _y, const double _z){
            this->setXYZ(_x, _y, _z);
        }

        // 整数座標からのコンストラクタ
        Position(const int _i, const int _j, const int _k){
            this->setIJK(_i, _j, _k);
        }

        // 整数座標からのコンストラクタ
        Position(const unsigned int _i, const unsigned int _j, const unsigned int _k){
            this->setIJK(static_cast<int>(_i), static_cast<int>(_j), static_cast<int>(_k));
        }

        Position(void) {}

        //! Particle 用のコンストラクタはinlineで書けない
        Position(const Particle&);
        ~Position(){}

        //! コピー演算はデフォルトで良い
        Position(const Position&) = default;
        Position& operator=(const Position&) = default;
        //! ムーブ演算
        Position(Position&&) = default;
        Position& operator=(Position&&) = default;

        void setXYZ(const double _x, const double _y, const double _z){
            x = _x;
            y = _y;
            z = _z;

            i = static_cast<int>(std::floor(x)) + 1;
            j = static_cast<int>(std::floor(y)) + 1;
            k = static_cast<int>(std::floor(z)) + 1;

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

        const Position operator-(const Position& rhs) {
            return Position{x - rhs.x, y - rhs.y, z - rhs.z};
        }

        Position getReferencePosition(const Position& old) const {
            //! i が floor(x) + 1 としているため、Umeda et al. 2003 と若干形が違っていることに注意
            //! 実数座標は正しいので問題ない
            return Position{
                std::min(static_cast<double>(std::min(i, old.i)), std::max(static_cast<double>(std::max(i, old.i)) - 1.0, 0.5 * (x + old.x))),
                std::min(static_cast<double>(std::min(j, old.j)), std::max(static_cast<double>(std::max(j, old.j)) - 1.0, 0.5 * (y + old.y))),
                std::min(static_cast<double>(std::min(k, old.k)), std::max(static_cast<double>(std::max(k, old.k)) - 1.0, 0.5 * (z + old.z)))
            };
        }

        friend std::ostream& operator<<(std::ostream&, const Position&);
};

class Velocity {
    public:
        double vx, vy, vz;

        Velocity(const double _vx, const double _vy, const double _vz){
            this->set(_vx, _vy, _vz);
        }

        Velocity(const Particle&);

        void set(const double _vx, const double _vy, const double _vz){
            vx = _vx;
            vy = _vy;
            vz = _vz;
        }

        Velocity() {}
        ~Velocity(){}

        //! コピー演算
        Velocity(const Velocity&) = default;
        Velocity& operator=(const Velocity&) = default;

        //! ムーブ演算
        Velocity(Velocity&&) = default;
        Velocity& operator=(Velocity&&) = default;

        double powered(void) const {
            return vx*vx + vy*vy + vz*vz;
        }

        double abs(void) const {
            return std::sqrt(vx*vx + vy*vy + vz*vz);
        }

        friend std::ostream& operator<<(std::ostream&, const Velocity&);
};
#endif
