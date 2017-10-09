#ifndef __TDPIC_PARTICLE_H_INCLUDED__
#define __TDPIC_PARTICLE_H_INCLUDED__

#include "global.hpp"
#include "environment.hpp"
#include "particle_type.hpp"

class Position;
class Velocity;

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

        //! コピー演算
        Particle(const Particle&) = default;
        Particle& operator=(const Particle&) = default;

        //! ムーブ演算
        Particle(Particle&&) = default;
        Particle& operator=(Particle&&) = default;

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

        //! isValidの操作をカプセル化
        void makeValid(void) { isValid = 1; }
        void makeInvalid(void) { isValid = 0; }

        //! Position生成
        Position getPosition(void) const;
        Position getOldPosition(void) const;
        Position getNewPosition(void) const;

        //! グリッド上の現在位置が旧位置からどの程度離れているかを返す(衝突面判定用)
        double getXMoveRatio() const;
        double getYMoveRatio() const;
        double getZMoveRatio() const;

        //! グリッドを横切る際の座標を羅列して返す
        std::vector<Position> computeCrossPoints() const;
        Position getNextXCrossPoint() const;
        Position getNextYCrossPoint() const;
        Position getNextZCrossPoint() const;

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

        //! ParticleTypeのメンバ関数呼び出しを中継
        double getCharge(void) const { return Environment::getParticleType(typeId)->getCharge(); }
        double getChargeOfSuperParticle(void) const { return Environment::getParticleType(typeId)->getChargeOfSuperParticle(); }

        friend std::ostream& operator<<(std::ostream&, Particle const&);
};
#endif
