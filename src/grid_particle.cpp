#include "global.hpp"
#include "position.hpp"
#include "grid.hpp"
#include "field.hpp"
#include "particle.hpp"
#include "environment.hpp"
#include "mpiw.hpp"
#include "dataio.hpp"
#include "normalizer.hpp"
#include <random>
#include <algorithm>

void Grid::updateParticleVelocity(void) {
    if (Environment::solver_type == "ES") {
        this->updateParticleVelocityES();
    } else {
        this->updateParticleVelocityEM();
    }
}

void Grid::updateParticleVelocityES(void) {
    tdArray& exref = field->getExRef();
    tdArray& eyref = field->getEyRef();
    tdArray& ezref = field->getEzRef();

    for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
        double qm = 0.5 * (Environment::ptype[pid].getCharge()) / (Environment::ptype[pid].getMass());

        for(int pnum = 0; pnum < particles[pid].size(); ++pnum){
            Particle& p = particles[pid][pnum];
            if(p.isValid) {
                //! 毎回生成するよりコピーのが早い？
                Position pos(p);

                int i = pos.i, j = pos.j, k = pos.k;
                double v1 = qm * pos.dx2 * pos.dy2 * pos.dz2;
                double v2 = qm * pos.dx1 * pos.dy2 * pos.dz2;
                double v3 = qm * pos.dx2 * pos.dy1 * pos.dz2;
                double v4 = qm * pos.dx1 * pos.dy1 * pos.dz2;
                double v5 = qm * pos.dx2 * pos.dy2 * pos.dz1;
                double v6 = qm * pos.dx1 * pos.dy2 * pos.dz1;
                double v7 = qm * pos.dx2 * pos.dy1 * pos.dz1;
                double v8 = qm * pos.dx1 * pos.dy1 * pos.dz1;

                const double ex =  v1*exref[i][j][k]
                    + v2*exref[i+1][j][k]
                    + v3*exref[i][j+1][k]
                    + v4*exref[i+1][j+1][k]
                    + v5*exref[i][j][k+1]
                    + v6*exref[i+1][j][k+1]
                    + v7*exref[i][j+1][k+1]
                    + v8*exref[i+1][j+1][k+1];
                const double ey =  v1*eyref[i][j][k]
                    + v2*eyref[i+1][j][k]
                    + v3*eyref[i][j+1][k]
                    + v4*eyref[i+1][j+1][k]
                    + v5*eyref[i][j][k+1]
                    + v6*eyref[i+1][j][k+1]
                    + v7*eyref[i][j+1][k+1]
                    + v8*eyref[i+1][j+1][k+1];
                const double ez =  v1*ezref[i][j][k]
                    + v2*ezref[i+1][j][k]
                    + v3*ezref[i][j+1][k]
                    + v4*ezref[i+1][j+1][k]
                    + v5*ezref[i][j][k+1]
                    + v6*ezref[i+1][j][k+1]
                    + v7*ezref[i][j+1][k+1]
                    + v8*ezref[i+1][j+1][k+1];

                const double bx = 0.0;
                const double by = 0.0;
                const double bz = 0.0;
                double boris = 2.0/(1.0 + (bx*bx+by*by+bz*bz));

                double vx1 = p.vx + ex;
                double vy1 = p.vy + ey;
                double vz1 = p.vz + ez;

                double vxt = vx1 + vy1*bz - vz1 * by;
                double vyt = vy1 + vz1*bx - vx1 * bz;
                double vzt = vz1 + vx1*by - vy1 * bx;

                p.vx = vx1 + ex + boris*(vyt*bz - vzt*by);
                p.vy = vy1 + ey + boris*(vzt*bx - vxt*bz);
                p.vz = vz1 + ez + boris*(vxt*by - vyt*bx);
            }
        }
    }
}

void Grid::updateParticleVelocityEM(void) {
    tdArray& exref = field->getExRef();
    tdArray& eyref = field->getEyRef();
    tdArray& ezref = field->getEzRef();
    tdArray& bxref = field->getBxRef();
    tdArray& byref = field->getByRef();
    tdArray& bzref = field->getBzRef();

    for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
        double qm = 0.5 * (Environment::ptype[pid].getCharge()) / (Environment::ptype[pid].getMass());

        for(int pnum = 0; pnum < particles[pid].size(); ++pnum){
            Particle& p = particles[pid][pnum];
            if(p.isValid) {
                //! 毎回生成するよりコピーのが早い？
                Position pos(p);

                int i = pos.i, j = pos.j, k = pos.k;
                double v1 = qm * pos.dx2 * pos.dy2 * pos.dz2;
                double v2 = qm * pos.dx1 * pos.dy2 * pos.dz2;
                double v3 = qm * pos.dx2 * pos.dy1 * pos.dz2;
                double v4 = qm * pos.dx1 * pos.dy1 * pos.dz2;
                double v5 = qm * pos.dx2 * pos.dy2 * pos.dz1;
                double v6 = qm * pos.dx1 * pos.dy2 * pos.dz1;
                double v7 = qm * pos.dx2 * pos.dy1 * pos.dz1;
                double v8 = qm * pos.dx1 * pos.dy1 * pos.dz1;

                const double ex =  v1*exref[i][j][k]
                    + v2*exref[i+1][j][k]
                    + v3*exref[i][j+1][k]
                    + v4*exref[i+1][j+1][k]
                    + v5*exref[i][j][k+1]
                    + v6*exref[i+1][j][k+1]
                    + v7*exref[i][j+1][k+1]
                    + v8*exref[i+1][j+1][k+1];
                const double ey =  v1*eyref[i][j][k]
                    + v2*eyref[i+1][j][k]
                    + v3*eyref[i][j+1][k]
                    + v4*eyref[i+1][j+1][k]
                    + v5*eyref[i][j][k+1]
                    + v6*eyref[i+1][j][k+1]
                    + v7*eyref[i][j+1][k+1]
                    + v8*eyref[i+1][j+1][k+1];
                const double ez =  v1*ezref[i][j][k]
                    + v2*ezref[i+1][j][k]
                    + v3*ezref[i][j+1][k]
                    + v4*ezref[i+1][j+1][k]
                    + v5*ezref[i][j][k+1]
                    + v6*ezref[i+1][j][k+1]
                    + v7*ezref[i][j+1][k+1]
                    + v8*ezref[i+1][j+1][k+1];
                const double bx =  v1*bxref[i][j][k]
                    + v2*bxref[i+1][j][k]
                    + v3*bxref[i][j+1][k]
                    + v4*bxref[i+1][j+1][k]
                    + v5*bxref[i][j][k+1]
                    + v6*bxref[i+1][j][k+1]
                    + v7*bxref[i][j+1][k+1]
                    + v8*bxref[i+1][j+1][k+1];
                const double by =  v1*byref[i][j][k]
                    + v2*byref[i+1][j][k]
                    + v3*byref[i][j+1][k]
                    + v4*byref[i+1][j+1][k]
                    + v5*byref[i][j][k+1]
                    + v6*byref[i+1][j][k+1]
                    + v7*byref[i][j+1][k+1]
                    + v8*byref[i+1][j+1][k+1];
                const double bz =  v1*bzref[i][j][k]
                    + v2*bzref[i+1][j][k]
                    + v3*bzref[i][j+1][k]
                    + v4*bzref[i+1][j+1][k]
                    + v5*bzref[i][j][k+1]
                    + v6*bzref[i+1][j][k+1]
                    + v7*bzref[i][j+1][k+1]
                    + v8*bzref[i+1][j+1][k+1];
                double boris = 2.0/(1.0 + (bx*bx+by*by+bz*bz));

                double vx1 = p.vx + ex;
                double vy1 = p.vy + ey;
                double vz1 = p.vz + ez;

                double vxt = vx1 + vy1*bz - vz1 * by;
                double vyt = vy1 + vz1*bx - vx1 * bz;
                double vzt = vz1 + vx1*by - vy1 * bx;

                p.vx = vx1 + ex + boris*(vyt*bz - vzt*by);
                p.vy = vy1 + ey + boris*(vzt*bx - vxt*bz);
                p.vz = vz1 + ez + boris*(vxt*by - vyt*bx);
            }
        }
    }
}

void Grid::updateParticlePosition(void) {
    if (Environment::solver_type == "ES") {
        this->updateParticlePositionES();
    } else {
        this->updateParticlePositionEM();
    }
}

void Grid::updateParticlePositionES(void) {
    std::vector< std::vector<Particle> > pbuff(6);
    std::vector< std::vector<Particle> > pbuffRecv(6);

    const double slx = dx * nx;
    const double sly = dx * ny;
    const double slz = dx * nz;

    for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
        for(auto& p : particles[pid]) {
            if(p.isValid) {
                p.updatePosition();
                checkXBoundary(pbuff, p, slx);
                checkYBoundary(pbuff, p, sly);
                checkZBoundary(pbuff, p, slz);
            }
        }
    }

    //! 対象の方向のプロセス数が 1 かつ 周期境界でない場合には通信しなくてよい
    if( (Environment::proc_x > 1) || Environment::isNotBoundary(AXIS::x, AXIS_SIDE::low)) {
        MPIw::Environment::sendRecvParticlesX(pbuff, pbuffRecv);

        for(int axis = 0; axis < 2; ++axis) {
            for(int i = 0; i < pbuffRecv[axis].size(); ++i){
                Particle& p = pbuffRecv[axis][i];

                // 0方向(下側)からやってきた時、領域長分だけ座標をずらす
                if(axis == 0) {
                    p.x -= slx;
                } else {
                    p.x += slx;
                }

                checkYBoundary(pbuff, p, sly);
                checkZBoundary(pbuff, p, slz);
            }
        }
    }

    if( (Environment::proc_y > 1) || Environment::isNotBoundary(AXIS::y, AXIS_SIDE::low)) {
        MPIw::Environment::sendRecvParticlesY(pbuff, pbuffRecv);

        for(int axis = 2; axis < 4; ++axis) {
            for(int i = 0; i < pbuffRecv[axis].size(); ++i){
                Particle& p = pbuffRecv[axis][i];

                // 0方向(下側)からやってきた時、領域長分だけ座標をずらす
                if(axis == 2) {
                    p.y -= sly;
                } else {
                    p.y += sly;
                }
                checkZBoundary(pbuff, p, slz);
            }
        }
    }

    if( (Environment::proc_z > 1) || Environment::isNotBoundary(AXIS::z, AXIS_SIDE::low)) {
        MPIw::Environment::sendRecvParticlesZ(pbuff, pbuffRecv);

        for(int axis = 4; axis < 6; ++axis) {
            for(int i = 0; i < pbuffRecv[axis].size(); ++i){
                Particle& p = pbuffRecv[axis][i];

                // 0方向(下側)からやってきた時、領域長分だけ座標をずらす
                if(axis == 4) {
                    p.z -= slz;
                } else {
                    p.z += slz;
                }
            }
        }
    }

    for(int i = 0; i < 6; ++i) {
        for(auto& p : pbuffRecv[i]){
#ifdef DEBUG
            if( particles[ p.typeId ].capacity() == particles[ p.typeId ].size() ) {
                cout << format("[WARNING] The size of %s array is full.: capacity = %d, size = %d") % Environment::ptype[ p.typeId ].getName() % particles[ p.typeId ].capacity() % particles[ p.typeId ].capacity()<< endl;
            }
#endif
            if(p.isValid) {
                particles[ p.typeId ].push_back( std::move(p) );
            }
        }
    }

    this->injectParticles();
}

void Grid::updateParticlePositionEM(void) {
    std::vector< std::vector<Particle> > pbuff(6);
    std::vector< std::vector<Particle> > pbuffRecv(6);

    const double slx = dx * nx;
    const double sly = dx * ny;
    const double slz = dx * nz;

    tdArray& jx = field->getJx();
    tdArray& jy = field->getJy();
    tdArray& jz = field->getJz();

    // 電流配列を初期化
    field->initializeCurrent(dt);

    for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
        //! note: 1/dxdydz 分を係数としてかけておく必要がある？
        const double q_per_dt = Environment::ptype[pid].getCharge() / dt;

        for(int i = 0; i < particles[pid].size(); ++i){
            Particle& p = particles[pid][i];
            if(p.isValid) {
                p.distributeCurrentAtOldPosition(q_per_dt, jx, jy, jz);
                p.updatePosition();

                checkXBoundary(pbuff, p, slx);
                checkYBoundary(pbuff, p, sly);
                checkZBoundary(pbuff, p, slz);

                // check 後も valid であれば配分してOK
                if (p.isValid) {
                    p.distributeCurrentAtNewPosition(q_per_dt, jx, jy, jz);
                }
            }
        }
    }

    // 電流通信
    MPIw::Environment::sendRecvField(jx);
    MPIw::Environment::sendRecvField(jy);
    MPIw::Environment::sendRecvField(jz);

    if( (Environment::proc_x > 1) || Environment::isNotBoundary(AXIS::x, AXIS_SIDE::low)) {
        MPIw::Environment::sendRecvParticlesX(pbuff, pbuffRecv);

        for(int axis = 0; axis < 2; ++axis) {
            for(int i = 0; i < pbuffRecv[axis].size(); ++i){
                Particle& p = pbuffRecv[axis][i];

                // 0方向(下側)からやってきた時、領域長分だけ座標をずらす
                if(axis == 0) {
                    p.x -= slx;
                } else {
                    p.x += slx;
                }

                checkYBoundary(pbuff, p, sly);
                checkZBoundary(pbuff, p, slz);
            }
        }
    }

    if( (Environment::proc_y > 1) || Environment::isNotBoundary(AXIS::y, AXIS_SIDE::low)) {
        MPIw::Environment::sendRecvParticlesY(pbuff, pbuffRecv);

        for(int axis = 2; axis < 4; ++axis) {
            for(int i = 0; i < pbuffRecv[axis].size(); ++i){
                Particle& p = pbuffRecv[axis][i];

                // 0方向(下側)からやってきた時、領域長分だけ座標をずらす
                if(axis == 2) {
                    p.y -= sly;
                } else {
                    p.y += sly;
                }
                checkZBoundary(pbuff, p, slz);
            }
        }
    }

    if( (Environment::proc_z > 1) || Environment::isNotBoundary(AXIS::z, AXIS_SIDE::low)) {
        MPIw::Environment::sendRecvParticlesZ(pbuff, pbuffRecv);

        for(int axis = 4; axis < 6; ++axis) {
            for(int i = 0; i < pbuffRecv[axis].size(); ++i){
                Particle& p = pbuffRecv[axis][i];

                // 0方向(下側)からやってきた時、領域長分だけ座標をずらす
                if(axis == 4) {
                    p.z -= slz;
                } else {
                    p.z += slz;
                }
            }
        }
    }

    for(int i = 0; i < 6; ++i) {
        for(auto& p : pbuffRecv[i]){
#ifdef DEBUG
            if( particles[ p.typeId ].capacity() == particles[ p.typeId ].size() ) {
                cout << format("[WARNING] The size of %s array is full.: capacity = %d, size = %d") % Environment::ptype[ p.typeId ].getName() % particles[ p.typeId ].capacity() % particles[ p.typeId ].capacity()<< endl;
            }
#endif
            if(p.isValid) {
                const int pid = p.typeId;
                const double q_per_dt = Environment::ptype[pid].getCharge() / dt;
                p.distributeCurrentAtNewPosition(q_per_dt, jx, jy, jz);
                particles[ pid ].push_back( std::move(p) );
            }
        }
    }

    this->injectParticles();
}

//! 粒子の位置から電荷を空間電荷にする
//! 基本的にはroot_gridに対してのみ呼ぶ
void Grid::updateRho() {
    RhoArray& rho = field->getRho();

#ifdef CHARGE_CONSERVATION
    // 電荷保存則をcheckするため、古いrhoを保持する
    auto old_rho = rho;
#endif

    //! rhoを初期化
    Utils::initializeRhoArray(rho);

    //! 物体上の一時的な情報を初期化
    for(auto& obj : objects) {
        if(obj.isDefined()) {
            obj.resetCurrent();
        }
    }

    for(int pid = 0; pid < Environment::num_of_particle_types; ++pid){
        double q = Environment::ptype[pid].getCharge();
        const auto rho_idx = pid + 1;

        for(auto& p : particles[pid]) {
            if(p.isValid) {
                for(auto& obj : objects) {
                    //! 物体中にいた場合には自動的に invalid になる
                    obj.distributeInnerParticleCharge(p);
                }

                //! もし物体内でなければ
                if (p.isValid) {
                    const auto pos = p.getPosition();
                    const int i = pos.i, j = pos.j, k = pos.k;

                    rho[rho_idx][i  ][j  ][k] += pos.dx2 * pos.dy2 * pos.dz2 * q;
                    rho[rho_idx][i+1][j  ][k] += pos.dx1 * pos.dy2 * pos.dz2 * q;
                    rho[rho_idx][i  ][j+1][k] += pos.dx2 * pos.dy1 * pos.dz2 * q;
                    rho[rho_idx][i+1][j+1][k] += pos.dx1 * pos.dy1 * pos.dz2 * q;
                    rho[rho_idx][i  ][j  ][k+1] += pos.dx2 * pos.dy2 * pos.dz1 * q;
                    rho[rho_idx][i+1][j  ][k+1] += pos.dx1 * pos.dy2 * pos.dz1 * q;
                    rho[rho_idx][i  ][j+1][k+1] += pos.dx2 * pos.dy1 * pos.dz1 * q;
                    rho[rho_idx][i+1][j+1][k+1] += pos.dx1 * pos.dy1 * pos.dz1 * q;
                }
            }
        }
    }


    for(auto& obj : objects) {
        if (obj.isDefined()) {
            //! 物体に配分された電荷を現在の rho に印加する
            obj.applyCharge(rho);
        }
    }

    for(int pid = 1; pid < Environment::num_of_particle_types + 1; ++pid) {
        for(int i = 1; i < nx + 2; ++i) {
            for(int j = 1; j < ny + 2; ++j) {
                for(int k = 1; k < nz + 2; ++k) {
                    rho[0][i][j][k] += rho[pid][i][j][k];
                }
            }
        }
    }

    //! rho を隣に送る
    for(int pid = 0; pid < Environment::num_of_particle_types + 1; ++pid) {
        MPIw::Environment::sendRecvField(rho[pid]);
    }

    //! 物体が存在する場合、電荷再配分が必要になる
    if (objects.size() > 0) {
        //! 一度 Poisson を解いて phi を更新
        solvePoisson();

        for(auto& obj : objects) {
            if (obj.isDefined()) obj.redistributeCharge(rho, field->getPhi());
        }

        MPIw::Environment::sendRecvField(rho[0]);
    }

#ifdef CHARGE_CONSERVATION
    if (Environment::solver_type == "EM") {
        field->checkChargeConservation(old_rho, 1.0, dx);
    }
#endif

}

void Grid::injectParticles(void) {
    static std::vector< std::vector<double> > residual;
    static bool isFirstCall = true;

    //! staticな残余変数の初期化
    if(isFirstCall) {
        residual.resize(Environment::num_of_particle_types);

        for(int i = 0; i < residual.size(); ++i) {
            residual[i].resize(6);
            for(int j = 0; j < 6; ++j) {
                residual[i][j] = 0.0;
            }
        }

        isFirstCall = false;
    }

    //! - 粒子位置の上限を設定
    double max_x = static_cast<double>(Environment::cell_x);
    double max_y = static_cast<double>(Environment::cell_y);
    double max_z = static_cast<double>(Environment::cell_z);

    //! - 上側境界にいる場合は外側にはみ出した粒子を生成しないようにする
    if(!Environment::isNotBoundary(AXIS::x, AXIS_SIDE::up)) max_x -= 1.0;
    if(!Environment::isNotBoundary(AXIS::y, AXIS_SIDE::up)) max_y -= 1.0;
    if(!Environment::isNotBoundary(AXIS::z, AXIS_SIDE::up)) max_z -= 1.0;

    for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
        std::vector<double> flux = Environment::ptype[pid].calcFlux(*this);

        if(!Environment::isNotBoundary(AXIS::x, AXIS_SIDE::low)) {
            const int index = 0;
            const int inject_num = floor(dt * flux[index] + residual[pid][index]);
            residual[pid][index] += dt * flux[index] - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Particle p(pid);
                p.generateNewVelocity();

                //! 流入方向速度に変換
                //! 実際はフラックスを積分して割合を求める必要がある
                while(p.vx <= 0.0) p.generateNewVelocity();

                p.generateNewPosition(0.0, p.vx * dt, 0.0, max_y, 0.0, max_z);
                particles[pid].push_back( std::move(p) );
            }
        }

        if(!Environment::isNotBoundary(AXIS::x, AXIS_SIDE::up)) {
            const int index = 1;
            const int inject_num = floor(dt * flux[index] + residual[pid][index]);
            residual[pid][index] += dt * flux[index] - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Particle p(pid);
                p.generateNewVelocity();

                while(p.vx >= 0.0) p.generateNewVelocity();

                //! 負方向速度をxの最大値から引いた点までがありうる範囲
                p.generateNewPosition(max_x + p.vx * dt, max_x, 0.0, max_y, 0.0, max_z);
                particles[pid].push_back( std::move(p) );
            }
        }

        if(!Environment::isNotBoundary(AXIS::y, AXIS_SIDE::low)) {
            const int index = 2;
            const int inject_num = floor(dt * flux[index] + residual[pid][index]);
            residual[pid][index] += dt * flux[index] - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Particle p(pid);
                p.generateNewVelocity();

                while(p.vy <= 0.0) p.generateNewVelocity();

                p.generateNewPosition(0.0, max_x, 0.0, p.vy * dt, 0.0, max_z);
                particles[pid].push_back( std::move(p) );
            }
        }

        if(!Environment::isNotBoundary(AXIS::y, AXIS_SIDE::up)) {
            const int index = 3;
            const int inject_num = floor(dt * flux[index] + residual[pid][index]);
            residual[pid][index] += dt * flux[index] - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Particle p(pid);
                p.generateNewVelocity();

                while(p.vy >= 0.0) p.generateNewVelocity();

                p.generateNewPosition(0.0, max_x, max_y + p.vy * dt, max_y, 0.0, max_z);
                particles[pid].push_back( std::move(p) );
            }
        }

        if(!Environment::isNotBoundary(AXIS::z, AXIS_SIDE::low)) {
            const int index = 4;
            const int inject_num = floor(dt * flux[index] + residual[pid][index]);
            residual[pid][index] += dt * flux[index] - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Particle p(pid);
                p.generateNewVelocity();

                while(p.vz <= 0.0) p.generateNewVelocity();

                p.generateNewPosition(0.0, max_x, 0.0, max_y, 0.0, p.vz * dt);
                particles[pid].push_back( std::move(p) );
            }
        }

        if(!Environment::isNotBoundary(AXIS::z, AXIS_SIDE::up)) {
            const int index = 5;
            const int inject_num = floor(dt * flux[index] + residual[pid][index]);
            residual[pid][index] += dt * flux[index] - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Particle p(pid);
                p.generateNewVelocity();

                while(p.vz >= 0.0) p.generateNewVelocity();

                p.generateNewPosition(0.0, max_x, 0.0, max_y, max_z + p.vz * dt, max_z);
                particles[pid].push_back( std::move(p) );
            }
        }
    }
}

double Grid::getParticleEnergy(void) const {
    double res = 0.0;

    for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
        double eachEnergy = 0.0;
        for(int i = 0; i < particles[pid].size(); ++i){
            if(particles[pid][i].isValid) {
                eachEnergy += particles[pid][i].getSquaredMagnitudeOfVelocity();
            }
        }
        res += 0.5 * Environment::ptype[pid].getMass() * eachEnergy * Environment::ptype[pid].getSize();
    }

    return res;
}
