#include "global.hpp"
#include "position.hpp"
#include "grid.hpp"
#include "field.hpp"
#include "particle.hpp"
#include "environment.hpp"
#include "mpiw.hpp"
#include "dataio.hpp"
#include <silo.h>
#include <random>
#include <algorithm>

void Grid::updateParticleVelocity(void) {
    tdArray& exref = field->getExRef();
    tdArray& eyref = field->getEyRef();
    tdArray& ezref = field->getEzRef();

    for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
        double qm = 0.5 * (Environment::ptype[pid].getCharge()) * (Environment::ptype[pid].getMass());

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

                double ex =  v1*exref[i][j][k]
                    + v2*exref[i+1][j][k]
                    + v3*exref[i][j+1][k]
                    + v4*exref[i+1][j+1][k]
                    + v5*exref[i][j][k+1]
                    + v6*exref[i+1][j][k+1]
                    + v7*exref[i][j+1][k+1]
                    + v8*exref[i+1][j+1][k+1];
                double ey =  v1*eyref[i][j][k]
                    + v2*eyref[i+1][j][k]
                    + v3*eyref[i][j+1][k]
                    + v4*eyref[i+1][j+1][k]
                    + v5*eyref[i][j][k+1]
                    + v6*eyref[i+1][j][k+1]
                    + v7*eyref[i][j+1][k+1]
                    + v8*eyref[i+1][j+1][k+1];
                double ez =  v1*ezref[i][j][k]
                    + v2*ezref[i+1][j][k]
                    + v3*ezref[i][j+1][k]
                    + v4*ezref[i+1][j+1][k]
                    + v5*ezref[i][j][k+1]
                    + v6*ezref[i+1][j][k+1]
                    + v7*ezref[i][j+1][k+1]
                    + v8*ezref[i+1][j+1][k+1];
                double bx = 0.0;
                double by = 0.0;
                double bz = 0.0;
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
    std::vector< std::vector<Particle> > pbuff(6);
    std::vector< std::vector<Particle> > pbuffRecv(6);

    const double slx = dx * nx;
    const double sly = dx * ny;
    const double slz = dx * nz;

    for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
        for(int i = 0; i < particles[pid].size(); ++i){
            Particle& p = particles[pid][i];
            if(p.isValid) {
                p.updatePosition();
                checkXBoundary(pbuff, p, slx);
                checkYBoundary(pbuff, p, sly);
                checkZBoundary(pbuff, p, slz);
            }
        }
    }

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

    MPIw::Environment::sendRecvParticlesY(pbuff, pbuffRecv);

    for(int axis = 2; axis < 4; ++axis) {
        for(int i = 0; i < pbuffRecv[axis].size(); ++i){
            Particle& p = pbuffRecv[axis][i];

            // 0方向(下側)からやってきた時、領域長分だけ座標をずらす
            if(axis == 2) {
                p.x -= sly;
            } else {
                p.x += sly;
            }
            checkZBoundary(pbuff, p, slz);
        }
    }

    MPIw::Environment::sendRecvParticlesZ(pbuff, pbuffRecv);

    for(int axis = 4; axis < 6; ++axis) {
        for(int i = 0; i < pbuffRecv[axis].size(); ++i){
            Particle& p = pbuffRecv[axis][i];

            // 0方向(下側)からやってきた時、領域長分だけ座標をずらす
            if(axis == 4) {
                p.x -= slz;
            } else {
                p.x += slz;
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
                particles[ p.typeId ].push_back(p);
            }
        }
    }

    this->injectParticles();
}

void Grid::injectParticles(void) {
    static std::vector< std::vector<double> > residual;
    static bool isFirstCall = true;

    ParticleArray inject_parray;
    inject_parray.resize(2);

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

    //! (正規化された) dt = 1と仮定
    const double dt = 1.0;

    //! - 粒子位置の上限を設定
    double max_x = static_cast<double>(Environment::cell_x);
    double max_y = static_cast<double>(Environment::cell_y);
    double max_z = static_cast<double>(Environment::cell_z);
    //! - 上側境界にいる場合は外側にはみ出した粒子を生成しないようにする
    if(Environment::onHighXedge) max_x -= 1.0;
    if(Environment::onHighYedge) max_y -= 1.0;
    if(Environment::onHighZedge) max_z -= 1.0;

    for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
        std::vector<double> flux = Environment::ptype[pid].calcFlux(*this);

        if(Environment::onLowXedge) {
            const int index = 0;
            const int inject_num = floor(dt * flux[index] + residual[pid][index]);
            // cout << Environment::rankStr() << "inject from -x = " << inject_num << ", flux = " << flux[index] << ", residual = " << residual[pid][index] << endl;
            residual[pid][index] += dt * flux[index] - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Particle p(pid);
                p.generateNewVelocity();

                //! 流入方向速度に変換
                //! 実際はフラックスを積分して割合を求める必要がある
                while(p.vx <= 0.0) p.generateNewVelocity();
                // if(p.vx < 0.0) p.vx = -p.vx;

                p.generateNewPosition(0.0, p.vx * dt, 0.0, max_y, 0.0, max_z);
                particles[pid].push_back(p);
                inject_parray[pid].push_back(p);
            }
        }

        if(Environment::onHighXedge) {
            const int index = 1;
            const int inject_num = floor(dt * flux[index] + residual[pid][index]);
            // cout << Environment::rankStr() << "inject from +x = " << inject_num << ", flux = " << flux[index] << ", residual = " << residual[pid][index] << endl;
            residual[pid][index] += dt * flux[index] - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Particle p(pid);
                p.generateNewVelocity();

                //! 流入方向速度に変換
                //! 実際はフラックスを積分して割合を求める必要がある
                // if(p.vx > 0.0) p.vx = -p.vx;
                while(p.vx >= 0.0) p.generateNewVelocity();

                //! 負方向速度をxの最大値から引いた点までがありうる範囲
                p.generateNewPosition(max_x + p.vx * dt, max_x, 0.0, max_y, 0.0, max_z);
                particles[pid].push_back(p);
                inject_parray[pid].push_back(p);
            }
        }

        if(Environment::onLowYedge) {
            const int index = 2;
            const int inject_num = floor(dt * flux[index] + residual[pid][index]);
            // cout << Environment::rankStr() << "inject from -y = " << inject_num << ", flux = " << flux[index] << endl;
            residual[pid][index] += dt * flux[index] - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Particle p(pid);
                p.generateNewVelocity();

                //! 流入方向速度に変換
                //! 実際はフラックスを積分して割合を求める必要がある
                // if(p.vy < 0.0) p.vy = -p.vy;
                while(p.vy <= 0.0) p.generateNewVelocity();

                p.generateNewPosition(0.0, max_x, 0.0, p.vy * dt, 0.0, max_z);
                particles[pid].push_back(p);
                inject_parray[pid].push_back(p);
            }
        }

        if(Environment::onHighYedge) {
            const int index = 3;
            const int inject_num = floor(dt * flux[index] + residual[pid][index]);
            // cout << Environment::rankStr() << "inject from +y = " << inject_num << ", flux = " << flux[index] << endl;
            residual[pid][index] += dt * flux[index] - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Particle p(pid);
                p.generateNewVelocity();

                //! 流入方向速度に変換
                //! 実際はフラックスを積分して割合を求める必要がある
                // if(p.vy > 0.0) p.vy = -p.vy;
                while(p.vy >= 0.0) p.generateNewVelocity();

                p.generateNewPosition(0.0, max_x, max_y + p.vy * dt, max_y, 0.0, max_z);
                particles[pid].push_back(p);
                inject_parray[pid].push_back(p);
            }
        }

        if(Environment::onLowZedge) {
            const int index = 4;
            const int inject_num = floor(dt * flux[index] + residual[pid][index]);
            // cout << Environment::rankStr() << "inject from -z = " << inject_num << ", flux = " << flux[index] << endl;
            residual[pid][index] += dt * flux[index] - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Particle p(pid);
                p.generateNewVelocity();

                //! 流入方向速度に変換
                //! 実際はフラックスを積分して割合を求める必要がある
                // if(p.vz < 0.0) p.vz = -p.vz;
                while(p.vz <= 0.0) p.generateNewVelocity();

                p.generateNewPosition(0.0, max_x, 0.0, max_y, 0.0, p.vz * dt);
                particles[pid].push_back(p);
                inject_parray[pid].push_back(p);
            }
        }

        if(Environment::onHighZedge) {
            const int index = 5;
            const int inject_num = floor(dt * flux[index] + residual[pid][index]);
            // cout << Environment::rankStr() << "inject from +z = " << inject_num << ", flux = " << flux[index] << endl;
            residual[pid][index] += dt * flux[index] - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Particle p(pid);
                p.generateNewVelocity();

                //! 流入方向速度に変換
                //! 実際はフラックスを積分して割合を求める必要がある
                // if(p.vz > 0.0) p.vz = -p.vz;
                while(p.vz >= 0.0) p.generateNewVelocity();

                p.generateNewPosition(0.0, max_x, 0.0, max_y, max_z + p.vz * dt, max_z);
                particles[pid].push_back(p);
                inject_parray[pid].push_back(p);
            }
        }
    }

    IO::plotParticleEnergyDistribution(inject_parray, "inject_");
    IO::plotParticleVelocityDistribution(inject_parray, "inject_");
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
