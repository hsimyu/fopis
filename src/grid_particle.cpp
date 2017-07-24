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
        const double q_per_dt = Environment::ptype[pid].getCharge() / Utils::Normalizer::normalizeTime(Environment::dt);

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

    if(Environment::proc_x > 1) {
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

    if(Environment::proc_y > 1) {
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

    if(Environment::proc_z > 1) {
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
                particles[ p.typeId ].push_back(p);
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

    for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
        //! note: 1/dxdydz 分を係数としてかけておく必要がある？
        const double q_per_dt = Environment::ptype[pid].getCharge() / Utils::Normalizer::normalizeTime(Environment::dt);

        for(int i = 0; i < particles[pid].size(); ++i){
            Particle& p = particles[pid][i];
            if(p.isValid) {
                Position old_pos(p);
                p.updatePosition();
                Position new_pos(p);
                Position ref_pos = new_pos.getReferencePosition(old_pos);

                // charge flux の計算
                const double flux1[3] = {
                    q_per_dt * (ref_pos.x - old_pos.x),
                    q_per_dt * (ref_pos.y - old_pos.y),
                    q_per_dt * (ref_pos.z - old_pos.z),
                };

                const double flux2[3] = {
                    q_per_dt * (new_pos.x - ref_pos.x),
                    q_per_dt * (new_pos.y - ref_pos.y),
                    q_per_dt * (new_pos.z - ref_pos.z),
                };

                //! 1次の shape factor の計算
                //! 中点における shape factor を取る
                //! グリッドをまたいだ時 ref_position の dx1 は 1.0 となるので、
                const double average_w1[3] = {
                    0.5 * (ref_pos.x + old_pos.x) - (static_cast<double>(old_pos.i) - 1.0),
                    0.5 * (ref_pos.y + old_pos.y) - (static_cast<double>(old_pos.j) - 1.0),
                    0.5 * (ref_pos.z + old_pos.z) - (static_cast<double>(old_pos.k) - 1.0),
                };
                const double average_w2[3] = {
                    0.5 * (ref_pos.x + new_pos.x) - (static_cast<double>(new_pos.i) - 1.0),
                    0.5 * (ref_pos.y + new_pos.y) - (static_cast<double>(new_pos.j) - 1.0),
                    0.5 * (ref_pos.z + new_pos.z) - (static_cast<double>(new_pos.k) - 1.0),
                };

                int pos_i = old_pos.i;
                int pos_j = old_pos.j;
                int pos_k = old_pos.k;
                
                // Current Conserving Weighting, l = 1 (old cell)
                jx[pos_i][pos_j    ][pos_k    ] += flux1[0] * (1.0 - average_w1[1]) * (1.0 - average_w1[2]);
                jx[pos_i][pos_j + 1][pos_k    ] += flux1[0] * average_w1[1] * (1.0 - average_w1[2]);
                jx[pos_i][pos_j    ][pos_k + 1] += flux1[0] * (1.0 - average_w1[1]) * average_w1[2];
                jx[pos_i][pos_j + 1][pos_k + 1] += flux1[0] * average_w1[1] * average_w1[2];

                jy[pos_i    ][pos_j][pos_k    ] += flux1[1] * (1.0 - average_w1[0]) * (1.0 - average_w1[2]);
                jy[pos_i + 1][pos_j][pos_k    ] += flux1[1] * average_w1[0] * (1.0 - average_w1[2]);
                jy[pos_i    ][pos_j][pos_k + 1] += flux1[1] * (1.0 - average_w1[0]) * average_w1[2];
                jy[pos_i + 1][pos_j][pos_k + 1] += flux1[1] * average_w1[0] * average_w1[2];

                jz[pos_i    ][pos_j    ][pos_k] += flux1[2] * (1.0 - average_w1[0]) * (1.0 - average_w1[1]);
                jz[pos_i + 1][pos_j    ][pos_k] += flux1[2] * average_w1[0] * (1.0 - average_w1[1]);
                jz[pos_i    ][pos_j + 1][pos_k] += flux1[2] * (1.0 - average_w1[0]) * average_w1[1];
                jz[pos_i + 1][pos_j + 1][pos_k] += flux1[2] * average_w1[0] * average_w1[1];

                pos_i = new_pos.i;
                pos_j = new_pos.j;
                pos_k = new_pos.k;

                // Current Conserving Weighting, l = 2 (new cell)
                jx[pos_i][pos_j    ][pos_k    ] += flux2[0] * (1.0 - average_w2[1]) * (1.0 - average_w2[2]);
                jx[pos_i][pos_j + 1][pos_k    ] += flux2[0] * average_w2[1] * (1.0 - average_w2[2]);
                jx[pos_i][pos_j    ][pos_k + 1] += flux2[0] * (1.0 - average_w2[1]) * average_w2[2];
                jx[pos_i][pos_j + 1][pos_k + 1] += flux2[0] * average_w2[1] * average_w2[2];

                jy[pos_i    ][pos_j][pos_k    ] += flux2[1] * (1.0 - average_w2[0]) * (1.0 - average_w2[2]);
                jy[pos_i + 1][pos_j][pos_k    ] += flux2[1] * average_w2[0] * (1.0 - average_w2[2]);
                jy[pos_i    ][pos_j][pos_k + 1] += flux2[1] * (1.0 - average_w2[0]) * average_w2[2];
                jy[pos_i + 1][pos_j][pos_k + 1] += flux2[1] * average_w2[0] * average_w2[2];

                jz[pos_i    ][pos_j    ][pos_k] += flux2[2] * (1.0 - average_w2[0]) * (1.0 - average_w2[1]);
                jz[pos_i + 1][pos_j    ][pos_k] += flux2[2] * average_w2[0] * (1.0 - average_w2[1]);
                jz[pos_i    ][pos_j + 1][pos_k] += flux2[2] * (1.0 - average_w2[0]) * average_w2[1];
                jz[pos_i + 1][pos_j + 1][pos_k] += flux2[2] * average_w2[0] * average_w2[1];

                checkXBoundary(pbuff, p, slx);
                checkYBoundary(pbuff, p, sly);
                checkZBoundary(pbuff, p, slz);
            }
        }
    }

    if(Environment::proc_x > 1) {
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

    if(Environment::proc_y > 1) {
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

    if(Environment::proc_z > 1) {
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
    if(!Environment::isPeriodic(AXIS::x, AXIS_SIDE::up)) max_x -= 1.0;
    if(!Environment::isPeriodic(AXIS::y, AXIS_SIDE::up)) max_y -= 1.0;
    if(!Environment::isPeriodic(AXIS::z, AXIS_SIDE::up)) max_z -= 1.0;

    for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
        std::vector<double> flux = Environment::ptype[pid].calcFlux(*this);

        if(!Environment::isPeriodic(AXIS::x, AXIS_SIDE::low)) {
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
                particles[pid].push_back(p);
                inject_parray[pid].push_back(p);
            }
        }

        if(!Environment::isPeriodic(AXIS::x, AXIS_SIDE::up)) {
            const int index = 1;
            const int inject_num = floor(dt * flux[index] + residual[pid][index]);
            residual[pid][index] += dt * flux[index] - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Particle p(pid);
                p.generateNewVelocity();

                while(p.vx >= 0.0) p.generateNewVelocity();

                //! 負方向速度をxの最大値から引いた点までがありうる範囲
                p.generateNewPosition(max_x + p.vx * dt, max_x, 0.0, max_y, 0.0, max_z);
                particles[pid].push_back(p);
                inject_parray[pid].push_back(p);
            }
        }

        if(!Environment::isPeriodic(AXIS::y, AXIS_SIDE::low)) {
            const int index = 2;
            const int inject_num = floor(dt * flux[index] + residual[pid][index]);
            residual[pid][index] += dt * flux[index] - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Particle p(pid);
                p.generateNewVelocity();

                while(p.vy <= 0.0) p.generateNewVelocity();

                p.generateNewPosition(0.0, max_x, 0.0, p.vy * dt, 0.0, max_z);
                particles[pid].push_back(p);
                inject_parray[pid].push_back(p);
            }
        }

        if(!Environment::isPeriodic(AXIS::y, AXIS_SIDE::up)) {
            const int index = 3;
            const int inject_num = floor(dt * flux[index] + residual[pid][index]);
            residual[pid][index] += dt * flux[index] - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Particle p(pid);
                p.generateNewVelocity();

                while(p.vy >= 0.0) p.generateNewVelocity();

                p.generateNewPosition(0.0, max_x, max_y + p.vy * dt, max_y, 0.0, max_z);
                particles[pid].push_back(p);
                inject_parray[pid].push_back(p);
            }
        }

        if(!Environment::isPeriodic(AXIS::z, AXIS_SIDE::low)) {
            const int index = 4;
            const int inject_num = floor(dt * flux[index] + residual[pid][index]);
            residual[pid][index] += dt * flux[index] - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Particle p(pid);
                p.generateNewVelocity();

                while(p.vz <= 0.0) p.generateNewVelocity();

                p.generateNewPosition(0.0, max_x, 0.0, max_y, 0.0, p.vz * dt);
                particles[pid].push_back(p);
                inject_parray[pid].push_back(p);
            }
        }

        if(!Environment::isPeriodic(AXIS::z, AXIS_SIDE::up)) {
            const int index = 5;
            const int inject_num = floor(dt * flux[index] + residual[pid][index]);
            residual[pid][index] += dt * flux[index] - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Particle p(pid);
                p.generateNewVelocity();

                while(p.vz >= 0.0) p.generateNewVelocity();

                p.generateNewPosition(0.0, max_x, 0.0, max_y, max_z + p.vz * dt, max_z);
                particles[pid].push_back(p);
                inject_parray[pid].push_back(p);
            }
        }
    }

    // IO::plotParticleEnergyDistribution(inject_parray, "inject_");
    // IO::plotParticleVelocityDistribution(inject_parray, "inject_");
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
