#include "global.hpp"
#include "position.hpp"
#include "grid.hpp"
#include "field.hpp"
#include "particle.hpp"
#include "environment.hpp"
#include "mpiw.hpp"
#include "dataio.hpp"
#include "normalizer.hpp"
#include "utils.hpp"
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
        double qm = 0.5 * (Environment::getParticleType(pid)->getCharge()) / (Environment::getParticleType(pid)->getMass());

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
        double qm = 0.5 * (Environment::getParticleType(pid)->getCharge()) / (Environment::getParticleType(pid)->getMass());

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
    //! 先に子を2回分更新する
    for(auto& child : children) {
        child->updateParticlePosition();
        child->updateParticleVelocity();
        child->updateParticlePosition();
    }

    if (Environment::solver_type == "ES") {
        this->updateParticlePositionES();
    } else {
        this->updateParticlePositionEM();
    }
}

void RootGrid::updateParticlePositionES(void) {
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
            if( p.isValid ) {
                particles[ p.typeId ].push_back( std::move(p) );
            }
        }
    }

    //! 子グリッドへの移動をチェック
    if(this->getChildrenLength() > 0) this->checkParticlesMoveIntoChildren();
}

void RootGrid::updateParticlePositionEM(void) {
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
        const double q_per_dt = Environment::getParticleType(pid)->getCharge() / dt;

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
            if(p.isValid) {
                const int pid = p.typeId;
                const double q_per_dt = Environment::getParticleType(pid)->getChargeOfSuperParticle() / dt;
                p.distributeCurrentAtNewPosition(q_per_dt, jx, jy, jz);
                particles[ pid ].push_back( std::move(p) );
            }
        }
    }

    if(this->getChildrenLength() > 0) this->checkParticlesMoveIntoChildren();
}

void ChildGrid::updateParticlePositionES(void) {
    const double slx = dx * nx;
    const double sly = dx * ny;
    const double slz = dx * nz;

    for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
        for(auto& p : particles[pid]) {
            if(p.isValid) {
                p.updatePosition();
                checkBoundary(p, slx, sly, slz);
            }
        }
    }

    //! 子グリッドへの移動をチェック
    if(this->getChildrenLength() > 0) this->checkParticlesMoveIntoChildren();
}

void ChildGrid::updateParticlePositionEM(void) {
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
        const double q_per_dt = Environment::getParticleType(pid)->getCharge() / dt;

        for(int i = 0; i < particles[pid].size(); ++i){
            Particle& p = particles[pid][i];
            if(p.isValid) {
                p.distributeCurrentAtOldPosition(q_per_dt, jx, jy, jz);
                p.updatePosition();

                checkBoundary(p, slx, sly, slz);

                // check 後も valid であれば配分してOK
                if (p.isValid) {
                    p.distributeCurrentAtNewPosition(q_per_dt, jx, jy, jz);
                }
            }
        }
    }

    if(this->getChildrenLength() > 0) this->checkParticlesMoveIntoChildren();
}

void Grid::checkParticlesMoveIntoChildren() {
    for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
        for(auto& p : particles[pid]) {
            if (p.isValid) {
                int index = this->getChildIndexIfCovered(p.getPosition());

                if (index >= 0) {
                    this->moveParticleToChild(index, p);
                }
            }
        }
    }
}

void Grid::moveParticleToChild(int child_index, Particle& p) {
    // cout << format("It should be move to child[%d]!") % child_index << endl;
    // cout << p << endl;
    Particle new_particle = p; // コピー演算

    auto& child = children[child_index];
    new_particle.x = 2.0 * (new_particle.x - static_cast<double>(child->getFromIX()) + 1.0);
    new_particle.y = 2.0 * (new_particle.y - static_cast<double>(child->getFromIY()) + 1.0);
    new_particle.z = 2.0 * (new_particle.z - static_cast<double>(child->getFromIZ()) + 1.0);

    // cout << "...will move to new particle on child grid:" << endl;
    // cout << new_particle << endl;

    //! 親グリッド上の粒子をinvalidに
    p.makeInvalid();
    //! 子グリッド上の粒子をpush
    particles[p.typeId].push_back(std::move( new_particle ));
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
        res += 0.5 * Environment::getParticleType(pid)->getMass() * eachEnergy * Environment::getParticleType(pid)->getSize();
    }

    return res;
}

size_t Grid::getValidParticleNumber(const int pid) const {
    size_t count = 0;

    for(const auto& p : particles[pid]) {
        if (p.isValid) ++count;
    }

    return count;
}
