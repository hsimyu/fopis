#include "grid.hpp"
#include "field.hpp"
#include "utils.hpp"

//! child grid constructor
//! 渡されたGridを親とした子グリッドを生成する
ChildGrid::ChildGrid(Grid* g, const int _from_ix, const int _from_iy, const int _from_iz, const int _to_ix,   const int _to_iy,   const int _to_iz) : Grid() {
    constexpr double refineRatio = 2.0;

    parent = std::shared_ptr<Grid>(g);
    level = g->getLevel() + 1;

    //! @{
    //! 子グリッドの場合, from_ix, to_ix変数は純粋に親グリッドの何番目に乗っているかを表す
    //! Glue cell 分も考慮した方に乗る
    from_ix = _from_ix;
    from_iy = _from_iy;
    from_iz = _from_iz;
    to_ix = _to_ix;
    to_iy = _to_iy;
    to_iz = _to_iz;
    base_x = g->getBaseX() + g->getDx() * (from_ix - 1);
    base_y = g->getBaseY() + g->getDx() * (from_iy - 1);
    base_z = g->getBaseZ() + g->getDx() * (from_iz - 1);
    //! @}

    // patchの大きさを計算
    nx = (_to_ix - _from_ix) * 2 + 1;
    ny = (_to_iy - _from_iy) * 2 + 1;
    nz = (_to_iz - _from_iz) * 2 + 1;

    // refineRatioは2で固定
    dx = g->getDx() / refineRatio;
    dt = g->getDt() / refineRatio;

    this->checkGridValidness();

    // Field初期化
    this->initializeField();
}

void ChildGrid::checkGridValidness() {
    const double refineRatio = 2.0;

    bool isValid = true;

    if(from_ix == 0 || from_iy == 0 || from_iz == 0) {
        std::cerr << "[ERROR] Base index of child patch cannot be defined as 0 (0 is glue cell)." << endl;
        isValid = false;
    }

    // x extent
    if( to_ix > (parent->getNX() + 1)){
        std::cerr << "[ERROR] A child patch's x-extent exceeds the parent's extent. : " << to_ix << " > " << (parent->getNX() + 1)<< endl;
        isValid = false;
    }

    // y extent
    if( to_iy > (parent->getNY() + 1)){
        std::cerr << "[ERROR] A child patch's y-extent exceeds the parent's extent. : " << to_iy << " > " << (parent->getNY() + 1)<< endl;
        isValid = false;
    }

    // z extent
    if( to_iz > (parent->getNZ() + 1)){
        std::cerr << "[ERROR] A child patch's z-extent exceeds the parent's extent. : " << to_iz << " > " << (parent->getNZ() + 1)<< endl;
        isValid = false;
    }

    if(!isValid) MPIw::Environment::abort(1);
}

int ChildGrid::getXNodeSize(void) const {
    //! 周期境界の場合は上側境界と下側境界の間の空間も有効な空間となるので、
    //! 上側のノードを1つ増やす
    return (level == 0 && Environment::isNotBoundary(AXIS::x, AXIS_SIDE::up)) ? nx + 1 : nx;
}

int ChildGrid::getYNodeSize(void) const {
    return (level == 0 && Environment::isNotBoundary(AXIS::y, AXIS_SIDE::up)) ? ny + 1 : ny;
}

int ChildGrid::getZNodeSize(void) const {
    return (level == 0 && Environment::isNotBoundary(AXIS::z, AXIS_SIDE::up)) ? nz + 1 : nz;
}

void ChildGrid::solvePoisson(void) {
    const int DEFAULT_ITERATION_LOOP = 500;
    field->solvePoisson(DEFAULT_ITERATION_LOOP, dx);
}

void ChildGrid::updateEfield(void) {
    field->updateEfield(dx);
}

void ChildGrid::updateEfieldFDTD(void) {
    field->updateEfieldFDTD(dx, dt);
}

void ChildGrid::updateBfield(void) {
    field->updateBfield(dx, nx, ny, nz, dt);
}

void ChildGrid::checkXBoundary(ParticleArray& pbuff, Particle& p, const double slx) {
    if(p.x < 0.0) {
        if (Environment::isNotBoundary(AXIS::x, AXIS_SIDE::low)) {
            // 計算空間の端でない場合は粒子を隣へ送る
            // 計算空間の端にいるが、周期境界の場合も粒子を送る必要がある -> isNotBoundary()でまとめて判定できる
            pbuff[0].push_back(p);
        }
        p.makeInvalid();
    } else if (p.x > (slx - dx)) {
        if (Environment::isNotBoundary(AXIS::x, AXIS_SIDE::up)) {
            //! 計算空間の端でない場合、slx - dx から slx までの空間は下側の空間が担当するため、 slx を超えた場合のみ粒子を送信する
            //! また、計算空間の端にいるが、周期境界の場合も同様の処理でよいため、isNotBoundary()でまとめて判定できる
            if (p.x > slx) {
                pbuff[1].push_back(p);
                p.makeInvalid();
            }
        } else {
            p.makeInvalid();
        }
    }
}

void ChildGrid::checkYBoundary(ParticleArray& pbuff, Particle& p, const double sly) {
    if(p.y < 0.0) {
        if (Environment::isNotBoundary(AXIS::y, AXIS_SIDE::low)) pbuff[2].push_back(p);
        p.makeInvalid();
    } else if (p.y > (sly - dx)) {
        if (Environment::isNotBoundary(AXIS::y, AXIS_SIDE::up)) {
            if (p.y > sly) {
                pbuff[3].push_back(p);
                p.makeInvalid();
            }
        } else {
            p.makeInvalid();
        }
    }
}

void ChildGrid::checkZBoundary(ParticleArray& pbuff, Particle& p, const double slz) {
    if(p.z < 0.0) {
        if (Environment::isNotBoundary(AXIS::z, AXIS_SIDE::low)) pbuff[4].push_back(p);
        p.makeInvalid();
    } else if (p.z > (slz - dx)) {
        if (Environment::isNotBoundary(AXIS::z, AXIS_SIDE::up)) {
            if (p.z > slz) {
                pbuff[5].push_back(p);
                p.makeInvalid();
            }
        } else {
            p.makeInvalid();
        }
    }
}

//! 粒子の位置から電荷を空間電荷にする
void ChildGrid::updateRho() {
    RhoArray& rho = field->getRho();

#ifdef CHARGE_CONSERVATION
    // 電荷保存則をcheckするため、古いrhoを保持する
    auto old_rho = rho;
#endif

    //! rhoを初期化
    Utils::initializeRhoArray(rho);

    for(int pid = 0; pid < Environment::num_of_particle_types; ++pid){
        double q = Environment::getParticleType(pid)->getChargeOfSuperParticle();
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