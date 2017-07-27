#include "global.hpp"
#include "environment.hpp"
#include "field.hpp"
#include "utils.hpp"
#include "mpiw.hpp"
#include <algorithm>

//! Poissonソルバを呼び出す
void Field::solvePoisson(const int loopnum, const double dx) {
    this->solvePoissonPSOR(loopnum, dx);
    // this->solvePoissonJacobi(loopnum, dx);
}

//! @brief SOR法
void Field::solvePoissonPSOR(const int loopnum, const double dx) {
    const double omega = 2.0/(1.0 + sin(M_PI/(phi.shape()[0] - 2))); // spectral radius
    const double rho_coeff = pow(dx, 2) / Utils::Normalizer::normalizeEpsilon(eps0);

    const int cx_with_glue = phi.shape()[0];
    const int cy_with_glue = phi.shape()[1];
    const int cz_with_glue = phi.shape()[2];

    constexpr double required_error = 1.0e-7;

    this->setBoundaryConditionPhi();

    const bool is_periodic[6] = {
        Environment::isPeriodic(AXIS::x, AXIS_SIDE::low),
        Environment::isPeriodic(AXIS::x, AXIS_SIDE::up),
        Environment::isPeriodic(AXIS::y, AXIS_SIDE::low),
        Environment::isPeriodic(AXIS::y, AXIS_SIDE::up),
        Environment::isPeriodic(AXIS::z, AXIS_SIDE::low),
        Environment::isPeriodic(AXIS::z, AXIS_SIDE::up),
    };

    for(int loop = 1; loop <= loopnum; ++loop) {
        //! is_periodicは各方向の境界について
        //! 「その境界は計算空間の境界でない」か「その境界が計算空間の境界であり、周期境界である」場合にtrueとなるため
        //! 各方向の端要素の計算をIterationで計算する場合にチェックする必要がある
        for(int k = 1; k < cz_with_glue - 1; ++k){
            if((k != 1 || is_periodic[4]) && (k != cz_with_glue - 2 || is_periodic[5])) {
                for(int j = 1; j < cy_with_glue - 1; ++j){
                    if((j != 1 || is_periodic[2]) && (j != cy_with_glue - 2 || is_periodic[3])) {
                        for(int i = 1; i < cx_with_glue - 1; ++i){
                            if((i != 1 || is_periodic[0]) && (i != cx_with_glue - 2 || is_periodic[2])) {
                                phi[i][j][k] = (1.0 - omega) * phi[i][j][k] + omega*(phi[i+1][j][k] + phi[i-1][j][k] + phi[i][j+1][k] + phi[i][j-1][k] + phi[i][j][k+1] + phi[i][j][k-1] + rho_coeff * rho[i][j][k])/6.0;
                            }
                        }
                    }
                }
            }
        }

        //! @note: 実際には部分部分をSORで計算して送信というのを繰り返す方が収束効率がよい
        MPIw::Environment::sendRecvField(phi);

        if ( (loop % 10 == 0) && (this->checkPhiResidual() < required_error) ) break;
    }

    this->setBoundaryConditionPhi();
}

void Field::setBoundaryConditionPhi(void) {
    const std::vector< AXIS > axisArray{AXIS::x, AXIS::y, AXIS::z};
    const std::vector< AXIS_SIDE > luArray{AXIS_SIDE::low, AXIS_SIDE::up};

    for (auto axis : axisArray) {
        for (auto low_or_up : luArray) {
            std::string boundary = Environment::getBoundaryCondition(axis, low_or_up);

            if (boundary == "D") {
                this->setDirichletPhi(axis, low_or_up);
            } else if (boundary == "N") {
                this->setNeumannPhi(axis, low_or_up);
            } else if (boundary == "P") {
                // 周期境界の場合は特に何もしない
                // 代わりにiteration時にsetする
            } else {
                throw std::invalid_argument( 
                    (format("Unknown boundary condition type was used: boundary = %s") % boundary).str()
                );
            }
        }
    }
}

void Field::setNeumannPhi(const AXIS axis, const AXIS_SIDE low_or_up) {
    const int max_index = phi.shape()[ Utils::getAxisIndex(axis) ];

    const int boundary_index = (low_or_up == AXIS_SIDE::low) ? 1 : (max_index - 2);
    const int neighbor_element_index = (low_or_up == AXIS_SIDE::low) ? 2 : (max_index - 3);

    const int cx_with_glue = phi.shape()[0];
    const int cy_with_glue = phi.shape()[1];
    const int cz_with_glue = phi.shape()[2];

    if ( Environment::isOnEdge(axis, low_or_up) ) {
        if (axis == AXIS::x) {
            for(int k = 1; k < cz_with_glue - 1; ++k){
                for(int j = 1; j < cy_with_glue - 1; ++j){
                    phi[boundary_index][j][k] = phi[neighbor_element_index][j][k];
                }
            }
        } else if (axis == AXIS::y) {
            for(int k = 1; k < cz_with_glue - 1; ++k){
                for(int i = 1; i < cx_with_glue - 1; ++i){
                    phi[i][boundary_index][k] = phi[i][neighbor_element_index][k];
                }
            }
        } else {
            for(int j = 1; j < cy_with_glue - 1; ++j){
                for(int i = 1; i < cx_with_glue - 1; ++i){
                    phi[i][j][boundary_index] = phi[i][j][neighbor_element_index];
                }
            }
        }
    }
}

void Field::setDirichletPhi(const AXIS axis, const AXIS_SIDE low_or_up) {
    const int max_index = phi.shape()[ Utils::getAxisIndex(axis) ];
    const int boundary_index = (low_or_up == AXIS_SIDE::low) ? 1 : (max_index - 2);

    const int cx_with_glue = phi.shape()[0];
    const int cy_with_glue = phi.shape()[1];
    const int cz_with_glue = phi.shape()[2];

    if ( Environment::isOnEdge(axis, low_or_up) ) {
        if (axis == AXIS::x) {
            for(int k = 1; k < cz_with_glue - 1; ++k){
                for(int j = 1; j < cy_with_glue - 1; ++j){
                    phi[boundary_index][j][k] = 0.0;
                }
            }
        } else if (axis == AXIS::y) {
            for(int k = 1; k < cz_with_glue - 1; ++k){
                for(int i = 1; i < cx_with_glue - 1; ++i){
                    phi[i][boundary_index][k] = 0.0;
                }
            }
        } else {
            for(int j = 1; j < cy_with_glue - 1; ++j){
                for(int i = 1; i < cx_with_glue - 1; ++i){
                    phi[i][j][boundary_index] = 0.0;
                }
            }
        }
    }
}

//! 電位分布の残差の最大ノルムを返す
double Field::checkPhiResidual() {
    double residual = 0.0;
    double rho_max = 0.0;
    double normalized_eps = Utils::Normalizer::normalizeEpsilon(eps0);

    const int cx_with_glue = phi.shape()[0];
    const int cy_with_glue = phi.shape()[1];
    const int cz_with_glue = phi.shape()[2];

    for(int k = 1; k < cz_with_glue - 1; ++k){
        if((k != 1 || !Environment::onLowZedge) && (k != cz_with_glue - 2 || !Environment::onHighZedge)) {
            for(int j = 1; j < cy_with_glue - 1; ++j){
                if((j != 1 || !Environment::onLowYedge) && (j != cy_with_glue - 2 || !Environment::onHighYedge)) {
                    for(int i = 1; i < cx_with_glue - 1; ++i){
                        if((i != 1 || !Environment::onLowXedge) && (i != cx_with_glue - 2 || !Environment::onHighXedge)) {
                            double tmp_res = (phi[i-1][j][k] + phi[i+1][j][k] + phi[i][j-1][k] + phi[i][j+1][k] + phi[i][j][k-1] + phi[i][j][k+1] - 6.0*phi[i][j][k]) + rho[i][j][k]/normalized_eps;
                            residual = std::max(residual, fabs(tmp_res));
                            rho_max = std::max(rho_max, fabs(rho[i][j][k]));
                        }
                    }
                }
            }
        }
    }

    if(MPIw::Environment::numprocs > 1) {
        residual = MPIw::Environment::Comms["world"]->max(residual);
        rho_max =  MPIw::Environment::Comms["world"]->max(rho_max);
    }

    return residual/(rho_max/normalized_eps);
}

//! @brief 差分法で電場を更新する
//! e = - (p_+1 - p_+0)/dx
void Field::updateEfield(const double dx) {
    const int cx_with_glue = ex.shape()[0] + 1; // nx + 2
    const int cy_with_glue = ey.shape()[1] + 1;
    const int cz_with_glue = ez.shape()[2] + 1;
    const double dxm = 1.0/dx;

    //! @note:隣と通信しなくてもいい
    //! phiが通信してあるため、端の要素を通信なしで計算可能
    for(int i = 0; i < cx_with_glue; ++i){
        for(int j = 0; j < cy_with_glue; ++j){
            for(int k = 0; k < cz_with_glue; ++k){
                //! 各方向には1つ少ないのでcx-1まで
                if(i < cx_with_glue - 1) ex[i][j][k] = (phi[i][j][k] - phi[i + 1][j][k]) * dxm;
                if(j < cy_with_glue - 1) ey[i][j][k] = (phi[i][j][k] - phi[i][j + 1][k]) * dxm;
                if(k < cz_with_glue - 1) ez[i][j][k] = (phi[i][j][k] - phi[i][j][k + 1]) * dxm;
            }
        }
    }

    this->updateReferenceEfield();
}

//! Reference Efield (ノード上で定義される電場) を更新する
void Field::updateReferenceEfield() {
    const int cx_with_glue = ex.shape()[0] + 1; // nx + 2
    const int cy_with_glue = ey.shape()[1] + 1;
    const int cz_with_glue = ez.shape()[2] + 1;

    //! reference 更新
    for(int i = 1; i < cx_with_glue - 1; ++i){
        for(int j = 1; j < cy_with_glue - 1; ++j){
            for(int k = 1; k < cz_with_glue - 1; ++k){
                exref[i][j][k] = 0.5 * (ex[i-1][j][k] + ex[i][j][k]);
                eyref[i][j][k] = 0.5 * (ey[i][j-1][k] + ey[i][j][k]);
                ezref[i][j][k] = 0.5 * (ez[i][j][k-1] + ez[i][j][k]);
            }
        }
    }

    //! 外側境界の条件設定
    if (Environment::isBoundary(AXIS::x, AXIS_SIDE::low)) {
        for(int j = 1; j < cy_with_glue - 1; ++j){
            for(int k = 1; k < cz_with_glue - 1; ++k){
                exref[0][j][k] = 0.0;

                // Boundary である場合、更新 (これ正しい??)
                exref[1][j][k] = 0.5 * (phi[3][j][k] - 4.0 * phi[2][j][k] + 3.0 * phi[1][j][k]);
            }
        }
    }

    if (Environment::isBoundary(AXIS::x, AXIS_SIDE::up)) {
        for(int j = 1; j < cy_with_glue - 1; ++j){
            for(int k = 1; k < cz_with_glue - 1; ++k){
                exref[cx_with_glue - 1][j][k] = 0.0;
                exref[cx_with_glue - 2][j][k] = 0.5 * (-phi[cx_with_glue - 4][j][k] + 4.0 * phi[cx_with_glue - 3][j][k] - 3.0 * phi[cx_with_glue - 2][j][k]);
            }
        }
    }

    if (Environment::isBoundary(AXIS::y, AXIS_SIDE::low)) {
        for(int i = 1; i < cx_with_glue - 1; ++i){
            for(int k = 1; k < cz_with_glue - 1; ++k){
                eyref[i][0][k] = 0.0;
                eyref[i][1][k] = 0.5 * (phi[i][3][k] - 4.0 * phi[i][2][k] + 3.0 * phi[i][1][k]);
            }
        }
    }

    if (Environment::isBoundary(AXIS::y, AXIS_SIDE::up)) {
        for(int i = 1; i < cx_with_glue - 1; ++i){
            for(int k = 1; k < cz_with_glue - 1; ++k){
                eyref[i][cy_with_glue - 1][k] = 0.0;
                eyref[i][cy_with_glue - 2][k] = 0.5 * (-phi[i][cy_with_glue - 4][k] + 4.0 * phi[i][cy_with_glue - 3][k] - 3.0 * phi[i][cy_with_glue - 2][k]);
            }
        }
    }

    if (Environment::isBoundary(AXIS::z, AXIS_SIDE::low)) {
        for(int i = 1; i < cx_with_glue - 1; ++i){
            for(int j = 1; j < cy_with_glue - 1; ++j){
                ezref[i][j][0] = 0.0;
                ezref[i][j][1] = 0.5 * (phi[i][j][3] - 4.0 * phi[i][j][2] + 3.0 * phi[i][j][1]);
            }
        }
    }

    if (Environment::isBoundary(AXIS::z, AXIS_SIDE::up)) {
        for(int i = 1; i < cx_with_glue - 1; ++i){
            for(int j = 1; j < cy_with_glue - 1; ++j){
                ezref[i][j][cz_with_glue - 1] = 0.0;
                ezref[i][j][cz_with_glue - 2] = 0.5 * (-phi[i][j][cz_with_glue - 4] + 4.0 * phi[i][j][cz_with_glue - 3] - 3.0 * phi[i][j][cz_with_glue - 2]);
            }
        }
    }

    MPIw::Environment::sendRecvField(exref);
    MPIw::Environment::sendRecvField(eyref);
    MPIw::Environment::sendRecvField(ezref);
}

//! FDTD法ベースで電場を更新する
void Field::updateEfieldFDTD(const double dx) {
    const int cx_with_glue = ex.shape()[0] + 1; // nx + 2
    const int cy_with_glue = ey.shape()[1] + 1;
    const int cz_with_glue = ez.shape()[2] + 1;
    const double dt = Utils::Normalizer::normalizeTime(Environment::dt);
    const double dt_per_eps0 = dt / Utils::Normalizer::normalizeEpsilon(eps0);
    const double dt_per_mu0_eps0_dx = dt_per_eps0 / (Utils::Normalizer::normalizeMu(mu0) * dx);

    for(int i = 1; i < cx_with_glue - 1; ++i){
        for(int j = 1; j < cy_with_glue - 1; ++j){
            for(int k = 1; k < cz_with_glue - 1; ++k){
                //! 各方向には1つ少ないのでcx-1まで
                if(i < cx_with_glue - 2) {
                    ex[i][j][k] = ex[i][j][k] - jx[i][j][k] * dt_per_eps0 +
                        dt_per_mu0_eps0_dx * (bz[i][j][k] - bz[i][j - 1][k] - by[i][j][k] + by[i][j][k - 1]);
                }
                if(j < cy_with_glue - 2) {
                    ey[i][j][k] = ey[i][j][k] - jy[i][j][k] * dt_per_eps0 +
                        dt_per_mu0_eps0_dx * (bx[i][j][k] - bx[i][j][k - 1] - bz[i][j][k] + bz[i - 1][j][k]);
                }
                if(k < cz_with_glue - 2) {
                    ez[i][j][k] = ez[i][j][k] - jz[i][j][k] * dt_per_eps0 +
                        dt_per_mu0_eps0_dx * (by[i][j][k] - by[i - 1][j][k] - bx[i][j][k] + bx[i][j - 1][k]);
                }
            }
        }
    }

    // FDTDの場合は通信が必要になる
    MPIw::Environment::sendRecvField(ex);
    MPIw::Environment::sendRecvField(ey);
    MPIw::Environment::sendRecvField(ez);

    //! 境界条件設定
    this->setDampingBoundaryOnEfield();

    //! phi correction?

    //! Reference 更新
    this->updateReferenceEfield();
}

void Field::setDampingBoundaryOnEfield(void) {
    const int cx_with_glue = ex.shape()[0] + 1;
    const int cy_with_glue = ey.shape()[1] + 1;
    const int cz_with_glue = ez.shape()[2] + 1;

    if (!Environment::isPeriodic(AXIS::x, AXIS_SIDE::low)) {
        for(int k = 0; k < cz_with_glue; ++k){
            for(int j = 0; j < cy_with_glue; ++j){
                ex[0][j][k] = 0.0; // glue cell
                ex[1][j][k] = 0.0;
            }
        }
    }

    if (!Environment::isPeriodic(AXIS::x, AXIS_SIDE::up)) {
        for(int k = 0; k < cz_with_glue; ++k){
            for(int j = 0; j < cy_with_glue; ++j){
                //! ex は x方向に 1 小さいので、cx_with_glue - 2 が glue cell になる
                ex[cx_with_glue - 2][j][k] = 0.0; // glue cell
                ex[cx_with_glue - 3][j][k] = 0.0;
            }
        }
    }

    if (!Environment::isPeriodic(AXIS::y, AXIS_SIDE::low)) {
        for(int i = 0; i < cx_with_glue; ++i){
            for(int k = 0; k < cz_with_glue; ++k){
                ey[i][0][k] = 0.0; // glue cell
                ey[i][1][k] = 0.0;
            }
        }
    }

    if (!Environment::isPeriodic(AXIS::y, AXIS_SIDE::up)) {
        for(int i = 0; i < cx_with_glue; ++i){
            for(int k = 0; k < cz_with_glue; ++k){
                ey[i][cy_with_glue - 2][k] = 0.0; // glue cell
                ey[i][cy_with_glue - 3][k] = 0.0;
            }
        }
    }

        if (!Environment::isPeriodic(AXIS::z, AXIS_SIDE::low)) {
        for(int i = 0; i < cx_with_glue; ++i){
            for(int j = 0; j < cy_with_glue; ++j){
                ez[i][j][0] = 0.0; // glue cell
                ez[i][j][1] = 0.0;
            }
        }
    }

    if (!Environment::isPeriodic(AXIS::z, AXIS_SIDE::up)) {
        for(int i = 0; i < cx_with_glue; ++i){
            for(int j = 0; j < cy_with_glue; ++j){
                ez[i][j][cz_with_glue - 2] = 0.0; // glue cell
                ez[i][j][cz_with_glue - 3] = 0.0;
            }
        }
    }
}

//! @brief 磁場を更新する(FDTD)
//!
void Field::updateBfield(const double dx, const int nx, const int ny, const int nz) {
    const double dt_per_dx = Utils::Normalizer::normalizeTime(Environment::dt) / dx;

    // const double epsilon_r = 1.0; //! 比誘電率
    // const double sigma = 1.0; //! 導電率 (各Faceでの)
    // const double mu_r = 1.0; //! 透磁率 (各Faceでの)
    // const double sigma_m = 0.0; //! 導磁率?

    // 磁束密度更新時の係数
    // const double d1 = mu_r/(mu_r + sigma_m * dt / mu0);
    // const double d2 = dt/(mu_r + sigma_m * dt / mu0);

    const int cx_with_glue = nx + 1;
    const int cy_with_glue = nx + 1;
    const int cz_with_glue = nx + 1;

    //! 0とcy + 1, 0とcz + 1はglueなので更新しなくてよい
    for(int i = 1; i < cx_with_glue; ++i){
        for(int j = 1; j < cy_with_glue; ++j){
            for(int k = 1; k < cz_with_glue; ++k){
                if (j != cy_with_glue - 1 && k != cz_with_glue - 1) {
                    bx[i][j][k] = bx[i][j][k] - dt_per_dx * (ez[i][j - 1][k] - ez[i][j][k] - ey[i][j][k + 1] + ey[i][j][k]);
                }

                if (i != cx_with_glue - 1 && k != cz_with_glue - 1) {
                    by[i][j][k] = by[i][j][k] - dt_per_dx * (ex[i][j][k + 1] - ex[i][j][k] - ez[i + 1][j][k] + ez[i][j][k]);
                }

                if (i != cx_with_glue - 1 && j != cy_with_glue - 1) {
                    bz[i][j][k] = bz[i][j][k] - dt_per_dx * (ey[i + 1][j][k] - ey[i][j][k] - ex[i][j + 1][k] + ex[i][j][k]);
                }
            }
        }
    }
    MPIw::Environment::sendRecvField(bx);
    MPIw::Environment::sendRecvField(by);
    MPIw::Environment::sendRecvField(bz);

    //! reference 更新
    for(int i = 1; i < cx_with_glue - 1; ++i){
        for(int j = 1; j < cy_with_glue - 1; ++j){
            for(int k = 1; k < cz_with_glue - 1; ++k){
                bxref[i][j][k] = 0.25 * (bx[i][j][k] + bx[i][j-1][k] + bx[i][j][k-1] + bx[i][j-1][k-1]);
                byref[i][j][k] = 0.25 * (by[i][j][k] + by[i-1][j][k] + by[i][j][k-1] + by[i-1][j][k-1]);
                bzref[i][j][k] = 0.25 * (bz[i][j][k] + bz[i-1][j][k] + bz[i][j-1][k] + bz[i-1][j-1][k]);
            }
        }
    }

    MPIw::Environment::sendRecvField(bxref);
    MPIw::Environment::sendRecvField(byref);
    MPIw::Environment::sendRecvField(bzref);
}

double Field::getEfieldEnergy(void) const {
    const int cx_with_glue = ey.shape()[0];
    const int cy_with_glue = ex.shape()[1];
    const int cz_with_glue = ex.shape()[2];

    double energy = 0.0;

    for(int i = 1; i < cx_with_glue - 1; ++i){
        for(int j = 1; j < cy_with_glue - 1; ++j){
            for(int k = 1; k < cz_with_glue - 1; ++k){
                //! 各点のエネルギーを計算する(Yee格子内のエネルギーの計算方法は?)
                if (i < cx_with_glue - 2) energy += pow(ex[i][j][k], 2);
                if (j < cy_with_glue - 2) energy += pow(ey[i][j][k], 2);
                if (k < cz_with_glue - 2) energy += pow(ez[i][j][k], 2);
            }
        }
    }

    return 0.5 * Utils::Normalizer::normalizeEpsilon(eps0) * energy;
}

double Field::getBfieldEnergy(void) const {
    const int cx_with_glue = bx.shape()[0];
    const int cy_with_glue = by.shape()[1];
    const int cz_with_glue = bz.shape()[2];

    double energy = 0.0;

    for(int i = 1; i < cx_with_glue - 1; ++i){
        for(int j = 1; j < cy_with_glue - 1; ++j){
            for(int k = 1; k < cz_with_glue - 1; ++k){
                //! 各点のエネルギーを計算する(Yee格子内のエネルギーの計算方法は?)
                if ( (j < cy_with_glue - 2) && (k < cz_with_glue - 2) ) energy += pow(bx[i][j][k], 2);
                if ( (i < cx_with_glue - 2) && (k < cz_with_glue - 2) ) energy += pow(by[i][j][k], 2);
                if ( (i < cx_with_glue - 2) && (j < cy_with_glue - 2) ) energy += pow(bz[i][j][k], 2);
            }
        }
    }

    return 0.5 * energy / Utils::Normalizer::normalizeMu(mu0);
}

void Field::checkChargeConservation(const tdArray& old_rho, const double dt, const double dx) const {
    double residual1 = 0.0, residual2 = 0.0;
    double j1 = 0.0, j2 = 0.0, j3 = 0.0;
    double residual = 0.0;
    const int cx_with_glue = rho.shape()[0];
    const int cy_with_glue = rho.shape()[1];
    const int cz_with_glue = rho.shape()[2];

    for(int i = 1; i < cx_with_glue - 1; ++i) {
        for(int j = 1; j < cy_with_glue - 1; ++j) {
            for(int k = 1; k < cz_with_glue - 1; ++k) {
                residual1 = 0.0, residual2 = 0.0;
                j1 = 0.0, j2 = 0.0, j3 = 0.0;

                residual1 = (rho[i][j][k])/dt;
                residual2 = (old_rho[i][j][k])/dt;
                j1 = (jx[i][j][k] - jx[i - 1][j][k])/dx;
                j2 = (jy[i][j][k] - jy[i][j - 1][k])/dx;
                j3 = (jz[i][j][k] - jz[i][j][k - 1])/dx;

                // cout << "[TEST] Charge conservation residual1 = " << residual1 << endl;
                // cout << "[TEST] Charge conservation residual2 = " << residual2 << endl;
                // cout << "[TEST] Charge conservation j1 = " << j1 << endl;
                // cout << "[TEST] Charge conservation j2 = " << j2 << endl;
                // cout << "[TEST] Charge conservation j3 = " << j3 << endl;
                // cout << "[TEST] Charge conservation residual = " << (residual1 - residual2 + j1 + j2 + j3) << endl;

                residual += residual1 - residual2 + j1 + j2 + j3;
            }
        }
    }

    cout << "[TEST] Charge conservation residual sum = " << residual << endl;
}
