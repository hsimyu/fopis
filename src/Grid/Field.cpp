#include "global.hpp"
#include "environment.hpp"
#include "field.hpp"
#include "utils.hpp"
#include "mpiw.hpp"
#include "normalizer.hpp"
#include <cmath>
#include <algorithm>

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
    const size_t max_index = phi.shape()[ Utils::getAxisIndex(axis) ];

    const size_t boundary_index = (low_or_up == AXIS_SIDE::low) ? 1 : (max_index - 2);
    const size_t neighbor_element_index = (low_or_up == AXIS_SIDE::low) ? 2 : (max_index - 3);

    const size_t cx_with_glue = phi.shape()[0];
    const size_t cy_with_glue = phi.shape()[1];
    const size_t cz_with_glue = phi.shape()[2];

    if ( Environment::isOnEdge(axis, low_or_up) ) {
        if (axis == AXIS::x) {
            for(size_t k = 1; k < cz_with_glue - 1; ++k){
                for(size_t j = 1; j < cy_with_glue - 1; ++j){
                    phi[boundary_index][j][k] = phi[neighbor_element_index][j][k];
                }
            }
        } else if (axis == AXIS::y) {
            for(size_t k = 1; k < cz_with_glue - 1; ++k){
                for(size_t i = 1; i < cx_with_glue - 1; ++i){
                    phi[i][boundary_index][k] = phi[i][neighbor_element_index][k];
                }
            }
        } else {
            for(size_t j = 1; j < cy_with_glue - 1; ++j){
                for(size_t i = 1; i < cx_with_glue - 1; ++i){
                    phi[i][j][boundary_index] = phi[i][j][neighbor_element_index];
                }
            }
        }
    }
}

void Field::setDirichletPhi(const AXIS axis, const AXIS_SIDE low_or_up) {
    const size_t max_index = phi.shape()[ Utils::getAxisIndex(axis) ];
    const size_t boundary_index = (low_or_up == AXIS_SIDE::low) ? 1 : (max_index - 2);

    const size_t cx_with_glue = phi.shape()[0];
    const size_t cy_with_glue = phi.shape()[1];
    const size_t cz_with_glue = phi.shape()[2];

    if ( Environment::isOnEdge(axis, low_or_up) ) {
        if (axis == AXIS::x) {
            for(size_t k = 1; k < cz_with_glue - 1; ++k){
                for(size_t j = 1; j < cy_with_glue - 1; ++j){
                    phi[boundary_index][j][k] = 0.0;
                }
            }
        } else if (axis == AXIS::y) {
            for(size_t k = 1; k < cz_with_glue - 1; ++k){
                for(size_t i = 1; i < cx_with_glue - 1; ++i){
                    phi[i][boundary_index][k] = 0.0;
                }
            }
        } else {
            for(size_t j = 1; j < cy_with_glue - 1; ++j){
                for(size_t i = 1; i < cx_with_glue - 1; ++i){
                    phi[i][j][boundary_index] = 0.0;
                }
            }
        }
    }
}

//! FDTD法ベースで電場を更新する
void Field::updateEfieldFDTD(const double dx, const double dt) {
    const size_t cx_with_glue = ex.shape()[0] + 1; // nx + 2
    const size_t cy_with_glue = ey.shape()[1] + 1;
    const size_t cz_with_glue = ez.shape()[2] + 1;
    const double dt_per_eps0 = dt / Normalizer::eps0;
    const double dt_per_mu0_eps0_dx = dt_per_eps0 / (Normalizer::mu0 * dx);

    for(size_t i = 1; i < cx_with_glue - 1; ++i){
        for(size_t j = 1; j < cy_with_glue - 1; ++j){
            for(size_t k = 1; k < cz_with_glue - 1; ++k){
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
    // this->updateReferenceEfield();
}

void Field::setDampingBoundaryOnEfield(void) {
    const size_t cx_with_glue = ex.shape()[0] + 1;
    const size_t cy_with_glue = ey.shape()[1] + 1;
    const size_t cz_with_glue = ez.shape()[2] + 1;

    if (!Environment::isNotBoundary(AXIS::x, AXIS_SIDE::low)) {
        for(size_t k = 0; k < cz_with_glue; ++k){
            for(size_t j = 0; j < cy_with_glue; ++j){
                ex[0][j][k] = 0.0; // glue cell
                ex[1][j][k] = 0.0;
            }
        }
    }

    if (!Environment::isNotBoundary(AXIS::x, AXIS_SIDE::up)) {
        for(size_t k = 0; k < cz_with_glue; ++k){
            for(size_t j = 0; j < cy_with_glue; ++j){
                //! ex は x方向に 1 小さいので、cx_with_glue - 2 が glue cell になる
                ex[cx_with_glue - 2][j][k] = 0.0; // glue cell
                ex[cx_with_glue - 3][j][k] = 0.0;
            }
        }
    }

    if (!Environment::isNotBoundary(AXIS::y, AXIS_SIDE::low)) {
        for(size_t i = 0; i < cx_with_glue; ++i){
            for(size_t k = 0; k < cz_with_glue; ++k){
                ey[i][0][k] = 0.0; // glue cell
                ey[i][1][k] = 0.0;
            }
        }
    }

    if (!Environment::isNotBoundary(AXIS::y, AXIS_SIDE::up)) {
        for(size_t i = 0; i < cx_with_glue; ++i){
            for(size_t k = 0; k < cz_with_glue; ++k){
                ey[i][cy_with_glue - 2][k] = 0.0; // glue cell
                ey[i][cy_with_glue - 3][k] = 0.0;
            }
        }
    }

        if (!Environment::isNotBoundary(AXIS::z, AXIS_SIDE::low)) {
        for(size_t i = 0; i < cx_with_glue; ++i){
            for(size_t j = 0; j < cy_with_glue; ++j){
                ez[i][j][0] = 0.0; // glue cell
                ez[i][j][1] = 0.0;
            }
        }
    }

    if (!Environment::isNotBoundary(AXIS::z, AXIS_SIDE::up)) {
        for(size_t i = 0; i < cx_with_glue; ++i){
            for(size_t j = 0; j < cy_with_glue; ++j){
                ez[i][j][cz_with_glue - 2] = 0.0; // glue cell
                ez[i][j][cz_with_glue - 3] = 0.0;
            }
        }
    }
}

//! @brief 磁場を更新する(FDTD)
//!
void Field::updateBfield(const double dx, const int nx, const int ny, const int nz, const double dt) {
    const double dt_per_dx = dt / dx;

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

void Field::initializeCurrent(const double dt) {
    Utils::initialize3DArray(jx);
    Utils::initialize3DArray(jy);
    Utils::initialize3DArray(jz);

    const size_t cx_with_glue = jx.shape()[0] + 1; // Edge 要素は各方向に1少ないので +1 する
    const size_t cy_with_glue = jy.shape()[1] + 1;
    const size_t cz_with_glue = jz.shape()[2] + 1;

    //! 背景電流などがある場合にはここで設定する

    //! Jz 方向に振動する電流
    const auto now = dt * static_cast<double>(Environment::timestep);
    const double real_freq = 5e7; // Hz
    // cout << "[NOTICE] " << 1.0 / Normalizer::normalizeFrequency(real_freq) << " step で 1周期です." << endl;
    const double freq = 2.0 * M_PI * Normalizer::normalizeFrequency(real_freq); // Hz
    const double J0 = 10.0;
    const size_t half_x = cx_with_glue / 2;
    const size_t half_y = cy_with_glue / 2;
    for (size_t k = 0; k < cz_with_glue - 1; ++k) {
        jz[half_x][half_y][k] = J0 * std::sin(freq * now);
    }
}

double Field::getEfieldEnergy(void) const {
    const size_t cx_with_glue = ey.shape()[0];
    const size_t cy_with_glue = ex.shape()[1];
    const size_t cz_with_glue = ex.shape()[2];

    double energy = 0.0;

    for(size_t i = 1; i < cx_with_glue - 1; ++i){
        for(size_t j = 1; j < cy_with_glue - 1; ++j){
            for(size_t k = 1; k < cz_with_glue - 1; ++k){
                //! 各点のエネルギーを計算する(Yee格子内のエネルギーの計算方法は?)
                if (i < cx_with_glue - 2) energy += pow(ex[i][j][k], 2);
                if (j < cy_with_glue - 2) energy += pow(ey[i][j][k], 2);
                if (k < cz_with_glue - 2) energy += pow(ez[i][j][k], 2);
            }
        }
    }

    return 0.5 * Normalizer::eps0 * energy;
}

double Field::getBfieldEnergy(void) const {
    const size_t cx_with_glue = bx.shape()[0];
    const size_t cy_with_glue = by.shape()[1];
    const size_t cz_with_glue = bz.shape()[2];

    double energy = 0.0;

    for(size_t i = 1; i < cx_with_glue - 1; ++i){
        for(size_t j = 1; j < cy_with_glue - 1; ++j){
            for(size_t k = 1; k < cz_with_glue - 1; ++k){
                //! 各点のエネルギーを計算する(Yee格子内のエネルギーの計算方法は?)
                if ( (j < cy_with_glue - 2) && (k < cz_with_glue - 2) ) energy += pow(bx[i][j][k], 2);
                if ( (i < cx_with_glue - 2) && (k < cz_with_glue - 2) ) energy += pow(by[i][j][k], 2);
                if ( (i < cx_with_glue - 2) && (j < cy_with_glue - 2) ) energy += pow(bz[i][j][k], 2);
            }
        }
    }

    return 0.5 * energy / Normalizer::mu0;
}

void Field::checkChargeConservation(const RhoArray& old_rho, const double dt, const double dx) const {
    double residual1 = 0.0, residual2 = 0.0;
    double j1 = 0.0, j2 = 0.0, j3 = 0.0;
    double residual = 0.0;
    const size_t cx_with_glue = rho[0].shape()[0];
    const size_t cy_with_glue = rho[0].shape()[1];
    const size_t cz_with_glue = rho[0].shape()[2];

    for(size_t i = 1; i < cx_with_glue - 1; ++i) {
        for(size_t j = 1; j < cy_with_glue - 1; ++j) {
            for(size_t k = 1; k < cz_with_glue - 1; ++k) {
                residual1 = 0.0, residual2 = 0.0;
                j1 = 0.0, j2 = 0.0, j3 = 0.0;

                residual1 = (rho[0][i][j][k])/dt;
                residual2 = (old_rho[0][i][j][k])/dt;
                j1 = (jx[i][j][k] - jx[i - 1][j][k])/dx;
                j2 = (jy[i][j][k] - jy[i][j - 1][k])/dx;
                j3 = (jz[i][j][k] - jz[i][j][k - 1])/dx;

                residual += residual1 - residual2 + j1 + j2 + j3;
            }
        }
    }

    residual = MPIw::Environment::Comms["world"].sum(residual);
    if (Environment::isRootNode) cout << "[TEST] Charge conservation residual sum = " << residual << endl;
}
