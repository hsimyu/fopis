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
    const int one_freq = static_cast<int>(1.0 / Normalizer::normalizeFrequency(real_freq));

    if (Environment::isRootNode) {
        cout << "[NOTICE] " << one_freq << " step で 1周期です." << endl;
    }

    const double freq = 2.0 * M_PI * Normalizer::normalizeFrequency(real_freq); // Hz
    const double J0 = 100.0;

    const size_t half_x = cx_with_glue / 4;
    const size_t half_y = cy_with_glue / 2;

    if (Environment::isRootNode && Environment::timestep < one_freq) {
        for (size_t k = 0; k < cz_with_glue - 1; ++k) {
            jz[half_x][half_y][k] = J0 * std::sin(freq * now);
        }
    }
}

double Field::getEfieldEnergy(void) const {
    const size_t cx_with_glue = exref.shape()[0];
    const size_t cy_with_glue = exref.shape()[1];
    const size_t cz_with_glue = exref.shape()[2];

    double energy = 0.0;

    for(size_t i = 1; i < cx_with_glue - 1; ++i){
        for(size_t j = 1; j < cy_with_glue - 1; ++j){
            for(size_t k = 1; k < cz_with_glue - 1; ++k){
                //! 各エッジ上のエネルギーを計算する
                if (i < cx_with_glue - 2) energy += pow(ex[i][j][k], 2);
                if (j < cy_with_glue - 2) energy += pow(ey[i][j][k], 2);
                if (k < cz_with_glue - 2) energy += pow(ez[i][j][k], 2);
            }
        }
    }

    return 0.5 * Normalizer::eps0 * energy;
}

double Field::getBfieldEnergy(void) const {
    const size_t cx_with_glue = bxref.shape()[0];
    const size_t cy_with_glue = byref.shape()[1];
    const size_t cz_with_glue = bzref.shape()[2];

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
