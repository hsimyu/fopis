#include "grid.hpp"
#include "normalizer.hpp"
#include "field.hpp"
#include "utils.hpp"

void RootGrid::solvePoisson(void) {
    const auto PRE_LOOP_NUM = Environment::getOptions().getMaximumPoissonPreLoop();
    const auto POST_LOOP_NUM = Environment::getOptions().getMaximumPoissonPostLoop();

    if (this->getChildrenLength() > 0) {
        this->solvePoissonPSOR(PRE_LOOP_NUM);
        this->restrictPhiToChildrenBoundary();

        for(auto& child : children) {
            child->solvePoissonFromParent();
        }
    }

    this->solvePoissonPSOR(POST_LOOP_NUM);

    if (this->getChildrenLength() > 0) {
        this->correctChildrenPhi();

        for(auto& child : children) {
            child->solvePoisson();
        }
    }
}

void RootGrid::solvePoissonPSOR(const int loopnum) {
    auto& phi = field->getPhi();
    auto& poisson_error = field->getPoissonError();
    auto& rho = field->getRho();

    const double omega = 2.0/(1.0 + sin(M_PI/(phi.shape()[0] - 2))); // spectral radius
    const double rho_coeff = pow(dx, 2) / Normalizer::eps0;

    const size_t cx_with_glue = phi.shape()[0];
    const size_t cy_with_glue = phi.shape()[1];
    const size_t cz_with_glue = phi.shape()[2];

    constexpr double required_error = 1.0e-7;

    const bool is_not_boundary[6] = {
        Environment::isNotBoundary(AXIS::x, AXIS_SIDE::low),
        Environment::isNotBoundary(AXIS::x, AXIS_SIDE::up),
        Environment::isNotBoundary(AXIS::y, AXIS_SIDE::low),
        Environment::isNotBoundary(AXIS::y, AXIS_SIDE::up),
        Environment::isNotBoundary(AXIS::z, AXIS_SIDE::low),
        Environment::isNotBoundary(AXIS::z, AXIS_SIDE::up),
    };

    auto time_counter = Utils::TimeCounter::getInstance();

    time_counter->begin("solvePoisson/updatePoissonErrorPre");
    poisson_error = phi;

    size_t loop = 1;
    for(; loop <= loopnum; ++loop) {
        time_counter->switchTo("solvePoisson/mainLoop");
        //! is_not_boundaryは各方向の境界について
        //! 「その境界は計算空間の境界でない」か「その境界が計算空間の境界であり、周期境界である」場合にtrueとなるため
        //! 各方向の端要素の計算をIterationで計算する場合にチェックする必要がある

        #pragma omp parallel
        {
            //! 奇数グリッド更新
            #pragma omp for
            for(int i = 1; i < cx_with_glue - 1; i += 2){
                if((i != 1 || is_not_boundary[0]) && (i != cx_with_glue - 2 || is_not_boundary[1])) {
                    for(int j = 1; j < cy_with_glue - 1; ++j){
                        if((j != 1 || is_not_boundary[2]) && (j != cy_with_glue - 2 || is_not_boundary[3])) {
                        for(int k = 1; k < cz_with_glue - 1; ++k){
                            if((k != 1 || is_not_boundary[4]) && (k != cz_with_glue - 2 || is_not_boundary[5])) {
                                    phi[i][j][k] = (1.0 - omega) * phi[i][j][k] + omega*(phi[i+1][j][k] + phi[i-1][j][k] + phi[i][j+1][k] + phi[i][j-1][k] + phi[i][j][k+1] + phi[i][j][k-1] + rho_coeff * rho[0][i][j][k])/6.0;
                                }
                            }
                        }
                    }
                }
            }

            //! 偶数グリッド更新
            #pragma omp for
            for(int i = 2; i < cx_with_glue - 1; i += 2){
                if((i != 1 || is_not_boundary[0]) && (i != cx_with_glue - 2 || is_not_boundary[1])) {
                    for(int j = 1; j < cy_with_glue - 1; ++j){
                        if((j != 1 || is_not_boundary[2]) && (j != cy_with_glue - 2 || is_not_boundary[3])) {
                            for(int k = 1; k < cz_with_glue - 1; ++k){
                                if((k != 1 || is_not_boundary[4]) && (k != cz_with_glue - 2 || is_not_boundary[5])) {
                                    phi[i][j][k] = (1.0 - omega) * phi[i][j][k] + omega*(phi[i+1][j][k] + phi[i-1][j][k] + phi[i][j+1][k] + phi[i][j-1][k] + phi[i][j][k+1] + phi[i][j][k-1] + rho_coeff * rho[0][i][j][k])/6.0;
                                }
                            }
                        }
                    }
                }
            }
        }

        time_counter->switchTo("solvePoisson/sendPhi");
        //! @note: 実際には部分部分をSORで計算して送信というのを繰り返す方が収束効率がよい
        MPIw::Environment::sendRecvField(phi);

        time_counter->switchTo("solvePoisson/checkResidual");
        if ( (loop % 10 == 0) && (this->checkPhiResidual() < required_error) ) {
            break;
        }
    }

    if (Environment::isRootNode) {
        cout << "[INFO] solve poisson: performed " << (loop - 1) << " iterations." << endl;
    }

    time_counter->switchTo("solvePoisson/setBoundaryCondition");
    field->setBoundaryConditionPhi();

    //! 全グリッド上のエラーを更新
    time_counter->switchTo("solvePoisson/updatePoissonErrorPost");
    #pragma omp parallel for
    for(int i = 1; i < cx_with_glue - 1; ++i){
        for(int j = 1; j < cy_with_glue - 1; ++j){
            for(int k = 1; k < cz_with_glue - 1; ++k){
                poisson_error[i][j][k] = phi[i][j][k] - poisson_error[i][j][k];
            }
        }
    }
    time_counter->end();
}

//! 電位分布の残差の最大ノルムを返す
double RootGrid::checkPhiResidual() {
    double residual = 0.0;
    const double normalized_eps = Normalizer::eps0;

    const bool is_not_boundary[6] = {
        Environment::isNotBoundary(AXIS::x, AXIS_SIDE::low),
        Environment::isNotBoundary(AXIS::x, AXIS_SIDE::up),
        Environment::isNotBoundary(AXIS::y, AXIS_SIDE::low),
        Environment::isNotBoundary(AXIS::y, AXIS_SIDE::up),
        Environment::isNotBoundary(AXIS::z, AXIS_SIDE::low),
        Environment::isNotBoundary(AXIS::z, AXIS_SIDE::up),
    };

    auto& phi = field->getPhi();
    auto& rho = field->getRho();
    auto& poisson_residual = field->getPoissonResidual();

    const size_t cx_with_glue = phi.shape()[0];
    const size_t cy_with_glue = phi.shape()[1];
    const size_t cz_with_glue = phi.shape()[2];

    #pragma omp parallel for reduction(max: residual)
    for(size_t i = 1; i < cx_with_glue - 1; ++i){
        if((i != 1 || is_not_boundary[0]) && (i != cx_with_glue - 2 || is_not_boundary[1])) {
            for(size_t j = 1; j < cy_with_glue - 1; ++j){
                if((j != 1 || is_not_boundary[2]) && (j != cy_with_glue - 2 || is_not_boundary[3])) {
                    for(size_t k = 1; k < cz_with_glue - 1; ++k){
                        if((k != 1 || is_not_boundary[4]) && (k != cz_with_glue - 2 || is_not_boundary[5])) {
                            double source_value = rho[0][i][j][k]/normalized_eps;

                            if (fabs(source_value) > 0.0) {
                                double tmp_res = field->poissonOperator(phi, i, j, k) + source_value;
                                poisson_residual[i][j][k] = tmp_res;

                                if (fabs(tmp_res / source_value) > residual) {
                                    residual = fabs(tmp_res / source_value);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if(MPIw::Environment::numprocs > 1) {
        residual = MPIw::Environment::Comms["world"].max(residual);
    }

    return residual;
}

void RootGrid::solvePoissonCorrection(void) {
    constexpr int POST_LOOP_NUM = 250;
    this->solvePoissonCorrectionPSOR(POST_LOOP_NUM);
}

//! FDTD用にsolvePoissonの代わりに誤差計算を行う
void RootGrid::solvePoissonCorrectionPSOR(const int loopnum) {
    //! 残差を更新
    this->checkPhiResidual();

    auto& residual = field->getPoissonResidual();
    auto& poisson_error = field->getPoissonError();

    const double omega = 2.0/(1.0 + sin(M_PI/(residual.shape()[0] - 2))); // spectral radius
    const double coeff = pow(dx, 2);

    const size_t cx_with_glue = residual.shape()[0];
    const size_t cy_with_glue = residual.shape()[1];
    const size_t cz_with_glue = residual.shape()[2];

    constexpr double required_error = 1.0e-7;

    const bool is_not_boundary[6] = {
        Environment::isNotBoundary(AXIS::x, AXIS_SIDE::low),
        Environment::isNotBoundary(AXIS::x, AXIS_SIDE::up),
        Environment::isNotBoundary(AXIS::y, AXIS_SIDE::low),
        Environment::isNotBoundary(AXIS::y, AXIS_SIDE::up),
        Environment::isNotBoundary(AXIS::z, AXIS_SIDE::low),
        Environment::isNotBoundary(AXIS::z, AXIS_SIDE::up),
    };

    for(int loop = 1; loop <= loopnum; ++loop) {
        #pragma omp parallel shared(poisson_error,residual)
        {
            //! 奇数グリッド更新
            #pragma omp for
            for(size_t k = 1; k < cz_with_glue - 1; k += 2){
                if((k != 1 || is_not_boundary[4]) && (k != cz_with_glue - 2 || is_not_boundary[5])) {
                    for(size_t j = 1; j < cy_with_glue - 1; ++j){
                        if((j != 1 || is_not_boundary[2]) && (j != cy_with_glue - 2 || is_not_boundary[3])) {
                            for(size_t i = 1; i < cx_with_glue - 1; ++i){
                                if((i != 1 || is_not_boundary[0]) && (i != cx_with_glue - 2 || is_not_boundary[1])) {
                                    poisson_error[i][j][k] = (1.0 - omega) * poisson_error[i][j][k] + omega * (poisson_error[i+1][j][k] + poisson_error[i-1][j][k] + poisson_error[i][j+1][k] + poisson_error[i][j-1][k] + poisson_error[i][j][k+1] + poisson_error[i][j][k-1] + coeff * residual[i][j][k])/6.0;
                                }
                            }
                        }
                    }
                }
            }

            //! 偶数グリッド更新
            #pragma omp for
            for(size_t k = 2; k < cz_with_glue - 1; k += 2){
                if((k != 1 || is_not_boundary[4]) && (k != cz_with_glue - 2 || is_not_boundary[5])) {
                    for(size_t j = 1; j < cy_with_glue - 1; ++j){
                        if((j != 1 || is_not_boundary[2]) && (j != cy_with_glue - 2 || is_not_boundary[3])) {
                            for(size_t i = 1; i < cx_with_glue - 1; ++i){
                                if((i != 1 || is_not_boundary[0]) && (i != cx_with_glue - 2 || is_not_boundary[1])) {
                                    poisson_error[i][j][k] = (1.0 - omega) * poisson_error[i][j][k] + omega * (poisson_error[i+1][j][k] + poisson_error[i-1][j][k] + poisson_error[i][j+1][k] + poisson_error[i][j-1][k] + poisson_error[i][j][k+1] + poisson_error[i][j][k-1] + coeff * residual[i][j][k])/6.0;
                                }
                            }
                        }
                    }
                }
            }
        }

        MPIw::Environment::sendRecvField(poisson_error);

        if ( (loop % 10 == 0) && (this->checkPhiCorrectionResidual() < required_error) ) {
            if (Environment::isRootNode) {
                cout << "[INFO] solve poisson correction: performed " << loop << " iterations." << endl;
            }
            break;
        }
    }

    auto& phi = field->getPhi();

    #pragma omp parallel for shared(poisson_error, phi)
    for (int i = 0; i < cx_with_glue; ++i) {
        for (int j = 0; j < cy_with_glue; ++j) {
            for (int k = 0; k < cz_with_glue; ++k) {
                phi[i][j][k] += poisson_error[i][j][k];
            }
        }
    }
}

//! 電位分布の残差の残差の最大ノルムを返す
double RootGrid::checkPhiCorrectionResidual() {
    double residual = 0.0;

    const bool is_not_boundary[6] = {
        Environment::isNotBoundary(AXIS::x, AXIS_SIDE::low),
        Environment::isNotBoundary(AXIS::x, AXIS_SIDE::up),
        Environment::isNotBoundary(AXIS::y, AXIS_SIDE::low),
        Environment::isNotBoundary(AXIS::y, AXIS_SIDE::up),
        Environment::isNotBoundary(AXIS::z, AXIS_SIDE::low),
        Environment::isNotBoundary(AXIS::z, AXIS_SIDE::up),
    };

    auto& poisson_error = field->getPoissonError();
    auto& poisson_residual = field->getPoissonResidual();

    const size_t cx_with_glue = poisson_residual.shape()[0];
    const size_t cy_with_glue = poisson_residual.shape()[1];
    const size_t cz_with_glue = poisson_residual.shape()[2];

    #pragma omp parallel for reduction(max: residual)
    for(size_t k = 1; k < cz_with_glue - 1; ++k){
        if((k != 1 || is_not_boundary[4]) && (k != cz_with_glue - 2 || is_not_boundary[5])) {
            for(size_t j = 1; j < cy_with_glue - 1; ++j){
                if((j != 1 || is_not_boundary[2]) && (j != cy_with_glue - 2 || is_not_boundary[3])) {
                    for(size_t i = 1; i < cx_with_glue - 1; ++i){
                        if((i != 1 || is_not_boundary[0]) && (i != cx_with_glue - 2 || is_not_boundary[1])) {
                            double source_value = poisson_residual[i][j][k];
                            double tmp_res = field->poissonOperator(poisson_error, i, j, k) + source_value;

                            if (fabs(tmp_res / source_value) > residual) {
                                residual = fabs(tmp_res / source_value);
                            }
                        }
                    }
                }
            }
        }
    }

    if (MPIw::Environment::numprocs > 1) {
        residual = MPIw::Environment::Comms["world"].max(residual);
    }
    return residual;
}

//! @brief 差分法で電場を更新する
//! e = - (p_+1 - p_+0)/dx
void RootGrid::updateEfield() {
    auto& ex = field->getEx();
    auto& ey = field->getEy();
    auto& ez = field->getEz();
    auto& phi = field->getPhi();

    const size_t cx_with_glue = ex.shape()[0] + 1; // nx + 2
    const size_t cy_with_glue = ey.shape()[1] + 1;
    const size_t cz_with_glue = ez.shape()[2] + 1;
    const double per_dx = 1.0 / dx;

    //! @note:隣と通信しなくてもいい
    //! phiが通信してあるため、端の要素を通信なしで計算可能
    #pragma omp parallel for shared(ex, ey, ez)
    for(size_t i = 0; i < cx_with_glue; ++i){
        for(size_t j = 0; j < cy_with_glue; ++j){
            for(size_t k = 0; k < cz_with_glue; ++k){
                //! 各方向には1つ少ないのでcx-1まで
                if(i < cx_with_glue - 1) ex[i][j][k] = (phi[i][j][k] - phi[i + 1][j][k]) * per_dx;
                if(j < cy_with_glue - 1) ey[i][j][k] = (phi[i][j][k] - phi[i][j + 1][k]) * per_dx;
                if(k < cz_with_glue - 1) ez[i][j][k] = (phi[i][j][k] - phi[i][j][k + 1]) * per_dx;
            }
        }
    }

    this->updateReferenceEfield();

    for(auto& child : children) {
        child->updateEfield();
    }
}

//! Reference Efield (ノード上で定義される電場) を更新する
void RootGrid::updateReferenceEfield() {
    auto& ex = field->getEx();
    auto& ey = field->getEy();
    auto& ez = field->getEz();
    auto& phi = field->getPhi();
    auto& exref = field->getExRef();
    auto& eyref = field->getEyRef();
    auto& ezref = field->getEzRef();
    const size_t cx_with_glue = ex.shape()[0] + 1; // nx + 2
    const size_t cy_with_glue = ey.shape()[1] + 1;
    const size_t cz_with_glue = ez.shape()[2] + 1;

    //! reference 更新
    #pragma omp parallel shared(ex, ey, ez, exref, eyref, ezref)
    {
        #pragma omp for
        for(size_t i = 1; i < cx_with_glue - 1; ++i){
            for(size_t j = 1; j < cy_with_glue - 1; ++j){
                for(size_t k = 1; k < cz_with_glue - 1; ++k){
                    exref[i][j][k] = 0.5 * (ex[i-1][j][k] + ex[i][j][k]);
                    eyref[i][j][k] = 0.5 * (ey[i][j-1][k] + ey[i][j][k]);
                    ezref[i][j][k] = 0.5 * (ez[i][j][k-1] + ez[i][j][k]);
                }
            }
        }

        //! 外側境界の条件設定
        if (Environment::isBoundary(AXIS::x, AXIS_SIDE::low)) {
            #pragma omp for
            for(size_t j = 1; j < cy_with_glue - 1; ++j){
                for(size_t k = 1; k < cz_with_glue - 1; ++k){
                    exref[0][j][k] = 0.0;

                    // Boundary である場合、更新 (これ正しい??)
                    exref[1][j][k] = 0.5 * (phi[3][j][k] - 4.0 * phi[2][j][k] + 3.0 * phi[1][j][k]);
                }
            }
        }

        if (Environment::isBoundary(AXIS::x, AXIS_SIDE::up)) {
            #pragma omp for
            for(size_t j = 1; j < cy_with_glue - 1; ++j){
                for(size_t k = 1; k < cz_with_glue - 1; ++k){
                    exref[cx_with_glue - 1][j][k] = 0.0;
                    exref[cx_with_glue - 2][j][k] = 0.5 * (-phi[cx_with_glue - 4][j][k] + 4.0 * phi[cx_with_glue - 3][j][k] - 3.0 * phi[cx_with_glue - 2][j][k]);
                }
            }
        }

        if (Environment::isBoundary(AXIS::y, AXIS_SIDE::low)) {
            #pragma omp for
            for(size_t i = 1; i < cx_with_glue - 1; ++i){
                for(size_t k = 1; k < cz_with_glue - 1; ++k){
                    eyref[i][0][k] = 0.0;
                    eyref[i][1][k] = 0.5 * (phi[i][3][k] - 4.0 * phi[i][2][k] + 3.0 * phi[i][1][k]);
                }
            }
        }

        if (Environment::isBoundary(AXIS::y, AXIS_SIDE::up)) {
            #pragma omp for
            for(size_t i = 1; i < cx_with_glue - 1; ++i){
                for(size_t k = 1; k < cz_with_glue - 1; ++k){
                    eyref[i][cy_with_glue - 1][k] = 0.0;
                    eyref[i][cy_with_glue - 2][k] = 0.5 * (-phi[i][cy_with_glue - 4][k] + 4.0 * phi[i][cy_with_glue - 3][k] - 3.0 * phi[i][cy_with_glue - 2][k]);
                }
            }
        }

        if (Environment::isBoundary(AXIS::z, AXIS_SIDE::low)) {
            #pragma omp for
            for(size_t i = 1; i < cx_with_glue - 1; ++i){
                for(size_t j = 1; j < cy_with_glue - 1; ++j){
                    ezref[i][j][0] = 0.0;
                    ezref[i][j][1] = 0.5 * (phi[i][j][3] - 4.0 * phi[i][j][2] + 3.0 * phi[i][j][1]);
                }
            }
        }

        if (Environment::isBoundary(AXIS::z, AXIS_SIDE::up)) {
            #pragma omp for
            for(size_t i = 1; i < cx_with_glue - 1; ++i){
                for(size_t j = 1; j < cy_with_glue - 1; ++j){
                    ezref[i][j][cz_with_glue - 1] = 0.0;
                    ezref[i][j][cz_with_glue - 2] = 0.5 * (-phi[i][j][cz_with_glue - 4] + 4.0 * phi[i][j][cz_with_glue - 3] - 3.0 * phi[i][j][cz_with_glue - 2]);
                }
            }
        }
    }

    MPIw::Environment::sendRecvField(exref);
    MPIw::Environment::sendRecvField(eyref);
    MPIw::Environment::sendRecvField(ezref);
}

void RootGrid::updateEfieldFDTD() {
    this->updateEfieldFDTDMur1();
}

void RootGrid::updateEfieldFDTDMur1() {
    auto& ex = field->getEx();
    auto& ey = field->getEy();
    auto& ez = field->getEz();
    auto& bx = field->getBx();
    auto& by = field->getBy();
    auto& bz = field->getBz();
    auto& jx = field->getJx();
    auto& jy = field->getJy();
    auto& jz = field->getJz();
    const size_t cx_with_glue = ex.shape()[0] + 1; // nx + 2
    const size_t cy_with_glue = ey.shape()[1] + 1;
    const size_t cz_with_glue = ez.shape()[2] + 1;
    const double dt_per_eps0 = dt / Normalizer::eps0;
    const double dt_per_mu0_eps0_dx = dt_per_eps0 / (Normalizer::mu0 * dx);

    static const bool is_boundary[6] = {
        Environment::isBoundary(AXIS::x, AXIS_SIDE::low),
        Environment::isBoundary(AXIS::x, AXIS_SIDE::up),
        Environment::isBoundary(AXIS::y, AXIS_SIDE::low),
        Environment::isBoundary(AXIS::y, AXIS_SIDE::up),
        Environment::isBoundary(AXIS::z, AXIS_SIDE::low),
        Environment::isBoundary(AXIS::z, AXIS_SIDE::up),
    };

    const double mur_coeff = (Normalizer::c * dt - dx) / (Normalizer::c * dt + dx);

    //! とりあえず Mur 1次
    for(size_t i = 1; i < cx_with_glue - 1; ++i){
        for(size_t j = 1; j < cy_with_glue - 1; ++j){
            for(size_t k = 1; k < cz_with_glue - 1; ++k){
                //! 各方向には1つ少ないのでcx-1まで
                if(i < cx_with_glue - 2) {
                    if (i == 1 && is_boundary[0]) {
                        const double new_value = ex[i][j][k] - jx[i][j][k] * dt_per_eps0 +
                            dt_per_mu0_eps0_dx * (bz[i][j][k] - bz[i][j - 1][k] - by[i][j][k] + by[i][j][k - 1]);

                        ex[i - 1][j][k] = ex[i][j][k] + mur_coeff * (new_value - ex[i - 1][j][k]);
                        ex[i][j][k] = new_value;
                    } else if (i == cx_with_glue - 2 && is_boundary[1]) {
                        const double new_value = ex[i][j][k] - jx[i][j][k] * dt_per_eps0 +
                            dt_per_mu0_eps0_dx * (bz[i][j][k] - bz[i][j - 1][k] - by[i][j][k] + by[i][j][k - 1]);

                        ex[i + 1][j][k] = ex[i][j][k] + mur_coeff * (new_value - ex[i + 1][j][k]);
                        ex[i][j][k] = new_value;
                    } else {
                        ex[i][j][k] = ex[i][j][k] - jx[i][j][k] * dt_per_eps0 +
                            dt_per_mu0_eps0_dx * (bz[i][j][k] - bz[i][j - 1][k] - by[i][j][k] + by[i][j][k - 1]);
                    }
                }
                if(j < cy_with_glue - 2) {
                    if (j == 1 && is_boundary[2]) {
                        const double new_value = ey[i][j][k] - jy[i][j][k] * dt_per_eps0 +
                            dt_per_mu0_eps0_dx * (bx[i][j][k] - bx[i][j][k - 1] - bz[i][j][k] + bz[i - 1][j][k]);

                        ey[i][j - 1][k] = ey[i][j][k] + mur_coeff * (new_value - ey[i][j - 1][k]);
                        ey[i][j][k] = new_value;
                    } else if (j == cy_with_glue - 2 && is_boundary[3]) {
                        const double new_value = ey[i][j][k] - jy[i][j][k] * dt_per_eps0 +
                            dt_per_mu0_eps0_dx * (bx[i][j][k] - bx[i][j][k - 1] - bz[i][j][k] + bz[i - 1][j][k]);

                        ey[i][j + 1][k] = ey[i][j][k] + mur_coeff * (new_value - ey[i][j + 1][k]);
                        ey[i][j][k] = new_value;
                    } else {
                        ey[i][j][k] = ey[i][j][k] - jy[i][j][k] * dt_per_eps0 +
                            dt_per_mu0_eps0_dx * (bx[i][j][k] - bx[i][j][k - 1] - bz[i][j][k] + bz[i - 1][j][k]);
                    }
                }
                if(k < cz_with_glue - 2) {
                    if (k == 1 && is_boundary[4]) {
                        const double new_value = ez[i][j][k] - jz[i][j][k] * dt_per_eps0 +
                            dt_per_mu0_eps0_dx * (by[i][j][k] - by[i - 1][j][k] - bx[i][j][k] + bx[i][j - 1][k]);

                        ez[i][j][k - 1] = ez[i][j][k] + mur_coeff * (new_value - ez[i][j][k - 1]);
                        ez[i][j][k] = new_value;
                    } else if (k == cz_with_glue - 2 && is_boundary[5]) {
                        const double new_value = ez[i][j][k] - jz[i][j][k] * dt_per_eps0 +
                            dt_per_mu0_eps0_dx * (by[i][j][k] - by[i - 1][j][k] - bx[i][j][k] + bx[i][j - 1][k]);

                        ez[i][j][k + 1] = ez[i][j][k] + mur_coeff * (new_value - ez[i][j][k + 1]);
                        ez[i][j][k] = new_value;
                    } else {
                        ez[i][j][k] = ez[i][j][k] - jz[i][j][k] * dt_per_eps0 +
                            dt_per_mu0_eps0_dx * (by[i][j][k] - by[i - 1][j][k] - bx[i][j][k] + bx[i][j - 1][k]);
                    }
                }
            }
        }
    }

    // FDTDの場合は通信が必要になる
    MPIw::Environment::sendRecvField(ex);
    MPIw::Environment::sendRecvField(ey);
    MPIw::Environment::sendRecvField(ez);

    //! Reference 更新
    this->updateReferenceEfield();
}

void RootGrid::updateBfield() {
    this->updateBfieldMur1();
}

void RootGrid::updateBfieldMur1() {
    const double dt_per_dx = dt / dx;

    // const double epsilon_r = 1.0; //! 比誘電率
    // const double sigma = 1.0; //! 導電率 (各Faceでの)
    // const double mu_r = 1.0; //! 透磁率 (各Faceでの)
    // const double sigma_m = 0.0; //! 導磁率?

    // 磁束密度更新時の係数
    // const double d1 = mu_r/(mu_r + sigma_m * dt / mu0);
    // const double d2 = dt/(mu_r + sigma_m * dt / mu0);

    auto& bx = field->getBx();
    auto& by = field->getBy();
    auto& bz = field->getBz();
    auto& ex = field->getEx();
    auto& ey = field->getEy();
    auto& ez = field->getEz();

    const size_t cx_with_glue = bx.shape()[0];
    const size_t cy_with_glue = by.shape()[1];
    const size_t cz_with_glue = bz.shape()[2];

    static const bool is_boundary[6] = {
        Environment::isBoundary(AXIS::x, AXIS_SIDE::low),
        Environment::isBoundary(AXIS::x, AXIS_SIDE::up),
        Environment::isBoundary(AXIS::y, AXIS_SIDE::low),
        Environment::isBoundary(AXIS::y, AXIS_SIDE::up),
        Environment::isBoundary(AXIS::z, AXIS_SIDE::low),
        Environment::isBoundary(AXIS::z, AXIS_SIDE::up),
    };

    const double mur_coeff = (Normalizer::c * dt - dx) / (Normalizer::c * dt + dx);

    //! 0とcy + 1, 0とcz + 1はglueなので更新しなくてよい
    for(int i = 1; i < cx_with_glue - 1; ++i){
        for(int j = 1; j < cy_with_glue - 1; ++j){
            for(int k = 1; k < cz_with_glue - 1; ++k){
                if (j != cy_with_glue - 2 && k != cz_with_glue - 2) {
                    if (i == 1 && is_boundary[0]) {
                        const double new_value = bx[i][j][k] + dt_per_dx * ((ey[i][j][k + 1] - ey[i][j][k]) - (ez[i][j + 1][k] - ez[i][j][k]));

                        bx[i - 1][j][k] = bx[i][j][k] + mur_coeff * (new_value - bx[i - 1][j][k]);
                        bx[i][j][k] = new_value;
                    } else if (i == cx_with_glue - 2 && is_boundary[1]) {
                        const double new_value = bx[i][j][k] + dt_per_dx * ((ey[i][j][k + 1] - ey[i][j][k]) - (ez[i][j + 1][k] - ez[i][j][k]));

                        bx[i + 1][j][k] = bx[i][j][k] + mur_coeff * (new_value - bx[i + 1][j][k]);
                        bx[i][j][k] = new_value;
                    } else {
                        bx[i][j][k] = bx[i][j][k] + dt_per_dx * ((ey[i][j][k + 1] - ey[i][j][k]) - (ez[i][j + 1][k] - ez[i][j][k]));
                    }
                }

                if (i != cx_with_glue - 2 && k != cz_with_glue - 2) {
                    if (j == 1 && is_boundary[2]) {
                        const double new_value = by[i][j][k] + dt_per_dx * ((ez[i + 1][j][k] - ez[i][j][k]) - (ex[i][j][k + 1] - ex[i][j][k]));

                        by[i][j - 1][k] = by[i][j][k] + mur_coeff * (new_value - by[i][j - 1][k]);
                        by[i][j][k] = new_value;
                    } else if (j == cy_with_glue - 2 && is_boundary[3]) {
                        const double new_value = by[i][j][k] + dt_per_dx * ((ez[i + 1][j][k] - ez[i][j][k]) - (ex[i][j][k + 1] - ex[i][j][k]));

                        by[i][j + 1][k] = by[i][j][k] + mur_coeff * (new_value - by[i][j + 1][k]);
                        by[i][j][k] = new_value;
                    } else {
                        by[i][j][k] = by[i][j][k] + dt_per_dx * ((ez[i + 1][j][k] - ez[i][j][k]) - (ex[i][j][k + 1] - ex[i][j][k]));
                    }
                }

                if (i != cx_with_glue - 2 && j != cy_with_glue - 2) {
                    if (k == 1 && is_boundary[4]) {
                        const double new_value = bz[i][j][k] + dt_per_dx * ((ex[i][j + 1][k] - ex[i][j][k]) - (ey[i + 1][j][k] - ey[i][j][k]));

                        bz[i][j][k - 1] = bz[i][j][k] + mur_coeff * (new_value - bz[i][j][k - 1]);
                        bz[i][j][k] = new_value;
                    } else if (k == cz_with_glue - 2 && is_boundary[5]) {
                        const double new_value = bz[i][j][k] + dt_per_dx * ((ex[i][j + 1][k] - ex[i][j][k]) - (ey[i + 1][j][k] - ey[i][j][k]));

                        bz[i][j][k + 1] = bz[i][j][k] + mur_coeff * (new_value - bz[i][j][k + 1]);
                        bz[i][j][k] = new_value;
                    } else {
                        bz[i][j][k] = bz[i][j][k] + dt_per_dx * ((ex[i][j + 1][k] - ex[i][j][k]) - (ey[i + 1][j][k] - ey[i][j][k]));
                    }
                }
            }
        }
    }

    MPIw::Environment::sendRecvField(bx);
    MPIw::Environment::sendRecvField(by);
    MPIw::Environment::sendRecvField(bz);

    this->updateReferenceBfield();
}

void RootGrid::updateReferenceBfield() {
    auto& bx = field->getBx();
    auto& by = field->getBy();
    auto& bz = field->getBz();
    auto& bxref = field->getBxRef();
    auto& byref = field->getByRef();
    auto& bzref = field->getBzRef();
    const size_t cx_with_glue = bxref.shape()[0];
    const size_t cy_with_glue = bxref.shape()[1];
    const size_t cz_with_glue = bxref.shape()[2];

    static const bool is_boundary[6] = {
        Environment::isBoundary(AXIS::x, AXIS_SIDE::low),
        Environment::isBoundary(AXIS::x, AXIS_SIDE::up),
        Environment::isBoundary(AXIS::y, AXIS_SIDE::low),
        Environment::isBoundary(AXIS::y, AXIS_SIDE::up),
        Environment::isBoundary(AXIS::z, AXIS_SIDE::low),
        Environment::isBoundary(AXIS::z, AXIS_SIDE::up),
    };

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

void RootGrid::updateDensity(void) {
    auto& density = field->getDensity();
    Utils::initializeRhoArray(density);

    for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
        const auto& parray = particles[pid];
        const auto size = Environment::getParticleType(pid)->getSize();

        for(int pnum = 0; pnum < parray.size(); ++pnum){
            const Particle& p = parray[pnum];

            if (p.isValid) {
                Position pos(p);
                density[pid][pos.i][pos.j][pos.k] += size;
            }
        }
    }

    for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
        MPIw::Environment::sendRecvCellScalar(density[pid]);
    }
}

//! 粒子の位置から電荷を空間電荷にする
void RootGrid::updateRho() {
    RhoArray& rho = field->getRho();

#ifdef CHARGE_CONSERVATION
    // 電荷保存則をcheckするため、古いrhoを保持する
    auto old_rho = rho;
#endif

    auto time_counter = Utils::TimeCounter::getInstance();

    //! rhoを初期化
    time_counter->switchTo("updateRho/initializeRhoArray");
    Utils::initializeRhoArray(rho);

    time_counter->switchTo("updateRho/mainLoop");
    #pragma omp parallel for
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

    time_counter->switchTo("updateRho/sumWholeCharge");
    for(auto& obj : objects) {
        if (obj.isDefined()) {
            obj.sumWholeCharge();
        }
    }

    time_counter->switchTo("updateRho/totalCharge Computation");
    for(int pid = 1; pid < Environment::num_of_particle_types + 1; ++pid) {
        for(int i = 0; i < nx + 2; ++i) {
            for(int j = 0; j < ny + 2; ++j) {
                for(int k = 0; k < nz + 2; ++k) {
                    rho[0][i][j][k] += rho[pid][i][j][k];
                }
            }
        }
    }

    //! 子の電荷の対応部分をRootにコピー(いらないかも)
    time_counter->switchTo("updateRho/childRho");
    for(auto& child : children) {
        child->updateRho();
        child->copyRhoToParent();
    }

    //! rho を隣に送る
    time_counter->switchTo("updateRho/sendRhoAll");
    for(int pid = 0; pid < Environment::num_of_particle_types + 1; ++pid) {
        MPIw::Environment::sendRecvNodeScalar(rho[pid]);
    }

    //! 物体が存在する場合、電荷再配分が必要になる
    if (objects.size() > 0) {
        //! 一度 Poisson を解いて phi を更新
        time_counter->switchTo("solvePoisson");
        solvePoisson();

        time_counter->switchTo("updateRho/redistributeCharge");
        for(auto& obj : objects) {
            if (obj.isDefined()) obj.redistributeCharge(rho, field->getPhi());
        }

        time_counter->switchTo("updateRho/sendRhoZero");
        //! 既に隣のFieldとRhoが同期しているため、再配分後のrhoの送信はField形式でいい
        MPIw::Environment::sendRecvField(rho[0]);
    }

#ifdef CHARGE_CONSERVATION
    if (Environment::solver_type == "EM") {
        field->checkChargeConservation(old_rho, 1.0, dx);
    }
#endif

    time_counter->end();
}