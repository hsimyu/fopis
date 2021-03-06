#include "grid.hpp"
#include "normalizer.hpp"
#include "field.hpp"
#include "utils.hpp"
#include "dataio.hpp"

void RootGrid::initializeField() {
    Grid::initializeField();
    
    //! RootGridかつEMの場合、Dampingのための領域を
    //! EとBに確保する必要がある
    if (Environment::isEMMode()) {
        const int cx = nx + 2;
        const int cy = ny + 2;
        const int cz = nz + 2;

        auto& damping_length = Environment::getDampingLength();
        damping_length.L_dx = nx;
        damping_length.L_dy = ny;
        damping_length.L_dz = nz;

        int x_low = (Environment::isBoundary(AXIS::x, AXIS_SIDE::low)) ? -(damping_length.L_dx) : 0;
        int x_up = (Environment::isBoundary(AXIS::x, AXIS_SIDE::up)) ? cx + damping_length.L_dx : cx;
        int y_low = (Environment::isBoundary(AXIS::y, AXIS_SIDE::low)) ? -(damping_length.L_dy) : 0;
        int y_up = (Environment::isBoundary(AXIS::y, AXIS_SIDE::up)) ? cy + damping_length.L_dy : cy;
        int z_low = (Environment::isBoundary(AXIS::z, AXIS_SIDE::low)) ? -(damping_length.L_dz) : 0;
        int z_up = (Environment::isBoundary(AXIS::z, AXIS_SIDE::up)) ? cz + damping_length.L_dz : cz;

        tdArray::extent_gen extents;
        using range = tdArray::extent_range;
        field->getEx().resize(extents[range(x_low, x_up - 1)][range(y_low, y_up    )][range(z_low, z_up    )]);
        field->getEy().resize(extents[range(x_low, x_up    )][range(y_low, y_up - 1)][range(z_low, z_up    )]);
        field->getEz().resize(extents[range(x_low, x_up    )][range(y_low, y_up    )][range(z_low, z_up - 1)]);

        field->getBx().resize(extents[range(x_low, x_up    )][range(y_low, y_up - 1)][range(z_low, z_up - 1)]);
        field->getBy().resize(extents[range(x_low, x_up - 1)][range(y_low, y_up    )][range(z_low, z_up - 1)]);
        field->getBz().resize(extents[range(x_low, x_up - 1)][range(y_low, y_up - 1)][range(z_low, z_up    )]);
    }
}

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
    constexpr double required_error = 1.0e-7;

    auto time_counter = Utils::TimeCounter::getInstance();
    auto& phi = field->getPhi();
    auto& poisson_error = field->getPoissonError();

    time_counter->begin("solvePoisson/updatePoissonErrorPre");
    poisson_error = phi;

    size_t loop = 1;
    for(; loop <= loopnum; ++loop) {
        time_counter->switchTo("solvePoisson/mainLoop");
        this->doPSOR();

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
    for(int i = 1; i < phi.shape()[0] - 1; ++i){
        for(int j = 1; j < phi.shape()[1] - 1; ++j){
            for(int k = 1; k < phi.shape()[2] - 1; ++k){
                poisson_error[i][j][k] = phi[i][j][k] - poisson_error[i][j][k];
            }
        }
    }
    time_counter->end();
}

void RootGrid::doPSOR() {
    auto time_counter = Utils::TimeCounter::getInstance();

    auto& phi = field->getPhi();
    const int minimum_radius = std::min({Environment::nx, Environment::ny, Environment::nz});
    const double omega = 2.0/(1.0 + sin(M_PI/minimum_radius)); // spectral radius
    const double rho_coeff = pow(dx, 2) / Normalizer::eps0;

    size_t i_begin = 1; size_t i_end = 1;
    size_t j_begin = 1; size_t j_end = 1;
    size_t k_begin = 1; size_t k_end = 1;

    //! 最初の部分を計算して
    time_counter->switchTo("solvePoisson/mainLoop/doSORPartial");
    this->doSORPartial(omega, rho_coeff, i_begin, i_end, j_begin, j_end, k_begin, k_end);

    //! prev-xに送る
    time_counter->switchTo("solvePoisson/mainLoop/sendPhi");
    MPIw::Environment::sendRecvPartialPhi(phi, i_begin, i_end, j_begin, j_end, k_begin, k_end);

    const size_t cx_with_glue = phi.shape()[0];
    const size_t cy_with_glue = phi.shape()[1];
    const size_t cz_with_glue = phi.shape()[2];

    //! x-edgeを更新して
    i_begin = 2; i_end = cx_with_glue - 2;
    time_counter->switchTo("solvePoisson/mainLoop/doSORPartial");
    this->doSORPartial(omega, rho_coeff, i_begin, i_end, j_begin, j_end, k_begin, k_end);

    //! prev-yに送る
    time_counter->switchTo("solvePoisson/mainLoop/sendPhi");
    //! Glueノード分も送るので, x座標の始点を-1, 終点を+1はする
    MPIw::Environment::sendRecvPartialPhi(phi, i_begin - 1, i_end + 1, j_begin, j_end, k_begin, k_end);

    //! y-edgeを更新して
    i_begin = 1; i_end = 1;
    j_begin = 2; j_end = cy_with_glue - 2;
    time_counter->switchTo("solvePoisson/mainLoop/doSORPartial");
    this->doSORPartial(omega, rho_coeff, i_begin, i_end, j_begin, j_end, k_begin, k_end);

    //! prev-xに送る
    time_counter->switchTo("solvePoisson/mainLoop/sendPhi");
    MPIw::Environment::sendRecvPartialPhi(phi, i_begin, i_end, j_begin, j_end + 1, k_begin, k_end);

    //! xy面を更新して
    i_begin = 2; i_end = cx_with_glue - 2;
    time_counter->switchTo("solvePoisson/mainLoop/doSORPartial");
    this->doSORPartial(omega, rho_coeff, i_begin, i_end, j_begin, j_end, k_begin, k_end);

    //! prev-zに送る
    time_counter->switchTo("solvePoisson/mainLoop/sendPhi");
    MPIw::Environment::sendRecvPartialPhi(phi, i_begin - 1, i_end + 1, j_begin - 1, j_end + 1, k_begin, k_end);

    //! z-edgeを更新して
    i_begin = 1; i_end = 1;
    j_begin = 1; j_end = 1;
    k_begin = 2; k_end = cz_with_glue - 2;
    time_counter->switchTo("solvePoisson/mainLoop/doSORPartial");
    this->doSORPartial(omega, rho_coeff, i_begin, i_end, j_begin, j_end, k_begin, k_end);

    //! prev-xに送る
    time_counter->switchTo("solvePoisson/mainLoop/sendPhi");
    MPIw::Environment::sendRecvPartialPhi(phi, i_begin, i_end, j_begin, j_end, k_begin, k_end + 1);

    //! xz面を更新して
    i_begin = 2; i_end = cx_with_glue - 2;
    time_counter->switchTo("solvePoisson/mainLoop/doSORPartial");
    this->doSORPartial(omega, rho_coeff, i_begin, i_end, j_begin, j_end, k_begin, k_end);

    //! prev-yに送る
    time_counter->switchTo("solvePoisson/mainLoop/sendPhi");
    MPIw::Environment::sendRecvPartialPhi(phi, i_begin - 1, i_end + 1, j_begin, j_end, k_begin, k_end + 1);

    //! yz面を更新して
    i_begin = 1; i_end = 1;
    j_begin = 2; j_end = cy_with_glue - 2;
    time_counter->switchTo("solvePoisson/mainLoop/doSORPartial");
    this->doSORPartial(omega, rho_coeff, i_begin, i_end, j_begin, j_end, k_begin, k_end);

    //! prev-xに送る
    time_counter->switchTo("solvePoisson/mainLoop/sendPhi");
    MPIw::Environment::sendRecvPartialPhi(phi, i_begin, i_end, j_begin, j_end + 1, k_begin, k_end + 1);

    //! 内部の値を更新
    i_begin = 2; i_end = cx_with_glue - 2;
    time_counter->switchTo("solvePoisson/mainLoop/doSORPartial");
    #pragma omp parallel
    {
        this->doSORPartial(omega, rho_coeff, i_begin, i_end, j_begin, j_end, k_begin, k_end);
    }

    // xy面をnext-zへ
    i_begin = 0; i_end = cx_with_glue - 1;
    j_begin = 0; j_end = cy_with_glue - 1;
    k_begin = cz_with_glue - 2; k_end = cz_with_glue - 2;
    time_counter->switchTo("solvePoisson/mainLoop/sendPhi");
    MPIw::Environment::sendRecvPartialPhi(phi, i_begin, i_end, j_begin, j_end, k_begin, k_end);

    // xz面をnext-yへ
    j_begin = cy_with_glue - 2; j_end = cy_with_glue - 2;
    k_begin = 0; k_end = cz_with_glue - 1;
    time_counter->switchTo("solvePoisson/mainLoop/sendPhi");
    MPIw::Environment::sendRecvPartialPhi(phi, i_begin, i_end, j_begin, j_end, k_begin, k_end);

    // yz面をnext-xへ
    i_begin = cx_with_glue - 2; i_end = cx_with_glue - 2;
    j_begin = 0; j_end = cy_with_glue - 1;
    time_counter->switchTo("solvePoisson/mainLoop/sendPhi");
    MPIw::Environment::sendRecvPartialPhi(phi, i_begin, i_end, j_begin, j_end, k_begin, k_end);
}

void RootGrid::doSORPartial(const double omega, const double rho_coeff, const size_t i_begin, const size_t i_end, const size_t j_begin, const size_t j_end, const size_t k_begin, const size_t k_end) {
    auto& phi = field->getPhi();
    auto& rho = field->getRho();

    size_t i_begin_temp = 0; size_t i_end_temp = 0;
    size_t j_begin_temp = 0; size_t j_end_temp = 0;
    size_t k_begin_temp = 0; size_t k_end_temp = 0;

    //! 境界は除く
    if (Environment::isBoundary(AXIS::x, AXIS_SIDE::low) && i_begin == 1) i_begin_temp = 1;
    if (Environment::isBoundary(AXIS::x, AXIS_SIDE::up) && i_end == phi.shape()[0] - 2) i_end_temp = -1;

    if (Environment::isBoundary(AXIS::y, AXIS_SIDE::low) && j_begin == 1) j_begin_temp = 1;
    if (Environment::isBoundary(AXIS::y, AXIS_SIDE::up) && j_end == phi.shape()[1] - 2) j_end_temp = -1;

    if (Environment::isBoundary(AXIS::z, AXIS_SIDE::low) && k_begin == 1) k_begin_temp = 1;
    if (Environment::isBoundary(AXIS::z, AXIS_SIDE::up) && k_end == phi.shape()[2] - 2) k_end_temp = -1;

    size_t i_size = (i_end + i_end_temp) - (i_begin + i_begin_temp) + 1;

    if (i_size > 2) {
        // red-black
        #pragma omp for
        for(size_t i = i_begin + i_begin_temp; i <= i_end + i_end_temp; i += 2){
            for(size_t j = j_begin + j_begin_temp; j <= j_end + j_end_temp; ++j){
                for(size_t k = k_begin + k_begin_temp; k <= k_end + k_end_temp; ++k){
                    phi[i][j][k] = (1.0 - omega) * phi[i][j][k] + omega * (phi[i+1][j][k] + phi[i-1][j][k] + phi[i][j+1][k] + phi[i][j-1][k] + phi[i][j][k+1] + phi[i][j][k-1] + rho_coeff * rho[0][i][j][k])/6.0;
                }
            }
        }

        #pragma omp for
        for(size_t i = i_begin + i_begin_temp + 1; i <= i_end + i_end_temp; i += 2){
            for(size_t j = j_begin + j_begin_temp; j <= j_end + j_end_temp; ++j){
                for(size_t k = k_begin + k_begin_temp; k <= k_end + k_end_temp; ++k){
                    phi[i][j][k] = (1.0 - omega) * phi[i][j][k] + omega * (phi[i+1][j][k] + phi[i-1][j][k] + phi[i][j+1][k] + phi[i][j-1][k] + phi[i][j][k+1] + phi[i][j][k-1] + rho_coeff * rho[0][i][j][k])/6.0;
                }
            }
        }
    } else {
        for(size_t i = i_begin + i_begin_temp; i <= i_end + i_end_temp; ++i){
            for(size_t j = j_begin + j_begin_temp; j <= j_end + j_end_temp; ++j){
                for(size_t k = k_begin + k_begin_temp; k <= k_end + k_end_temp; ++k){
                    phi[i][j][k] = (1.0 - omega) * phi[i][j][k] + omega * (phi[i+1][j][k] + phi[i-1][j][k] + phi[i][j+1][k] + phi[i][j-1][k] + phi[i][j][k+1] + phi[i][j][k-1] + rho_coeff * rho[0][i][j][k])/6.0;
                }
            }
        }
    }
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

    const size_t cx_with_glue = phi.shape()[0]; // nx + 2
    const size_t cy_with_glue = phi.shape()[1];
    const size_t cz_with_glue = phi.shape()[2];
    const double per_dx = 1.0 / dx;

    //! @note:隣と通信しなくてもいい
    //! phiが通信してあるため、端の要素を通信なしで計算可能
    #pragma omp parallel
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
    const size_t cx_with_glue = exref.shape()[0]; // nx + 2
    const size_t cy_with_glue = eyref.shape()[1];
    const size_t cz_with_glue = ezref.shape()[2];

    //! reference 更新
    #pragma omp parallel
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
    // this->updateEfieldFDTDMur1();
    this->updateEfieldFDTDDamping();
}

void RootGrid::updateBfield() {
    // this->updateBfieldMur1();
    this->updateBfieldDamping();
}

void RootGrid::updateEfieldFDTDDamping() {
    auto& ex = field->getEx();
    auto& ey = field->getEy();
    auto& ez = field->getEz();
    auto& bx = field->getBx();
    auto& by = field->getBy();
    auto& bz = field->getBz();
    auto& jx = field->getJx();
    auto& jy = field->getJy();
    auto& jz = field->getJz();
    const double dt_per_eps0 = dt / Normalizer::eps0;
    const double dt_per_mu0_eps0_dx = dt_per_eps0 / (Normalizer::mu0 * dx);

    auto bases = ex.index_bases();
    auto shapes = ex.shape();

    //! 実領域上端側のGlueのindex + 1
    const int real_x_index = jx.shape()[0] + 1;
    const int real_y_index = jy.shape()[1] + 1;
    const int real_z_index = jz.shape()[2] + 1;

    //! Damping 領域も含んだ領域上限のindex + 1
    const int max_x_index = bases[0] + shapes[0] + 1; // ex.shape()を使っているのでx方向だけ+1する
    const int max_y_index = bases[1] + shapes[1];
    const int max_z_index = bases[2] + shapes[2];

    //! 1以上、rx以下なら実領域（Glueエッジも自動的に除かれる）
    auto is_real_region = [rx = real_x_index, ry = real_y_index, rz = real_z_index](const int i, const int j, const int k)-> bool {
        return (
            ((i > 0) && (i < rx - 1)) &&
            ((j > 0) && (j < ry - 1)) &&
            ((k > 0) && (k < rz - 1))
        );
    };

    constexpr double masking_parameter = 1.0;
    auto one_masking_factor = [](const double r, const double damping_length, const int real_index) -> double {
        double dist_r = (r <= 0) ? - r + 1: 
            (r >= real_index - 1) ? r - static_cast<double>(real_index) + 2.0 : 0.0;
        return 1.0 - std::pow(masking_parameter * dist_r / damping_length, 2);
    };
    auto x_masking_factor = [rx = real_x_index, dl = static_cast<double>(Environment::getDampingLength().L_dx), om = &one_masking_factor](const double x) -> double { return (*om)(x, dl, rx); };
    auto y_masking_factor = [ry = real_y_index, dl = static_cast<double>(Environment::getDampingLength().L_dy), om = &one_masking_factor](const double y) -> double { return (*om)(y, dl, ry); };
    auto z_masking_factor = [rz = real_z_index, dl = static_cast<double>(Environment::getDampingLength().L_dz), om = &one_masking_factor](const double z) -> double { return (*om)(z, dl, rz); };

    auto masking_factor = [xm = &x_masking_factor, ym = &y_masking_factor, zm = &z_masking_factor](const int x, const int y, const int z) -> double {
        return (*xm)(x) * (*ym)(y) * (*zm)(z);
    };

    for(int i = bases[0] + 1; i < max_x_index - 1; ++i){
        for(int j = bases[1] + 1; j < max_y_index - 1; ++j){
            for(int k = bases[2] + 1; k < max_z_index - 1; ++k){
                if (is_real_region(i, j, k)) {
                    //! Edge要素のため、各方向には1つ少ない = rx - 2 までしか要素がない = rx - 3 までの更新でいい
                    if(i < real_x_index - 2) {
                        ex[i][j][k] = ex[i][j][k] - jx[i][j][k] * dt_per_eps0 +
                            dt_per_mu0_eps0_dx * (bz[i][j][k] - bz[i][j - 1][k] - by[i][j][k] + by[i][j][k - 1]);
                    } else {
                        const double mask = 0.5 * (masking_factor(i, j, k) + masking_factor(i + 1, j, k));
                        ex[i][j][k] = mask * (
                            ex[i][j][k] + dt_per_mu0_eps0_dx * (bz[i][j][k] - bz[i][j - 1][k] - by[i][j][k] + by[i][j][k - 1])
                        );
                    }

                    if(j < real_y_index - 2) {
                        ey[i][j][k] = ey[i][j][k] - jy[i][j][k] * dt_per_eps0 +
                            dt_per_mu0_eps0_dx * (bx[i][j][k] - bx[i][j][k - 1] - bz[i][j][k] + bz[i - 1][j][k]);
                    } else {
                        const double mask = 0.5 * (masking_factor(i, j, k) + masking_factor(i, j + 1, k));
                        ey[i][j][k] = mask * (
                            ey[i][j][k] + dt_per_mu0_eps0_dx * (bx[i][j][k] - bx[i][j][k - 1] - bz[i][j][k] + bz[i - 1][j][k])
                        );
                    }

                    if(k < real_z_index - 2) {
                        ez[i][j][k] = ez[i][j][k] - jz[i][j][k] * dt_per_eps0 +
                            dt_per_mu0_eps0_dx * (by[i][j][k] - by[i - 1][j][k] - bx[i][j][k] + bx[i][j - 1][k]);
                    } else {
                        const double mask = 0.5 * (masking_factor(i, j, k) + masking_factor(i, j, k + 1));
                        ez[i][j][k] = mask * (
                            ez[i][j][k] + dt_per_mu0_eps0_dx * (by[i][j][k] - by[i - 1][j][k] - bx[i][j][k] + bx[i][j - 1][k])
                        );
                    }
                } else {
                    if(i < max_x_index - 2) {
                        const double mask = 0.5 * (masking_factor(i, j, k) + masking_factor(i + 1, j, k));
                        ex[i][j][k] = mask * (
                            ex[i][j][k] + dt_per_mu0_eps0_dx * (bz[i][j][k] - bz[i][j - 1][k] - by[i][j][k] + by[i][j][k - 1])
                        );
                    }

                    if(j < max_y_index - 2) {
                        const double mask = 0.5 * (masking_factor(i, j, k) + masking_factor(i, j + 1, k));
                        ey[i][j][k] = mask * (
                            ey[i][j][k] + dt_per_mu0_eps0_dx * (bx[i][j][k] - bx[i][j][k - 1] - bz[i][j][k] + bz[i - 1][j][k])
                        );
                    }

                    if(k < max_z_index - 2) {
                        const double mask = 0.5 * (masking_factor(i, j, k) + masking_factor(i, j, k + 1));
                        ez[i][j][k] = mask * (
                            ez[i][j][k] + dt_per_mu0_eps0_dx * (by[i][j][k] - by[i - 1][j][k] - bx[i][j][k] + bx[i][j - 1][k])
                        );
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

void RootGrid::updateBfieldDamping() {
    const double dt_per_dx = dt / dx;

    auto& bx = field->getBx();
    auto& by = field->getBy();
    auto& bz = field->getBz();
    auto& ex = field->getEx();
    auto& ey = field->getEy();
    auto& ez = field->getEz();

    auto bases = bx.index_bases();
    auto shapes = bx.shape();

    //! 実領域上端側のGlueのindex + 1
    const int real_x_index = field->getBxRef().shape()[0];
    const int real_y_index = field->getByRef().shape()[1];
    const int real_z_index = field->getBzRef().shape()[2];

    //! Damping 領域も含んだ領域上限のindex + 1
    const int max_x_index = bases[0] + shapes[0]; // bx.shape()を使っているのでyz方向を+1する
    const int max_y_index = bases[1] + shapes[1] + 1;
    const int max_z_index = bases[2] + shapes[2] + 1;

    //! 1以上、rx以下なら実領域（Glueエッジも自動的に除かれる）
    auto is_real_region = [rx = real_x_index, ry = real_y_index, rz = real_z_index](const int i, const int j, const int k)-> bool {
        return (
            ((i > 0) && (i < rx - 1)) &&
            ((j > 0) && (j < ry - 1)) &&
            ((k > 0) && (k < rz - 1))
        );
    };

    constexpr double masking_parameter = 1.0;
    auto one_masking_factor = [](const double r, const double damping_length, const int real_index) -> double {
        double dist_r = (r <= 0) ? - r + 1: 
            (r >= real_index - 1) ? r - static_cast<double>(real_index) + 2.0 : 0.0;
        return 1.0 - std::pow(masking_parameter * dist_r / damping_length, 2);
    };
    auto x_masking_factor = [rx = real_x_index, dl = static_cast<double>(Environment::getDampingLength().L_dx), om = &one_masking_factor](const double x) -> double { return (*om)(x, dl, rx); };
    auto y_masking_factor = [ry = real_y_index, dl = static_cast<double>(Environment::getDampingLength().L_dy), om = &one_masking_factor](const double y) -> double { return (*om)(y, dl, ry); };
    auto z_masking_factor = [rz = real_z_index, dl = static_cast<double>(Environment::getDampingLength().L_dz), om = &one_masking_factor](const double z) -> double { return (*om)(z, dl, rz); };

    auto masking_factor = [xm = &x_masking_factor, ym = &y_masking_factor, zm = &z_masking_factor](const int x, const int y, const int z) -> double {
        return (*xm)(x) * (*ym)(y) * (*zm)(z);
    };

    for(int i = bases[0] + 1; i < max_x_index - 1; ++i){
        for(int j = bases[1] + 1; j < max_y_index - 1; ++j){
            for(int k = bases[2] + 1; k < max_z_index - 1; ++k){

                if (is_real_region(i, j, k)) {
                    if (j != real_y_index - 2 && k != real_z_index - 2) {
                        bx[i][j][k] = bx[i][j][k] + dt_per_dx * ((ey[i][j][k + 1] - ey[i][j][k]) - (ez[i][j + 1][k] - ez[i][j][k]));
                    } else {
                        const double mask = 0.25 * (masking_factor(i, j, k) + masking_factor(i, j + 1, k) + masking_factor(i, j, k + 1) + masking_factor(i, j + 1, k + 1));
                        bx[i][j][k] = mask * (
                            bx[i][j][k] + dt_per_dx * ((ey[i][j][k + 1] - ey[i][j][k]) - (ez[i][j + 1][k] - ez[i][j][k]))
                        );
                    }

                    if (i != real_x_index - 2 && k != real_z_index - 2) {
                        by[i][j][k] = by[i][j][k] + dt_per_dx * ((ez[i + 1][j][k] - ez[i][j][k]) - (ex[i][j][k + 1] - ex[i][j][k]));
                    } else {
                        const double mask = 0.25 * (masking_factor(i, j, k) + masking_factor(i + 1, j, k) + masking_factor(i, j, k + 1) + masking_factor(i + 1, j, k + 1));
                        by[i][j][k] = mask * (
                            by[i][j][k] + dt_per_dx * ((ez[i + 1][j][k] - ez[i][j][k]) - (ex[i][j][k + 1] - ex[i][j][k]))
                        );
                    }

                    if (i != real_x_index - 2 && j != real_y_index - 2) {
                        bz[i][j][k] = bz[i][j][k] + dt_per_dx * ((ex[i][j + 1][k] - ex[i][j][k]) - (ey[i + 1][j][k] - ey[i][j][k]));
                    } else {
                        const double mask = 0.25 * (masking_factor(i, j, k) + masking_factor(i + 1, j, k) + masking_factor(i, j + 1, k) + masking_factor(i + 1, j + 1, k));
                        bz[i][j][k] = mask * (
                            bz[i][j][k] + dt_per_dx * ((ex[i][j + 1][k] - ex[i][j][k]) - (ey[i + 1][j][k] - ey[i][j][k]))
                        );
                    }
                } else {
                    if (j != max_y_index - 2 && k != max_z_index - 2) {
                        const double mask = 0.25 * (masking_factor(i, j, k) + masking_factor(i, j + 1, k) + masking_factor(i, j, k + 1) + masking_factor(i, j + 1, k + 1));
                        bx[i][j][k] = mask * (
                            bx[i][j][k] + dt_per_dx * ((ey[i][j][k + 1] - ey[i][j][k]) - (ez[i][j + 1][k] - ez[i][j][k]))
                        );
                    }

                    if (i != max_x_index - 2 && k != max_z_index - 2) {
                        const double mask = 0.25 * (masking_factor(i, j, k) + masking_factor(i + 1, j, k) + masking_factor(i, j, k + 1) + masking_factor(i + 1, j, k + 1));
                        by[i][j][k] = mask * (
                            by[i][j][k] + dt_per_dx * ((ez[i + 1][j][k] - ez[i][j][k]) - (ex[i][j][k + 1] - ex[i][j][k]))
                        );
                    }

                    if (i != max_x_index - 2 && j != max_y_index - 2) {
                        const double mask = 0.25 * (masking_factor(i, j, k) + masking_factor(i + 1, j, k) + masking_factor(i, j + 1, k) + masking_factor(i + 1, j + 1, k));
                        bz[i][j][k] = mask * (
                            bz[i][j][k] + dt_per_dx * ((ex[i][j + 1][k] - ex[i][j][k]) - (ey[i + 1][j][k] - ey[i][j][k]))
                        );
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

    //! Mur 1次
    //! x 下部境界
    if (is_boundary[0]) {
        int i = 1; int i_neighbor = i + 1;
        for(size_t j = 1; j < cy_with_glue - 1; ++j){
            for(size_t k = 1; k < cz_with_glue - 1; ++k){
                const double new_neighbor_value = ex[i_neighbor][j][k] - jx[i_neighbor][j][k] * dt_per_eps0 +
                    dt_per_mu0_eps0_dx * (bz[i_neighbor][j][k] - bz[i_neighbor][j - 1][k] - by[i_neighbor][j][k] + by[i_neighbor][j][k - 1]);
                ex[i][j][k] = ex[i_neighbor][j][k] + mur_coeff * (new_neighbor_value - ex[i][j][k]);
            }
        }
    }

    //! x 上部境界
    if (is_boundary[1]) {
        int i = cx_with_glue - 3; int i_neighbor = i - 1;
        for(size_t j = 1; j < cy_with_glue - 1; ++j){
            for(size_t k = 1; k < cz_with_glue - 1; ++k){
                const double new_neighbor_value = ex[i_neighbor][j][k] - jx[i_neighbor][j][k] * dt_per_eps0 +
                    dt_per_mu0_eps0_dx * (bz[i_neighbor][j][k] - bz[i_neighbor][j - 1][k] - by[i_neighbor][j][k] + by[i_neighbor][j][k - 1]);
                ex[i][j][k] = ex[i_neighbor][j][k] + mur_coeff * (new_neighbor_value - ex[i][j][k]);
            }
        }
    }

    //! y下部境界
    if (is_boundary[2]) {
        int j = 1; int j_neighbor = j + 1;
        for(size_t i = 1; i < cx_with_glue - 1; ++i){
            for(size_t k = 1; k < cz_with_glue - 1; ++k){
                const double new_neighbor_value = ey[i][j_neighbor][k] - jy[i][j_neighbor][k] * dt_per_eps0 +
                    dt_per_mu0_eps0_dx * (bx[i][j_neighbor][k] - bx[i][j_neighbor][k - 1] - bz[i][j_neighbor][k] + bz[i - 1][j_neighbor][k]);

                ey[i][j][k] = ey[i][j_neighbor][k] + mur_coeff * (new_neighbor_value - ey[i][j][k]);
            }
        }
    }

    //! y上部境界
    if (is_boundary[3]) {
        int j = cy_with_glue - 3; int j_neighbor = j - 1;
        for(size_t i = 1; i < cx_with_glue - 1; ++i){
            for(size_t k = 1; k < cz_with_glue - 1; ++k){
                const double new_neighbor_value = ey[i][j_neighbor][k] - jy[i][j_neighbor][k] * dt_per_eps0 +
                    dt_per_mu0_eps0_dx * (bx[i][j_neighbor][k] - bx[i][j_neighbor][k - 1] - bz[i][j_neighbor][k] + bz[i - 1][j_neighbor][k]);

                ey[i][j][k] = ey[i][j_neighbor][k] + mur_coeff * (new_neighbor_value - ey[i][j][k]);
            }
        }
    }

    //! z下部境界
    if (is_boundary[4]) {
        int k = 1; int k_neighbor = k + 1;
        for(size_t i = 1; i < cx_with_glue - 1; ++i){
            for(size_t j = 1; j < cy_with_glue - 1; ++j){
                const double new_neighbor_value = ez[i][j][k_neighbor] - jz[i][j][k_neighbor] * dt_per_eps0 +
                    dt_per_mu0_eps0_dx * (by[i][j][k_neighbor] - by[i - 1][j][k_neighbor] - bx[i][j][k_neighbor] + bx[i][j - 1][k_neighbor]);

                ez[i][j][k] = ez[i][j][k_neighbor] + mur_coeff * (new_neighbor_value - ez[i][j][k]);
            }
        }
    }

    //! z上部境界
    if (is_boundary[5]) {
        int k = cz_with_glue - 3; int k_neighbor = k - 1;
        for(size_t i = 1; i < cx_with_glue - 1; ++i){
            for(size_t j = 1; j < cy_with_glue - 1; ++j){
                const double new_neighbor_value = ez[i][j][k_neighbor] - jz[i][j][k_neighbor] * dt_per_eps0 +
                    dt_per_mu0_eps0_dx * (by[i][j][k_neighbor] - by[i - 1][j][k_neighbor] - bx[i][j][k_neighbor] + bx[i][j - 1][k_neighbor]);

                ez[i][j][k] = ez[i][j][k_neighbor] + mur_coeff * (new_neighbor_value - ez[i][j][k]);
            }
        }
    }

    for(size_t i = 1; i < cx_with_glue - 1; ++i){
        for(size_t j = 1; j < cy_with_glue - 1; ++j){
            for(size_t k = 1; k < cz_with_glue - 1; ++k){
                //! Edge要素のため、各方向には1つ少ない = cx - 2 までしか要素がない = cx - 3 までの更新でいい(Glue Edgeを除く)
                if(i < cx_with_glue - 2) {
                    if (!( (i == 1 && is_boundary[0]) || (i == cx_with_glue - 3 && is_boundary[1]) )) {
                        ex[i][j][k] = ex[i][j][k] - jx[i][j][k] * dt_per_eps0 +
                            dt_per_mu0_eps0_dx * (bz[i][j][k] - bz[i][j - 1][k] - by[i][j][k] + by[i][j][k - 1]);
                    }
                }

                if(j < cy_with_glue - 2) {
                    if (!( (j == 1 && is_boundary[2]) && ((j == cy_with_glue - 3 && is_boundary[3])) )) {
                        ey[i][j][k] = ey[i][j][k] - jy[i][j][k] * dt_per_eps0 +
                            dt_per_mu0_eps0_dx * (bx[i][j][k] - bx[i][j][k - 1] - bz[i][j][k] + bz[i - 1][j][k]);
                    }
                }

                if(k < cz_with_glue - 2) {
                    if (!( (k == 1 && is_boundary[4]) && ((k == cz_with_glue - 3 && is_boundary[5])) )) {
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

    //! x下部境界
    /*
    if (is_boundary[0]) {
        int i = 1; int i_neighbor = i + 1;
        for(int j = 1; j < cy_with_glue - 2; ++j){
            for(int k = 1; k < cz_with_glue - 2; ++k){
                const double new_neighbor_value = bx[i_neighbor][j][k] + dt_per_dx * ((ey[i_neighbor][j][k + 1] - ey[i_neighbor][j][k]) - (ez[i_neighbor][j + 1][k] - ez[i_neighbor][j][k]));

                bx[i][j][k] = bx[i_neighbor][j][k] + mur_coeff * (new_neighbor_value - bx[i][j][k]);
            }
        }
    }

    //! x上部境界
    if (is_boundary[1]) {
        int i = cx_with_glue - 2; int i_neighbor = i - 1;
        for(int j = 1; j < cy_with_glue - 2; ++j){
            for(int k = 1; k < cz_with_glue - 2; ++k){
                const double new_neighbor_value = bx[i_neighbor][j][k] + dt_per_dx * ((ey[i_neighbor][j][k + 1] - ey[i_neighbor][j][k]) - (ez[i_neighbor][j + 1][k] - ez[i_neighbor][j][k]));

                bx[i][j][k] = bx[i_neighbor][j][k] + mur_coeff * (new_neighbor_value - bx[i][j][k]);
            }
        }
    }

    //! y下部境界
    if (is_boundary[2]) {
        int j = 1; int j_neighbor = j + 1;
        for(int i = 1; i < cx_with_glue - 1; ++i){
            for(int k = 1; k < cz_with_glue - 1; ++k){
                const double new_neighbor_value = by[i][j_neighbor][k] + dt_per_dx * ((ez[i + 1][j_neighbor][k] - ez[i][j_neighbor][k]) - (ex[i][j_neighbor][k + 1] - ex[i][j_neighbor][k]));

                by[i][j][k] = by[i][j_neighbor][k] + mur_coeff * (new_neighbor_value - by[i][j][k]);
            }
        }
    }

    //! y上部境界
    if (is_boundary[3]) {
        int j = cy_with_glue - 2; int j_neighbor = j - 1;
        for(int i = 1; i < cx_with_glue - 1; ++i){
            for(int k = 1; k < cz_with_glue - 1; ++k){
                const double new_neighbor_value = by[i][j_neighbor][k] + dt_per_dx * ((ez[i + 1][j_neighbor][k] - ez[i][j_neighbor][k]) - (ex[i][j_neighbor][k + 1] - ex[i][j_neighbor][k]));

                by[i][j][k] = by[i][j_neighbor][k] + mur_coeff * (new_neighbor_value - by[i][j][k]);
            }
        }
    }

    //! z下部境界
    if (is_boundary[4]) {
        int k = 1; int k_neighbor = k + 1;
        for(int i = 1; i < cx_with_glue - 1; ++i){
            for(int j = 1; j < cy_with_glue - 1; ++j){
                const double new_neighbor_value = bz[i][j][k_neighbor] + dt_per_dx * ((ex[i][j + 1][k_neighbor] - ex[i][j][k_neighbor]) - (ey[i + 1][j][k_neighbor] - ey[i][j][k_neighbor]));

                bz[i][j][k] = bz[i][j][k_neighbor] + mur_coeff * (new_neighbor_value - bz[i][j][k]);
            }
        }
    }

    //! z上部境界
    if (is_boundary[5]) {
        int k = cz_with_glue - 2; int k_neighbor = k - 1;
        for(int i = 1; i < cx_with_glue - 1; ++i){
            for(int j = 1; j < cy_with_glue - 1; ++j){
                const double new_neighbor_value = bz[i][j][k_neighbor] + dt_per_dx * ((ex[i][j + 1][k_neighbor] - ex[i][j][k_neighbor]) - (ey[i + 1][j][k_neighbor] - ey[i][j][k_neighbor]));

                bz[i][j][k] = bz[i][j][k_neighbor] + mur_coeff * (new_neighbor_value - bz[i][j][k]);
            }
        }
    }
    */

    for(int i = 1; i < cx_with_glue - 1; ++i){
        for(int j = 1; j < cy_with_glue - 1; ++j){
            for(int k = 1; k < cz_with_glue - 1; ++k){
                if (j != cy_with_glue - 2 && k != cz_with_glue - 2) {
                    if (!( (i == 1 && is_boundary[0]) && (i == cx_with_glue - 2 && is_boundary[1]) )) {
                        bx[i][j][k] = bx[i][j][k] + dt_per_dx * ((ey[i][j][k + 1] - ey[i][j][k]) - (ez[i][j + 1][k] - ez[i][j][k]));
                    }
                }

                if (i != cx_with_glue - 2 && k != cz_with_glue - 2) {
                    if (!( (j == 1 && is_boundary[2]) && (j == cy_with_glue - 2 && is_boundary[3]) )) {
                        by[i][j][k] = by[i][j][k] + dt_per_dx * ((ez[i + 1][j][k] - ez[i][j][k]) - (ex[i][j][k + 1] - ex[i][j][k]));
                    }
                }

                if (i != cx_with_glue - 2 && j != cy_with_glue - 2) {
                    if (!( (k == 1 && is_boundary[4]) && (k == cz_with_glue - 2 && is_boundary[5]) )) {
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

//! Reference Current (ノード上で定義される電流) を更新する
//! 計算には直接必要ないが、データ出力のため
void RootGrid::updateReferenceCurrent() {
    auto& jx = field->getJx();
    auto& jy = field->getJy();
    auto& jz = field->getJz();
    auto& jxref = field->getJxRef();
    auto& jyref = field->getJyRef();
    auto& jzref = field->getJzRef();
    const size_t cx_with_glue = jx.shape()[0] + 1; // nx + 2
    const size_t cy_with_glue = jy.shape()[1] + 1;
    const size_t cz_with_glue = jz.shape()[2] + 1;

    //! reference 更新
    #pragma omp parallel
    {
        #pragma omp for
        for(size_t i = 1; i < cx_with_glue - 1; ++i){
            for(size_t j = 1; j < cy_with_glue - 1; ++j){
                for(size_t k = 1; k < cz_with_glue - 1; ++k){
                    jxref[i][j][k] = 0.5 * (jx[i-1][j][k] + jx[i][j][k]);
                    jyref[i][j][k] = 0.5 * (jy[i][j-1][k] + jy[i][j][k]);
                    jzref[i][j][k] = 0.5 * (jz[i][j][k-1] + jz[i][j][k]);
                }
            }
        }
    }

    MPIw::Environment::sendRecvField(jxref);
    MPIw::Environment::sendRecvField(jyref);
    MPIw::Environment::sendRecvField(jzref);
}
