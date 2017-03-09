#include "global.hpp"
#include "environment.hpp"
#include "field.hpp"
#include "utils.hpp"
#include "mpiw.hpp"
#include <algorithm>

void Field::initializePoisson(const int cx, const int cy, const int cz){
#ifdef DEBUG
    if( Environment::isRootNode ) cout << "--  Poisson Initializing  --" << endl;
#endif

    psn = new Poisson();
    // mesh interval なので -1 する
    int nx = cx + 2 - 1;
    int ny = cy + 2 - 1;
    int nz = cz + 2 - 1;

    psn->nx = nx;
    psn->ny = ny;
    psn->nz = nz;

    psn->ipar = new MKL_INT[128];
    for(int i = 0; i < 128; ++i) {
        psn->ipar[i] = 0;
    }

    psn->dpar = new double[5*(nx*ny)/2 +9];

    double zero = 0.0;
    double upper_x = cx;
    double q = 0.0; // 0 for Poisson and Laplace

    d_init_Helmholtz_3D(&zero, &upper_x, &zero, &upper_x, &zero, &upper_x, &nx, &ny, &nz, Environment::boundary.c_str(), &q, psn->ipar, psn->dpar, &(psn->stat));

    if(psn->stat != 0) cout << Environment::rankStr() << format("[Poisson Init] stat == %d") % psn->stat << endl;

    psn->b_lx = new double[(ny + 1) * (nz + 1)];
    for(int i = 0; i < (ny +1) * (nz + 1); ++i) psn->b_lx[i] = 0.0;

    psn->b_ly = new double[(nx + 1) * (nz + 1)];
    for(int i = 0; i < (nx +1) * (nz + 1); ++i) psn->b_ly[i] = 0.0;

    psn->b_lz = new double[(nx + 1) * (ny + 1)];
    for(int i = 0; i < (nx +1) * (ny + 1); ++i) psn->b_lz[i] = 0.0;

    psn->xhandle = new DFTI_DESCRIPTOR_HANDLE();
    psn->yhandle = new DFTI_DESCRIPTOR_HANDLE();

    psn->rho1D = new double[(nx + 1)*(ny + 1)*(nz + 1)];

#ifdef DEBUG
    if( Environment::isRootNode ) cout << "--  End Poisson Initializing  --" << endl;
#endif
}

//! @brief MKLのPoisson Solverを呼び出す
//! もしPoisson用構造体のポインタがnullptrだった場合、初期化を行う
void Field::solvePoissonMKL(const int cx, const int cy, const int cz) {
    if(psn == nullptr) initializePoisson(cx, cy, cz);

    //! space charge を phi 配列に copy
    //! multi_array& operator=(const multi_array& x);
    //! "This performs an element-wise copy of x into the current multi_array."
    phi = rho;

    //! Commit solver
    //! @note Is it required?
    d_commit_Helmholtz_3D(phi.data(), psn->b_lx, psn->b_lx, psn->b_ly, psn->b_ly, psn->b_lz, psn->b_lz, psn->xhandle, psn->yhandle, psn->ipar, psn->dpar, &(psn->stat));
    if(psn->stat != 0) cout << Environment::rankStr() << format("[Poisson Commit] stat == %d") % psn->stat << endl;

    //! rho(1D) -> phi(1D)
    d_Helmholtz_3D(phi.data(), psn->b_lx, psn->b_lx, psn->b_ly, psn->b_ly, psn->b_lz, psn->b_lz, psn->xhandle, psn->yhandle, psn->ipar, psn->dpar, &(psn->stat));
    if(psn->stat != 0) cout << Environment::rankStr() << format("[Poisson Solve] stat == %d") % psn->stat << endl;
}

//! @brief SOR法
void Field::solvePoissonPSOR(const int loopnum, const double dx) {
    double omega = 2.0/(1.0 + sqrt(1.0 - pow(cos(M_PI/(phi.shape()[0] - 2)), 2))); // spectral radius
    double rho_coeff = 6.0 * dx / Utils::Normalizer::normalizeEpsilon(eps0);

    const int cx_with_glue = phi.shape()[0];
    const int cy_with_glue = phi.shape()[1];
    const int cz_with_glue = phi.shape()[2];

    this->setDirichletPhi();

    const double required_error = 1.0e-3;

    for(int loop = 1; loop <= loopnum; ++loop) {
        for(int k = 1; k < cz_with_glue - 1; ++k){
            if((k != 1 || !Environment::onLowZedge) && (k != cz_with_glue - 2 || !Environment::onHighZedge)) {
                for(int j = 1; j < cy_with_glue - 1; ++j){
                    if((j != 1 || !Environment::onLowYedge) && (j != cy_with_glue - 2 || !Environment::onHighYedge)) {
                        for(int i = 1; i < cx_with_glue - 1; ++i){
                            if((i != 1 || !Environment::onLowXedge) && (i != cx_with_glue - 2 || !Environment::onHighXedge)) {
                                phi[i][j][k] = (1.0 - omega) * phi[i][j][k] + omega*(phi[i+1][j][k] + phi[i-1][j][k] + phi[i][j+1][k] + phi[i][j-1][k] + phi[i][j][k+1] + phi[i][j][k-1] + rho_coeff*rho[i][j][k])/6.0;
                            }
                        }
                    }
                }
            }
        }

        if(MPIw::Environment::numprocs > 1) {
            //! @note: 実際には部分部分をSORで計算して送信というのを繰り返す方が収束効率がよい
            MPIw::Environment::sendRecvField(phi);
        }

        if(loop % 10 == 0) {
            double residual = this->checkPhiResidual();

            if (Environment::isRootNode) cout << format("[ITR %d] residual = %s") % loop % residual << endl;
            if (residual < required_error) break;
        }
    }
}

//! @brief Jacobi法
void Field::solvePoissonJacobi(const int loopnum, const double dx) {
    double rho_coeff = dx / Utils::Normalizer::normalizeEpsilon(eps0);

    const int cx_with_glue = phi.shape()[0];
    const int cy_with_glue = phi.shape()[1];
    const int cz_with_glue = phi.shape()[2];

    this->setDirichletPhi();

    cout << format("initial residual = %s") % this->checkPhiResidual() << endl;
    const double required_error = 1.0e-3;

    for(int loop = 0; loop < loopnum; ++loop) {
        for(int k = 1; k < cz_with_glue - 1; ++k){
            if((k != 1 || !Environment::onLowZedge) && (k != cz_with_glue - 2 || !Environment::onHighZedge)) {
                for(int j = 1; j < cy_with_glue - 1; ++j){
                    if((j != 1 || !Environment::onLowYedge) && (j != cy_with_glue - 2 || !Environment::onHighYedge)) {
                        for(int i = 1; i < cx_with_glue - 1; ++i){
                            if((i != 1 || !Environment::onLowXedge) && (i != cx_with_glue - 2 || !Environment::onHighXedge)) {
                                phi[i][j][k] = (phi[i+1][j][k] + phi[i-1][j][k] + phi[i][j+1][k] + phi[i][j-1][k] + phi[i][j][k+1] + phi[i][j][k-1])/6.0 + rho_coeff*rho[i][j][k];
                            }
                        }
                    }
                }
            }
        }

        double residual = this->checkPhiResidual();

        cout << format("[ITR %d] residual = %s") % loop % residual << endl;
        if (residual < required_error) break;
        //! ここで通信する必要がある
        // MPIw::Environment::sendRecvPhi(phi);
        //! @note: 実際には部分部分をSORで計算して送信というのを繰り返す方が収束効率がよい
    }
}

void Field::setDirichletPhi(void){
    const int cx_with_glue = phi.shape()[0];
    const int cy_with_glue = phi.shape()[1];
    const int cz_with_glue = phi.shape()[2];

    //! Dirichlet境界条件に値を設定 (V=0V)
    if(Environment::onLowZedge) {
        for(int j = 1; j < cy_with_glue - 1; ++j){
            for(int i = 1; i < cx_with_glue - 1; ++i){
                phi[i][j][1] = 0.0;
            }
        }
    }

    if(Environment::onHighZedge) {
        int k = cz_with_glue - 2;
        for(int j = 1; j < cy_with_glue - 1; ++j){
            for(int i = 1; i < cx_with_glue - 1; ++i){
                phi[i][j][k] = 0.0;
            }
        }
    }

    if(Environment::onLowYedge) {
        for(int k = 1; k < cz_with_glue - 1; ++k){
            for(int i = 1; i < cx_with_glue - 1; ++i){
                phi[i][1][k] = 0.0;
            }
        }
    }

    if(Environment::onHighYedge) {
        int j = cy_with_glue - 2;
        for(int k = 1; k < cz_with_glue - 1; ++k){
            for(int i = 1; i < cx_with_glue - 1; ++i){
                phi[i][j][k] = 0.0;
            }
        }
    }

    if(Environment::onLowXedge) {
        for(int k = 1; k < cz_with_glue - 1; ++k){
            for(int j = 1; j < cy_with_glue - 1; ++j){
                phi[1][j][k] = 0.0;
            }
        }
    }

    if(Environment::onHighXedge) {
        int i = cx_with_glue - 2;
        for(int k = 1; k < cz_with_glue - 1; ++k){
            for(int j = 1; j < cy_with_glue - 1; ++j){
                phi[i][j][k] = 0.0;
            }
        }
    }
}

//! Poissonソルバを呼び出す
void Field::solvePoisson(const int loopnum, const double dx) {
    // this->solvePoissonMKL(cx, cy, cz);
    this->solvePoissonPSOR(loopnum, dx);
    // this->solvePoissonJacobi(loopnum, dx);
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
                            double tmp_res = (phi[i-1][j][k] + phi[i+1][j][k] + phi[i][j-1][k] + phi[i][j+1][k] + phi[i][j][k-1] + phi[i][j][k+1] - 6.0*phi[i][j][k])/6.0 + rho[i][j][k]/normalized_eps;
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

//! @brief 電場を更新する
//! e = - (p_+1 - p_+0)/dx
void Field::updateEfield(const double dx) {
    const int cx_with_glue = ex.shape()[0] + 1; // nx + 2
    const int cy_with_glue = ey.shape()[1] + 1;
    const int cz_with_glue = ez.shape()[2] + 1;
    const double dxm = 1.0/dx;

    //! phiは通信してあるとする -> 0番目のedgeも計算可能
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

    //! @note:隣と通信しなくてもいい？？

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

    if(MPIw::Environment::numprocs > 1) {
        MPIw::Environment::sendRecvField(exref);
        MPIw::Environment::sendRecvField(eyref);
        MPIw::Environment::sendRecvField(ezref);
    }
}

//! @brief 磁場を更新する(FDTD)
//!
void Field::updateBfield(const double dx, const int nx, const int ny, const int nz) {
    const double dt = Utils::Normalizer::normalizeTime(Environment::dt);

    const double mu0 = 4.0 * M_PI * 1e-7; // H/m
    const double epsilon_r = 1.0; //! 比誘電率
    const double sigma = 1.0; //! 導電率 (各Faceでの)
    const double mu_r = 1.0; //! 透磁率 (各Faceでの)
    const double sigma_m = 0.0; //! 導磁率?
    //
    // // 電場への係数
    // const double c1 = epsilon_r/(epsilon_r + eta * sigma * c * dt);
    // const double c2 = 1.0/(epsilon_r + eta * sigma * c * dt);

    // 磁束密度更新時の係数
    const double d1 = mu_r/(mu_r + sigma_m * dt / mu0);
    const double d2 = dt/(mu_r + sigma_m * dt / mu0);

    // cout << "c1 = " << c1 << ", c2 = " << c2 << endl;
    cout << "d1 = " << d1 << ", d2 = " << d2 << endl;

    const int cx_with_glue = nx + 1;
    const int cy_with_glue = nx + 1;
    const int cz_with_glue = nx + 1;

    //! 0とcy + 1, 0とcz + 1はglueなので更新しなくてよい
    for(int i = 1; i < cx_with_glue; ++i){
        for(int j = 1; j < cy_with_glue; ++j){
            for(int k = 1; k < cz_with_glue; ++k){
                if(j != cy_with_glue - 1 && k != cz_with_glue - 1)
                    bx[i][j][k] = d1 * bx[i][j][k] - d2 * (ez[i][j + 1][k] - ez[i][j][k] - ey[i][j][k + 1] + ey[i][j][k]) / dx;
                if(i != cx_with_glue - 1 && k != cz_with_glue - 1)
                    by[i][j][k] = d1 * by[i][j][k] - d2 * (ex[i][j][k + 1] - ex[i][j][k] - ez[i + 1][j][k] + ez[i][j][k]) / dx;
                if(i != cx_with_glue - 1 && j != cy_with_glue - 1)
                    bz[i][j][k] = d1 * bz[i][j][k] - d2 * (ey[i + 1][j][k] - ey[i][j][k] - ex[i][j + 1][k] + ex[i][j][k]) / dx;
            }
        }
    }

}

double Field::getEnergy(const int nx, const int ny, const int nz) {
    double energy = 0.0;
    double normalized_eps0 = Utils::Normalizer::normalizeEpsilon(eps0);

    for(int i = 1; i < nx - 1; ++i){
        for(int j = 1; j < ny - 1; ++j){
            for(int k = 1; k < nz - 1; ++k){
                //! 各点のエネルギーを計算する(Yee格子内のエネルギーの計算方法は?)
                energy += pow(exref[i][j][k],2) + pow(eyref[i][j][k], 2) + pow(ezref[i][j][k], 2);
            }
        }
    }

    return 0.5 * Utils::Normalizer::normalizeEpsilon(eps0) * energy;
}

