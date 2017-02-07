#include "global.hpp"
#include "environment.hpp"
#include "field.hpp"
#include "utils.hpp"

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
void Field::solvePoissonSOR(const int loopnum, const double dx) {
    double omega = 2.0/(1.0 + M_PI/phi.shape()[0]); // spectral radius
    double rho_coeff = 6.0 * dx / Utils::Normalizer::normalizeEpsilon(eps0);

    for(int loop = 0; loop < loopnum; ++loop) {
        for(int k = 1; k < phi.shape()[2] - 1; ++k){
            for(int j = 1; j < phi.shape()[1] - 1; ++j){
                for(int i = 1; i < phi.shape()[0] - 1; ++i){
                    phi[i][j][k] = (1.0 - omega) * phi[i][j][k] + omega*(phi[i+1][j][k] + phi[i-1][j][k] + phi[i][j+1][k] + phi[i][j-1][k] + phi[i][j][k+1] + phi[i][j][k-1] + rho_coeff*rho[i][j][k])/6.0;
                }
            }
        }
    }
}

//! Poissonソルバを呼び出す
void Field::solvePoisson(const int loopnum, const double dx) {
    // this->solvePoissonMKL(cx, cy, cz);
    this->solvePoissonSOR(loopnum, dx);
}

//! @brief 電場を更新する
//! e = - (p_+1 - p_+0)/dx
void Field::updateEfield(const int nx, const int ny, const int nz) {
    const int nx_with_glue = nx + 1; // (cx + 2) - 1
    const int ny_with_glue = ny + 1;
    const int nz_with_glue = nz + 1;

    //! phiは通信してあるとする -> 0番目のedgeも計算可能
    for(int i = 0; i < nx_with_glue; ++i){
        for(int j = 0; j < ny_with_glue; ++j){
            for(int k = 0; k < nz_with_glue; ++k){
                //! 各方向には1つ少ないのでcx-1まで
                if(i < nx_with_glue - 1) ex[i][j][k] = phi[i][j][k] - phi[i + 1][j][k];
                if(j < ny_with_glue - 1) ey[i][j][k] = phi[i][j][k] - phi[i][j + 1][k];
                if(k < nz_with_glue - 1) ez[i][j][k] = phi[i][j][k] - phi[i][j][k + 1];
            }
        }
    }

    //! @note:隣と通信しなくてもいい？？

    //! reference 更新
    for(int i = 1; i < nx_with_glue; ++i){
        for(int j = 1; j < ny_with_glue; ++j){
            for(int k = 1; k < nz_with_glue; ++k){
                exref[i][j][k] = 0.5 * (ex[i-1][j][k] + ex[i][j][k]);
                eyref[i][j][k] = 0.5 * (ey[i][j-1][k] + ey[i][j][k]);
                ezref[i][j][k] = 0.5 * (ez[i][j][k-1] + ez[i][j][k]);
            }
        }
    }
}

//! @brief 磁場を更新する(FDTD)
//!
void Field::updateBfield(const int nx, const int ny, const int nz) {
    const double dx = Environment::dx;
    const double dt = Environment::dt;
    //const double DT = 0.1 * DX/c
    const double eta = 1.0/(c * eps0);
    const double epsilon_r = 1.0;
    const double mu_r = 1.0;
    const double sigma = 1.0;
    const double sigma_m = 1.0;
    const double c1 = epsilon_r/(epsilon_r + eta * sigma * c * dt);
    const double c2 = 1.0/(epsilon_r + eta * sigma * c * dt);
    const double d1 = mu_r/(mu_r + sigma_m * c * dt / eta);
    const double d2 = 1.0/(mu_r + sigma_m * c * dt / eta);

    cout << "c1 = " << c1 << ", c2 = " << c2 << endl;
    cout << "d1 = " << d1 << ", d2 = " << d2 << endl;
    cout << "eta = " << eta << endl;

    const int cx_with_glue = nx + 1;
    const int cy_with_glue = nx + 1;
    const int cz_with_glue = nx + 1;

    //! 0とcy + 1, 0とcz + 1はglueなので更新しなくてよい
    for(int i = 1; i < cx_with_glue; ++i){
        for(int j = 1; j < cy_with_glue; ++j){
            for(int k = 1; k < cz_with_glue; ++k){
                // bx[i][j][k] = d1 * bx[i][j][k] - d2 * ((ez[i][j + 1][k] - ez[i][j][k]) / (dx/(c * dt)) - (ey[i][j][k + 1] - ey[i][j][k]) / (dx/(c * dt))) / eta;
                // by[i][j][k] = d1 * by[i][j][k] - d2 * ((ex[i][j][k + 1] - ex[i][j][k]) / (dx/(c * dt)) - (ez[i + 1][j][k] - ez[i][j][k]) / (dx/(c * dt))) / eta;
                // bz[i][j][k] = d1 * bz[i][j][k] - d2 * ((ey[i + 1][j][k] - ey[i][j][k]) / (dx/(c * dt)) - (ex[i][j + 1][k] - ex[i][j][k]) / (dx/(c * dt))) / eta;
                // hx[i][j][k] = d1 * hx[i][j][k] - d2 * ((ez[i][j + 1][k] - ez[i][j][k]) / (DY/(c * DT)) - (ey[i][j][k + 1] - ey[i][j][k]) / (DZ/(c * DT))) / eta
                // hy[i][j][k] = d1 * hy[i][j][k] - d2 * ((ex[i][j][k + 1] - ex[i][j][k]) / (DZ/(c * DT)) - (ez[i + 1][j][k] - ez[i][j][k]) / (DX/(c * DT))) / eta
                // hz[i][j][k] = d1 * hz[i][j][k] - d2 * ((ey[i + 1][j][k] - ey[i][j][k]) / (DX/(c * DT)) - (ex[i][j + 1][k] - ex[i][j][k]) / (DY/(c * DT))) / eta
            }
        }
    }

}

