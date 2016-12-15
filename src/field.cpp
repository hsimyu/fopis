#include <tdpic.h>

// potential
void Field::setPhi(threeDArray _phi){
    phi = _phi;
}
threeDArray Field::getPhi(){
    return phi;
}

// charge density
void Field::setRho(threeDArray _rho){
    rho = _rho;
}
threeDArray Field::getRho(){
    return rho;
}

// electric fields
void Field::setEx(threeDArray _ex){
    ex = _ex;
}
threeDArray Field::getEx(){
    return ex;
}

void Field::setEy(threeDArray _ey){
    ey = _ey;
}
threeDArray Field::getEy(){
    return ey;
}

void Field::setEz(threeDArray _ez){
    ez = _ez;
}
threeDArray Field::getEz(){
    return ez;
}

// magnetic fields
void Field::setBx(threeDArray _bx){
    bx = _bx;
}
threeDArray Field::getBx(){
    return bx;
}

void Field::setBy(threeDArray _by){
    by = _by;
}
threeDArray Field::getBy(){
    return by;
}

void Field::setBz(threeDArray _bz){
    bz = _bz;
}
threeDArray Field::getBz(){
    return bz;
}

void Field::initializePoisson(const Environment* env){

#ifdef DEBUG
    std::cout << "--  Poisson Initializing  --" << std::endl;
#endif

    psn = new Poisson();
    // mesh interval なので -1 する
    int nx = env->cell_x + 2 - 1;
    int ny = env->cell_y + 2 - 1;
    int nz = env->cell_z + 2 - 1;

    psn->nx = nx;
    psn->ny = ny;
    psn->nz = nz;

    psn->ipar = new MKL_INT[128];
    for(int i = 0; i < 128; ++i) {
        psn->ipar[i] = 0;
    }

    psn->dpar = new double[5*(nx*ny)/2 +9];

    double zero = 0.0;
    double upper_x = env->cell_x;
    double q = 0.0; // 0 for Poisson and Laplace

    d_init_Helmholtz_3D(&zero, &upper_x, &zero, &upper_x, &zero, &upper_x, &nx, &ny, &nz, env->boundary.c_str(), &q, psn->ipar, psn->dpar, &(psn->stat));
    if(psn->stat != 0) std::cout << "stat == " << psn->stat << std::endl;

    psn->b_lx = new double[(ny + 1) * (nz + 1)];
    for(int i = 0; i < (ny +1) * (nz + 1); ++i) psn->b_lx[i] = 0.0;

    psn->b_ly = new double[(nx + 1) * (nz + 1)];
    for(int i = 0; i < (nx +1) * (nz + 1); ++i) psn->b_ly[i] = 0.0;

    psn->b_lz = new double[(nx + 1) * (ny + 1)];
    for(int i = 0; i < (nx +1) * (ny + 1); ++i) psn->b_lz[i] = 0.0;

    psn->xhandle = new DFTI_DESCRIPTOR_HANDLE();
    psn->yhandle = new DFTI_DESCRIPTOR_HANDLE();

    psn->rho1D = new double[(nx + 1)*(ny + 1)*(nz + 1)];

    //! rho(3D) -> rho(1D)
#ifdef DEBUG
    std::cout << "--  End Poisson Initializing  --" << std::endl;
#endif
}

//! @brief MKLのPoisson Solverを呼び出す
//! もしPoisson用構造体のポインタがnullptrだった場合、初期化を行う
//!
//! @param Environment
void Field::solvePoisson(const Environment* env) {
    if(psn == nullptr) initializePoisson(env);

    //! rho(3D) -> rho(1D)
    Utils::convert3Dto1Darray(rho, psn->nx + 1, psn->ny + 1, psn->nz + 1, psn->rho1D);

    //! Commit solver
    //! @note Is it required?
    d_commit_Helmholtz_3D(psn->rho1D, psn->b_lx, psn->b_lx, psn->b_ly, psn->b_ly, psn->b_lz, psn->b_lz, psn->xhandle, psn->yhandle, psn->ipar, psn->dpar, &(psn->stat));
    if(psn->stat != 0) std::cout << "stat == " << psn->stat << std::endl;

    //! rho(1D) -> phi(1D)
    d_Helmholtz_3D(psn->rho1D, psn->b_lx, psn->b_lx, psn->b_ly, psn->b_ly, psn->b_lz, psn->b_lz, psn->xhandle, psn->yhandle, psn->ipar, psn->dpar, &(psn->stat));
    if( psn->stat != 0) std::cout << "stat == " << psn->stat << std::endl;

    //! phi(1D) -> phi(3D)
    Utils::convert1Dto3Darray(psn->rho1D, psn->nx + 1, psn->ny + 1, psn->nz + 1, phi);
}

//! @brief 電場を更新する
//! e = - (p_+1 - p_+0)/dx
//! @param const Environment*
void Field::updateEfield(const Environment* env) {

    const int nx = env->cell_x + 1;
    const int ny = env->cell_y + 1;
    const int nz = env->cell_z + 1;

    //! 0とcy + 1, 0とcz + 1はglueなので更新しなくてよい
    for(int i = 1; i < nx; ++i){
        for(int j = 1; j < ny; ++j){
            for(int k = 1; k < nz; ++k){
                //! 各方向には1つ少ないのでcx-1まで
                if(i < nx - 1) ex[i][j][k] = phi[i][j][k] - phi[i + 1][j][k];
                if(j < ny - 1) ey[i][j][k] = phi[i][j][k] - phi[i][j + 1][k];
                if(k < nz - 1) ez[i][j][k] = phi[i][j][k] - phi[i][j][k + 1];
            }
        }
    }
}

// destructor
Field::~Field(){
    Utils::delete3DArray(phi);
    Utils::delete3DArray(rho);
    Utils::delete3DArray(ex);
    Utils::delete3DArray(ey);
    Utils::delete3DArray(ez);
    Utils::delete3DArray(bx);
    Utils::delete3DArray(by);
    Utils::delete3DArray(bz);

    delete psn;
}

Poisson::~Poisson(){
    delete [] ipar;
    delete [] dpar;
    delete [] b_lx;
    delete [] b_ly;
    delete [] b_lz;
    delete [] xhandle;
    delete [] yhandle;
    delete [] rho1D;
}
