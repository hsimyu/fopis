#include "grid.hpp"
#include "field.hpp"
#include "utils.hpp"
#include "normalizer.hpp"

//! child grid constructor
//! 渡されたGridを親とした子グリッドを生成する
ChildGrid::ChildGrid(std::shared_ptr<Grid> g, const int _from_ix, const int _from_iy, const int _from_iz, const int _to_ix,   const int _to_iy,   const int _to_iz) : Grid() {
    constexpr double refineRatio = 2.0;

    parent = g;
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

int ChildGrid::getXNodeSize(void) const { return nx; }
int ChildGrid::getYNodeSize(void) const { return ny; } 
int ChildGrid::getZNodeSize(void) const { return nz; }

void ChildGrid::solvePoisson(void) {
    constexpr int PRE_LOOP_NUM = 10;
    constexpr int POST_LOOP_NUM = 250;

    if (this->getChildrenLength() > 0) {
        this->solvePoissonPSOR(PRE_LOOP_NUM);
        this->updateChildrenPhi();

        for(auto& child : children) {
            child->solvePoisson();
        }
    }

    this->solvePoissonPSOR(POST_LOOP_NUM);
    this->copyPhiToParent();

    if (this->getChildrenLength() > 0) {
        this->correctChildrenPhi();
    }
}

void ChildGrid::solvePoissonPSOR(const int loopnum) {
    auto& phi = field->getPhi();
    auto& poisson_error = field->getPoissonError();
    auto& rho = field->getRho();

    const double omega = 2.0/(1.0 + sin(M_PI/(phi.shape()[0] - 2))); // spectral radius
    const double rho_coeff = pow(dx, 2) / Normalizer::eps0;

    const int cx_with_glue = phi.shape()[0];
    const int cy_with_glue = phi.shape()[1];
    const int cz_with_glue = phi.shape()[2];

    constexpr double required_error = 1.0e-7;

    poisson_error = phi;

    for(int loop = 1; loop <= loopnum; ++loop) {
        for(int k = 2; k < cz_with_glue - 2; ++k){
            for(int j = 2; j < cy_with_glue - 2; ++j){
                for(int i = 2; i < cx_with_glue - 2; ++i){
                    phi[i][j][k] = (1.0 - omega) * phi[i][j][k] + omega*(phi[i+1][j][k] + phi[i-1][j][k] + phi[i][j+1][k] + phi[i][j-1][k] + phi[i][j][k+1] + phi[i][j][k-1] + rho_coeff * rho[0][i][j][k])/6.0;
                }
            }
        }

        if ( (loop % 10 == 0) && (this->checkPhiResidual() < required_error) ) {
            cout << "performed " << loop << " iterations." << endl;
            break;
        }
    }

    //! 全グリッド上のエラーを更新
    for(int k = 0; k < cz_with_glue; ++k){
        for(int j = 0; j < cy_with_glue; ++j){
            for(int i = 0; i < cx_with_glue; ++i){
                poisson_error[i][j][k] = phi[i][j][k] - poisson_error[i][j][k];
            }
        }
    }
}

double ChildGrid::checkPhiResidual() {
    double residual = 0.0;
    const double normalized_eps = Normalizer::eps0;
    const double per_dx2 = 1.0 / pow(dx, 2);

    auto& phi = field->getPhi();
    auto& rho = field->getRho();
    auto& poisson_residual = field->getPoissonResidual();

    const int cx_with_glue = phi.shape()[0];
    const int cy_with_glue = phi.shape()[1];
    const int cz_with_glue = phi.shape()[2];

    for(int k = 2; k < cz_with_glue - 2; ++k){
        for(int j = 2; j < cy_with_glue - 2; ++j){
            for(int i = 2; i < cx_with_glue - 2; ++i){
                double source_value = rho[0][i][j][k]/normalized_eps;
                double tmp_res = field->poissonOperator(phi, i, j, k) * per_dx2 + source_value;
                poisson_residual[i][j][k] = tmp_res;

                if (fabs(tmp_res / source_value) > residual) {
                    residual = fabs(tmp_res / source_value);
                }
            }
        }
    }
    return residual;
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
        //! 親へ送る
        pbuff[0].push_back(p);
        p.makeInvalid();
    } else if (p.x > slx) {
        pbuff[1].push_back(p);
        p.makeInvalid();
    }
}

void ChildGrid::checkYBoundary(ParticleArray& pbuff, Particle& p, const double sly) {
    if(p.y < 0.0) {
        pbuff[2].push_back(p);
        p.makeInvalid();
    } else if (p.y > sly) {
        pbuff[3].push_back(p);
        p.makeInvalid();
    }
}

void ChildGrid::checkZBoundary(ParticleArray& pbuff, Particle& p, const double slz) {
    if(p.z < 0.0) {
        pbuff[4].push_back(p);
        p.makeInvalid();
    } else if (p.z > slz) {
        pbuff[5].push_back(p);
        p.makeInvalid();
    }
}

void ChildGrid::updateRho() {}

void ChildGrid::copyPhiToParent(){
    tdArray& phi = field->getPhi();
    tdArray& residual = field->getPoissonResidual();
    tdArray& parentPhi = parent->getPhi();
    RhoArray& parentRho = parent->getRho();
    const double per_dx2 = 1.0 / pow(parent->getDx(), 2);
    const double rho_coeff = -Normalizer::eps0;

    for(int ix = from_ix; ix <= to_ix; ++ix){
        int i = 2 * (ix - from_ix) + 1;
        for(int iy = from_iy; iy <= to_iy; ++iy){
            int j = 2 * (iy - from_iy) + 1;
            for(int iz = from_iz; iz <= to_iz; ++iz){
                int k = 2 * (iz - from_iz) + 1;
                //! まず v^h -> v^2h にコピーする
                parentPhi[ix][iy][iz] = phi[i][j][k];
            }
        }
    }
    for(int ix = from_ix + 1; ix < to_ix; ++ix){
        int i = 2 * (ix - from_ix) + 1;
        for(int iy = from_iy + 1; iy < to_iy; ++iy){
            int j = 2 * (iy - from_iy) + 1;
            for(int iz = from_iz + 1; iz < to_iz; ++iz){
                int k = 2 * (iz - from_iz) + 1;
                //! ソース項が A^2h v^2h + r^2h となるように
                //! PoissonOperatorを適用し、現在のGridのResidualを足したものに-eps0をかけて新しいrhoとする
                parentRho[0][ix][iy][iz] = rho_coeff * (field->poissonOperator(parentPhi, ix, iy, iz) * per_dx2 + residual[i][j][k]);
            }
        }
    }
}

void ChildGrid::decrementSumOfChild() {
    --sumTotalNumOfChildGrids;
    parent->decrementSumOfChild();
}

void ChildGrid::incrementSumOfChild() {
    ++sumTotalNumOfChildGrids;
    parent->incrementSumOfChild();
}
