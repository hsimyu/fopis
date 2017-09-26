#include "grid.hpp"
#include "field.hpp"
#include "utils.hpp"
#include "normalizer.hpp"

#define USE_BOOST
#include "simple_vtk.hpp"

//! child grid constructor
//! 渡されたGridを親とした子グリッドを生成する
ChildGrid::ChildGrid(Grid* g, const int _from_ix, const int _from_iy, const int _from_iz, const int _to_ix,   const int _to_iy,   const int _to_iz) : Grid() {
    constexpr double refineRatio = 2.0;

    parent = g;
    level = g->getLevel() + 1;

    //! UniqueなIDをセット
    constexpr int minimum_id_offset = 10;
    id = minimum_id_offset * MPIw::Environment::rank + this->getNextID();

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

    this->initializeChildMap();

    // refineRatioは2で固定
    dx = g->getDx() / refineRatio;
    dt = g->getDt() / refineRatio;

    this->checkGridValidness();

    // Field初期化
    this->initializeField();

    //! - particle types 分だけresize
    particles.resize(Environment::num_of_particle_types);
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
        this->restrictPhiToChildrenBoundary();

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
            if (Environment::isRootNode) {
                cout << "performed " << loop << " iterations." << endl;
            }
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

//! @brief 差分法で電場を更新する
//! e = - (p_+1 - p_+0)/dx
void ChildGrid::updateEfield() {
    auto& ex = field->getEx();
    auto& ey = field->getEy();
    auto& ez = field->getEz();
    auto& phi = field->getPhi();

    const int cx_with_glue = ex.shape()[0] + 1; // nx + 2
    const int cy_with_glue = ey.shape()[1] + 1;
    const int cz_with_glue = ez.shape()[2] + 1;
    const double per_dx = 1.0 / dx;

    //! @note:隣と通信しなくてもいい
    //! phiが通信してあるため、端の要素を通信なしで計算可能
    for(int i = 0; i < cx_with_glue; ++i){
        for(int j = 0; j < cy_with_glue; ++j){
            for(int k = 0; k < cz_with_glue; ++k){
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
void ChildGrid::updateReferenceEfield() {
    auto& ex = field->getEx();
    auto& ey = field->getEy();
    auto& ez = field->getEz();
    auto& exref = field->getExRef();
    auto& eyref = field->getEyRef();
    auto& ezref = field->getEzRef();
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
}

void ChildGrid::updateEfieldFDTD(void) {
    field->updateEfieldFDTD(dx, dt);
}

void ChildGrid::updateBfield(void) {
    field->updateBfield(dx, nx, ny, nz, dt);
}

void ChildGrid::moveParticleToParent(Particle& p) {
    Particle new_particle = p; // コピー演算

    new_particle.x = (new_particle.x / 2.0 + static_cast<double>(from_ix) - 1.0);
    new_particle.y = (new_particle.y / 2.0 + static_cast<double>(from_iy) - 1.0);
    new_particle.z = (new_particle.z / 2.0 + static_cast<double>(from_iz) - 1.0);

    //! 親グリッド上の粒子をinvalidに
    p.makeInvalid();
    //! 子グリッド上の粒子をpush
    parent->addParticle( std::move(new_particle) );
}

//! 粒子の位置から電荷を空間電荷にする
void ChildGrid::updateRho() {
    RhoArray& rho = field->getRho();

    //! rhoを初期化
    Utils::initializeRhoArray(rho);

    for(int pid = 0; pid < Environment::num_of_particle_types; ++pid){
        double q = Environment::getParticleType(pid)->getChargeOfSuperParticle() * 8.0;
        const auto rho_idx = pid + 1;

        for(auto& p : particles[pid]) {
            if(p.isValid) {
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

    for(int pid = 1; pid < Environment::num_of_particle_types + 1; ++pid) {
        for(int i = 1; i < nx + 2; ++i) {
            for(int j = 1; j < ny + 2; ++j) {
                for(int k = 1; k < nz + 2; ++k) {
                    rho[0][i][j][k] += rho[pid][i][j][k];
                }
            }
        }
    }

    //! 子の電荷更新
    for(auto& child : children) {
        child->updateRho();
    }
}

void ChildGrid::copyPhiToParent(){
    tdArray& phi = field->getPhi();
    tdArray& residual = field->getPoissonResidual();
    tdArray& parentPhi = parent->getPhi();
    RhoArray& parentRho = parent->getRho();
    const double per_dx2 = 1.0 / pow(parent->getDx(), 2);
    const double rho_coeff = -Normalizer::eps0;
    constexpr double div = 1.0 / 12.0;

    for(int ix = from_ix; ix <= to_ix; ++ix){
        int i = 2 * (ix - from_ix) + 1;
        for(int iy = from_iy; iy <= to_iy; ++iy){
            int j = 2 * (iy - from_iy) + 1;
            for(int iz = from_iz; iz <= to_iz; ++iz){
                int k = 2 * (iz - from_iz) + 1;
                //! まず v^h -> v^2h にコピーする
                //! 7-point stencil restriction
                parentPhi[ix][iy][iz] = div * (6.0 * phi[i][j][k] + phi[i - 1][j][k] + phi[i + 1][j][k] + phi[i][j - 1][k] + phi[i][j + 1][k] + phi[i][j][k - 1] + phi[i][j][k + 1]);
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

void ChildGrid::insertAMRBlockInfo(SimpleVTK& vtk_gen, const std::string& data_type_name, const std::string& i_timestamp) const {
    vtk_gen.beginBlock(level);
    vtk_gen.beginDataSet(id);
    vtk_gen.setAMRBoxNodeFromParentIndex(parent->getID(), level - 1, from_ix - 1, to_ix - 1, from_iy - 1, to_iy - 1, from_iz - 1, to_iz - 1);
    vtk_gen.setFile("raw_data/" + data_type_name + "_id_" + std::to_string(id) + "_" + i_timestamp + ".vti");
    vtk_gen.endDataSet();
    vtk_gen.endBlock();

    for(const auto& child : children) {
        child->insertAMRBlockInfo(vtk_gen, data_type_name, i_timestamp);
    }
}
