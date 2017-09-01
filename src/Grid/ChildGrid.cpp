#include "grid.hpp"
#include "field.hpp"
#include "utils.hpp"

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
    const int DEFAULT_ITERATION_LOOP = 500;
    field->solvePoisson(DEFAULT_ITERATION_LOOP, dx);
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

void ChildGrid::copyScalarToParent(std::string varname){
    tdArray& tdValue = field->getScalar(varname);

    // @note: OpenMP
    tdArray& parentValue = parent->getScalar(varname);

    for(int ix = from_ix; ix <= to_ix; ++ix){
        int i = 2 * (ix - from_ix) + 1;
        for(int iy = from_iy; iy <= to_iy; ++iy){
            int j = 2 * (iy - from_iy) + 1;
            for(int iz = from_iz; iz <= to_iz; ++iz){
                int k = 2 * (iz - from_iz) + 1;

                // とりあえずダイレクトにコピーする
                parentValue[ix][iy][iz] = tdValue[i][j][k];
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
