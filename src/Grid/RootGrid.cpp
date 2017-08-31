#include "grid.hpp"

void RootGrid::solvePoisson(void) {
    const int DEFAULT_ITERATION_LOOP = 500;
    field->solvePoisson(DEFAULT_ITERATION_LOOP, dx);
}

void RootGrid::updateEfield(void) {
    field->updateEfield(dx);
}

void RootGrid::updateEfieldFDTD(void) {
    field->updateEfieldFDTD(dx, dt);
}

void RootGrid::updateBfield(void) {
    field->updateBfield(dx, nx, ny, nz, dt);
}
//! 場の resize を行う
void Grid::initializeField(void){
    tdArray::extent_gen tdExtents;

    const int cx = nx + 2;
    const int cy = ny + 2;
    const int cz = nz + 2;

    field->getPhi().resize(tdExtents[cx][cy][cz]);

    auto& rho = field->getRho();
    for(int i = 0; i < Environment::num_of_particle_types + 1; ++i) {
        //! 総和のtdArray + 粒子種毎のtdArrayを直接配置で生成する
        rho.emplace_back(tdExtents[cx][cy][cz], boost::fortran_storage_order());
    }

    field->getEx().resize(tdExtents[cx-1][cy][cz]);
    field->getEy().resize(tdExtents[cx][cy-1][cz]);
    field->getEz().resize(tdExtents[cx][cy][cz-1]);

    field->getBx().resize(tdExtents[cx][cy-1][cz-1]);
    field->getBy().resize(tdExtents[cx-1][cy][cz-1]);
    field->getBz().resize(tdExtents[cx-1][cy-1][cz]);

    // reference fields have the same size as nodal size
    field->getExRef().resize(tdExtents[cx][cy][cz]);
    field->getEyRef().resize(tdExtents[cx][cy][cz]);
    field->getEzRef().resize(tdExtents[cx][cy][cz]);
    field->getBxRef().resize(tdExtents[cx][cy][cz]);
    field->getByRef().resize(tdExtents[cx][cy][cz]);
    field->getBzRef().resize(tdExtents[cx][cy][cz]);

    //! 電流密度は Edge 要素なので Efield と同じ要素数を持つ
    field->getJx().resize(tdExtents[cx-1][cy][cz]);
    field->getJy().resize(tdExtents[cx][cy-1][cz]);
    field->getJz().resize(tdExtents[cx][cy][cz-1]);
}