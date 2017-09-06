#include "grid.hpp"

void Grid::makeChild(const int _from_ix, const int _from_iy, const int _from_iz, const int _to_ix, const int _to_iy, const int _to_iz) {
    this->addChild(
        std::make_unique<ChildGrid>(shared_from_this(), _from_ix, _from_iy, _from_iz, _to_ix, _to_iy, _to_iz)
    );
    incrementSumOfChild();
}

void Grid::addChild(std::unique_ptr<ChildGrid>&& child) {
    children.push_back( std::move(child) );
}

int Grid::getChildIndexIfCovered(const int i, const int j, const int k) const {
    for(int index = 0; index < this->getChildrenLength(); ++index) {
        auto& child = children[index];
        if (
            i >= child->getFromIX() && i < child->getToIX() &&
            j >= child->getFromIY() && j < child->getToIY() &&
            k >= child->getFromIZ() && k < child->getToIZ()
        ) return index;
    }
    return -1;
}

int Grid::getChildIndexIfCovered(const Position& pos) const {
    return this->getChildIndexIfCovered(pos.i, pos.j, pos.k);
}

int Grid::getChildIndexIfCovered(Position&& pos) const {
    return this->getChildIndexIfCovered(pos.i, pos.j, pos.k);
}

void Grid::correctChildrenPhi() {
    tdArray& poisson_error = field->getPoissonError();

    for(int chidx = 0; chidx < children.size(); ++chidx) {
        tdArray& childValue = children[chidx]->getPhi();

        int child_from_ix = children[chidx]->getFromIX();
        int child_from_iy = children[chidx]->getFromIY();
        int child_from_iz = children[chidx]->getFromIZ();
        int child_to_ix = children[chidx]->getToIX();
        int child_to_iy = children[chidx]->getToIY();
        int child_to_iz = children[chidx]->getToIZ();

        for(int ix = child_from_ix; ix <= child_to_ix; ++ix){
            int i = 2 * (ix - child_from_ix) + 1;
            for(int iy = child_from_iy; iy <= child_to_iy; ++iy){
                int j = 2 * (iy - child_from_iy) + 1;
                for(int iz = child_from_iz; iz <= child_to_iz; ++iz){
                    int k = 2 * (iz - child_from_iz) + 1;

                    childValue[i][j][k] += poisson_error[ix][iy][iz];

                    if(iz != child_to_iz) {
                        childValue[i][j][k + 1] += 0.5 * (poisson_error[ix][iy][iz] + poisson_error[ix][iy][iz + 1]);
                    }

                    if(iy != child_to_iy) {
                        childValue[i][j + 1][k] += 0.5 * (poisson_error[ix][iy][iz] + poisson_error[ix][iy + 1][iz]);

                        if(iz != child_to_iz) {
                            childValue[i][j + 1][k + 1] += 0.25 * (poisson_error[ix][iy][iz] + poisson_error[ix][iy][iz + 1] + poisson_error[ix][iy + 1][iz] + poisson_error[ix][iy + 1][iz + 1]);
                        }
                    }

                    if(ix != child_to_ix) {
                        childValue[i + 1][j][k] += 0.5 * (poisson_error[ix][iy][iz] + poisson_error[ix + 1][iy][iz]);

                        if(iz != child_to_iz) {
                            childValue[i + 1][j][k + 1] += 0.25 * (poisson_error[ix][iy][iz] + poisson_error[ix][iy][iz + 1] + poisson_error[ix + 1][iy][iz] + poisson_error[ix + 1][iy][iz + 1]);
                        }

                        if(iy != child_to_iy) {
                            childValue[i + 1][j + 1][k] += 0.25 * (poisson_error[ix][iy][iz] + poisson_error[ix + 1][iy][iz] + poisson_error[ix][iy + 1][iz] + poisson_error[ix + 1][iy + 1][iz]);

                            if(iz != child_to_iz) {
                                childValue[i + 1][j + 1][k + 1] += 0.125 *
                                    ( poisson_error[ix][iy][iz] + poisson_error[ix][iy][iz + 1] + poisson_error[ix][iy + 1][iz] + poisson_error[ix][iy + 1][iz + 1]
                                    + poisson_error[ix + 1][iy][iz] + poisson_error[ix + 1][iy][iz + 1] + poisson_error[ix + 1][iy + 1][iz] + poisson_error[ix + 1][iy + 1][iz + 1]);
                            }
                        }
                    }
                }
            }
        }
    }
}

void Grid::updateChildrenPhi() {
    tdArray& parent_phi = field->getPhi();

    for(int chidx = 0; chidx < children.size(); ++chidx) {
        tdArray& childValue = children[chidx]->getPhi();

        int child_from_ix = children[chidx]->getFromIX();
        int child_from_iy = children[chidx]->getFromIY();
        int child_from_iz = children[chidx]->getFromIZ();
        int child_to_ix = children[chidx]->getToIX();
        int child_to_iy = children[chidx]->getToIY();
        int child_to_iz = children[chidx]->getToIZ();

        for(int ix = child_from_ix; ix <= child_to_ix; ++ix){
            int i = 2 * (ix - child_from_ix) + 1;
            for(int iy = child_from_iy; iy <= child_to_iy; ++iy){
                int j = 2 * (iy - child_from_iy) + 1;
                for(int iz = child_from_iz; iz <= child_to_iz; ++iz){
                    int k = 2 * (iz - child_from_iz) + 1;

                    childValue[i][j][k] = parent_phi[ix][iy][iz];

                    if(iz != child_to_iz) {
                        childValue[i][j][k + 1] = 0.5 * (parent_phi[ix][iy][iz] + parent_phi[ix][iy][iz + 1]);
                    }

                    if(iy != child_to_iy) {
                        childValue[i][j + 1][k] = 0.5 * (parent_phi[ix][iy][iz] + parent_phi[ix][iy + 1][iz]);

                        if(iz != child_to_iz) {
                            childValue[i][j + 1][k + 1] = 0.25 * (parent_phi[ix][iy][iz] + parent_phi[ix][iy][iz + 1] + parent_phi[ix][iy + 1][iz] + parent_phi[ix][iy + 1][iz + 1]);
                        }
                    }

                    if(ix != child_to_ix) {
                        childValue[i + 1][j][k] = 0.5 * (parent_phi[ix][iy][iz] + parent_phi[ix + 1][iy][iz]);

                        if(iz != child_to_iz) {
                            childValue[i + 1][j][k + 1] = 0.25 * (parent_phi[ix][iy][iz] + parent_phi[ix][iy][iz + 1] + parent_phi[ix + 1][iy][iz] + parent_phi[ix + 1][iy][iz + 1]);
                        }

                        if(iy != child_to_iy) {
                            childValue[i + 1][j + 1][k] = 0.25 * (parent_phi[ix][iy][iz] + parent_phi[ix + 1][iy][iz] + parent_phi[ix][iy + 1][iz] + parent_phi[ix + 1][iy + 1][iz]);

                            if(iz != child_to_iz) {
                                childValue[i + 1][j + 1][k + 1] = 0.125 *
                                    ( parent_phi[ix][iy][iz] + parent_phi[ix][iy][iz + 1] + parent_phi[ix][iy + 1][iz] + parent_phi[ix][iy + 1][iz + 1]
                                    + parent_phi[ix + 1][iy][iz] + parent_phi[ix + 1][iy][iz + 1] + parent_phi[ix + 1][iy + 1][iz] + parent_phi[ix + 1][iy + 1][iz + 1]);
                            }
                        }
                    }
                }
            }
        }
    }
}

// 子グリッドへRhoをコピーする
void Grid::interpolateRhoValueToChildren() {
    RhoArray& tdValue = field->getRho();

    for(int chidx = 0; chidx < children.size(); ++chidx) {
        RhoArray& childValue = children[chidx]->getRho();

        int child_from_ix = children[chidx]->getFromIX();
        int child_from_iy = children[chidx]->getFromIY();
        int child_from_iz = children[chidx]->getFromIZ();
        int child_to_ix = children[chidx]->getToIX();
        int child_to_iy = children[chidx]->getToIY();
        int child_to_iz = children[chidx]->getToIZ();

        for(int ix = child_from_ix; ix <= child_to_ix; ++ix) {
            int i = 2 * (ix - child_from_ix) + 1;
            for(int iy = child_from_iy; iy <= child_to_iy; ++iy) {
                int j = 2 * (iy - child_from_iy) + 1;
                for(int iz = child_from_iz; iz <= child_to_iz; ++iz) {
                    int k = 2 * (iz - child_from_iz) + 1;

                    childValue[0][i][j][k] = tdValue[0][ix][iy][iz];

                    if(iz != child_to_iz) {
                        childValue[0][i][j][k + 1] = 0.5 * (tdValue[0][ix][iy][iz] + tdValue[0][ix][iy][iz + 1]);
                    }

                    if(iy != child_to_iy) {
                        childValue[0][i][j + 1][k] = 0.5 * (tdValue[0][ix][iy][iz] + tdValue[0][ix][iy + 1][iz]);

                        if(iz != child_to_iz) {
                            childValue[0][i][j + 1][k + 1] = 0.25 * (tdValue[0][ix][iy][iz] + tdValue[0][ix][iy][iz + 1] + tdValue[0][ix][iy + 1][iz] + tdValue[0][ix][iy + 1][iz + 1]);
                        }
                    }

                    if(ix != child_to_ix) {
                        childValue[0][i + 1][j][k] = 0.5 * (tdValue[0][ix][iy][iz] + tdValue[0][ix + 1][iy][iz]);

                        if(iz != child_to_iz) {
                            childValue[0][i + 1][j][k + 1] = 0.25 * (tdValue[0][ix][iy][iz] + tdValue[0][ix][iy][iz + 1] + tdValue[0][ix + 1][iy][iz] + tdValue[0][ix + 1][iy][iz + 1]);
                        }

                        if(iy != child_to_iy) {
                            childValue[0][i + 1][j + 1][k] = 0.25 * (tdValue[0][ix][iy][iz] + tdValue[0][ix + 1][iy][iz] + tdValue[0][ix][iy + 1][iz] + tdValue[0][ix + 1][iy + 1][iz]);

                            if(iz != child_to_iz) {
                                childValue[0][i + 1][j + 1][k + 1] = 0.125 *
                                    ( tdValue[0][ix][iy][iz] + tdValue[0][ix][iy][iz + 1] + tdValue[0][ix][iy + 1][iz] + tdValue[0][ix][iy + 1][iz + 1]
                                    + tdValue[0][ix + 1][iy][iz] + tdValue[0][ix + 1][iy][iz + 1] + tdValue[0][ix + 1][iy + 1][iz] + tdValue[0][ix + 1][iy + 1][iz + 1]);
                            }
                        }
                    }
                }
            }
        }

        if(children[chidx]->getChildrenLength() > 0) children[chidx]->interpolateRhoValueToChildren();
    }
}

// 子グリッドへPhiをコピーする
void Grid::restrictPhiValueToChildren() {
    tdArray& parentPhi = field->getPhi();

    for(int chidx = 0; chidx < children.size(); ++chidx) {
        auto& child = children[chidx];

        tdArray& childPhi = child->getPhi();
        int child_from_ix = child->getFromIX();
        int child_from_iy = child->getFromIY();
        int child_from_iz = child->getFromIZ();
        int child_to_ix = child->getToIX();
        int child_to_iy = child->getToIY();
        int child_to_iz = child->getToIZ();

        const int child_nx = child->getNX();
        const int child_ny = child->getNY();
        const int child_nz = child->getNZ();

        //! x面束縛
        //! child_nx は 必ず奇数なので-1すれば2で割れる
        for(int i = 1; i < child_nx + 1; i += child_nx - 1) {
            int ix = ((i - 1)/2) + child_from_ix;
            for(int iy = child_from_iy; iy <= child_to_iy; ++iy) {
                int j = 2 * (iy - child_from_iy) + 1;
                for(int iz = child_from_iz; iz <= child_to_iz; ++iz) {
                    int k = 2 * (iz - child_from_iz) + 1;
                    childPhi[i][j][k] = parentPhi[ix][iy][iz];

                    if(iz != child_to_iz) {
                        childPhi[i][j][k + 1] = 0.5 * (parentPhi[ix][iy][iz] + parentPhi[ix][iy][iz + 1]);
                    }

                    if(iy != child_to_iy) {
                        childPhi[i][j + 1][k] = 0.5 * (parentPhi[ix][iy][iz] + parentPhi[ix][iy + 1][iz]);

                        if(iz != child_to_iz) {
                            childPhi[i][j + 1][k + 1] = 0.25 * (parentPhi[ix][iy][iz] + parentPhi[ix][iy][iz + 1] + parentPhi[ix][iy + 1][iz] + parentPhi[ix][iy + 1][iz + 1]);
                        }
                    }
                }
            }
        }

        //! y面束縛
        for(int j = 1; j < child_ny + 1; j += child_ny - 1) {
            int iy = ((j - 1)/2) + child_from_iy;
            for(int ix = child_from_ix; ix <= child_to_ix; ++ix) {
                int i = 2 * (ix - child_from_ix) + 1;
                for(int iz = child_from_iz; iz <= child_to_iz; ++iz) {
                    int k = 2 * (iz - child_from_iz) + 1;
                    childPhi[i][j][k] = parentPhi[ix][iy][iz];

                    if(iz != child_to_iz) {
                        childPhi[i][j][k + 1] = 0.5 * (parentPhi[ix][iy][iz] + parentPhi[ix][iy][iz + 1]);
                    }

                    if(ix != child_to_ix) {
                        childPhi[i + 1][j][k] = 0.5 * (parentPhi[ix][iy][iz] + parentPhi[ix + 1][iy][iz]);

                        if(iz != child_to_iz) {
                            childPhi[i + 1][j][k + 1] = 0.25 * (parentPhi[ix][iy][iz] + parentPhi[ix][iy][iz + 1] + parentPhi[ix + 1][iy][iz] + parentPhi[ix + 1][iy][iz + 1]);
                        }
                    }
                }
            }
        }

        //! z面束縛
        for(int k = 1; k < child_nz + 1; k += child_nz - 1) {
            int iz = ((k - 1)/2) + child_from_iz;
            for(int ix = child_from_ix; ix <= child_to_ix; ++ix) {
                int i = 2 * (ix - child_from_ix) + 1;
                for(int iy = child_from_iy; iy <= child_to_iy; ++iy) {
                    int j = 2 * (iy - child_from_iy) + 1;
                    childPhi[i][j][k] = parentPhi[ix][iy][iz];

                    if(iy != child_to_iy) {
                        childPhi[i][j + 1][k] = 0.5 * (parentPhi[ix][iy][iz] + parentPhi[ix][iy + 1][iz]);
                    }

                    if(ix != child_to_ix) {
                        childPhi[i + 1][j][k] = 0.5 * (parentPhi[ix][iy][iz] + parentPhi[ix + 1][iy][iz]);

                        if(iy != child_to_iy) {
                            childPhi[i + 1][j + 1][k] = 0.25 * (parentPhi[ix][iy][iz] + parentPhi[ix][iy + 1][iz] + parentPhi[ix + 1][iy][iz] + parentPhi[ix + 1][iy + 1][iz]);
                        }
                    }
                }
            }
        }
    }
}
