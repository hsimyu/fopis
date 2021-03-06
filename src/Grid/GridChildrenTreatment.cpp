#include "grid.hpp"

void Grid::makeChild(const int _from_ix, const int _from_iy, const int _from_iz, const int _to_ix, const int _to_iy, const int _to_iz) {
    std::unique_ptr<ChildGrid> new_child_ptr(new ChildGrid(this, _from_ix, _from_iy, _from_iz, _to_ix, _to_iy, _to_iz));

    this->addChild( std::move(new_child_ptr) );

    incrementSumOfChild();
}

void Grid::addChild(std::unique_ptr<ChildGrid>&& child) {
    children.push_back(std::move(child));

    //! 粒子を内部に移動
    const int new_index = static_cast<int>(children.size() - 1);
    this->moveParticlesIntoSpecifiedChild(new_index);

    //! Childrenの存在場所を塗り潰す
    this->mapWithNewChild(new_index);
}

//! ChildMap の resize と 初期化
void Grid::initializeChildMap(void) {
    ChildDefinedMapInt::extent_gen mapExtentGen;
    child_map.resize(mapExtentGen[nx + 2][ny + 2][nz + 2]);
}

void Grid::mapWithNewChild(const int child_index) {
    auto& child = children[child_index];

    const auto child_f_ix = child->getFromIX();
    const auto child_t_ix = child->getToIX();
    const auto child_f_iy = child->getFromIY();
    const auto child_t_iy = child->getToIY();
    const auto child_f_iz = child->getFromIZ();
    const auto child_t_iz = child->getToIZ();

    for(int i = child_f_ix; i <= child_t_ix; ++i) {
        for(int j = child_f_iy; j <= child_t_iy; ++j) {
            for(int k = child_f_iz; k <= child_t_iz; ++k) {
                if (i == child_f_ix || i == child_t_ix || j == child_f_iy || j == child_t_iy || k == child_f_iz || k == child_t_iz) {
                    child_map[i][j][k].tag = CHILD_MAP_TAG::EDGE;
                    child_map[i][j][k].child_indices.push_back(child_index);
                } else {
                    child_map[i][j][k].tag = CHILD_MAP_TAG::INNER;
                    child_map[i][j][k].child_indices.push_back(child_index);
                }
            }
        }
    }
}

//! Childを削除する前に存在していた場所を塗り戻す
void Grid::resetChildMapWithSpecifiedChild(const int child_index) {
    auto& child = children[child_index];

    for(int i = child->getFromIX(); i <= child->getToIX(); ++i) {
        for(int j = child->getFromIY(); j <= child->getToIY(); ++j) {
            for(int k = child->getFromIZ(); k <= child->getToIZ(); ++k) {
                child_map[i][j][k].tag = CHILD_MAP_TAG::NOT_EXIST;
                child_map[i][j][k].child_indices.clear();
            }
        }
    }
}

bool Grid::checkSpecifiedChildDoesCoverThisPosition(const int index, const int i, const int j, const int k) const {
    auto& child = children[index];
    if (
        i >= child->getFromIX() && i < child->getToIX() &&
        j >= child->getFromIY() && j < child->getToIY() &&
        k >= child->getFromIZ() && k < child->getToIZ()
    ) return true;
    return false;
}

bool Grid::checkSpecifiedChildDoesCoverThisPosition(const int index, const Position& pos) const {
    return checkSpecifiedChildDoesCoverThisPosition(index, pos.i, pos.j, pos.k);
}

int Grid::getChildIndexIfCovered(const int i, const int j, const int k) const {
    for(int di = 0; di < 2; ++di) {
        for(int dj = 0; dj < 2; ++dj) {
            for(int dk = 0; dk < 2; ++dk) {
                if (child_map[i + di][j + dj][k + dk].tag == CHILD_MAP_TAG::INNER) {
                    return child_map[i + di][j + dj][k + dk].child_indices[0];
                }
            }
        }
    }
    return -1;
}

int Grid::getChildIndexIfCovered(const Position& pos) const {
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

// 子グリッドの境界へPhiをコピーする
void Grid::restrictPhiToChildrenBoundary() {
    tdArray& parentPhi = field->getPhi();

    for(int chidx = 0; chidx < children.size(); ++chidx) {
        auto& child = children[chidx];

        tdArray& childPhi = child->getPhi();
        int child_from_ix = child->getFromIX();
        int child_from_iy = child->getFromIY();
        int child_from_iz = child->getFromIZ();

        const int child_nx = child->getNX();
        const int child_ny = child->getNY();
        const int child_nz = child->getNZ();

        //! x面束縛
        for(int i = 1; i < child_nx + 1; i += child_nx - 1) {
            int ix = ((i - 1) / 2) + child_from_ix;
            for(int j = 1; j < child_ny + 1; ++j) {
                int iy = (j / 2) + child_from_iy;
                for(int k = 1; k < child_nz + 1; ++k) {
                    int iz = (k / 2) + child_from_iz;
                    const bool j_is_odd = (j % 2 == 1);
                    const bool k_is_odd = (k % 2 == 1);

                    if (j_is_odd && k_is_odd) {
                        childPhi[i][j][k] = parentPhi[ix][iy][iz];
                    } else if (j_is_odd) {
                        childPhi[i][j][k] = 0.5 * (parentPhi[ix][iy][iz] + parentPhi[ix][iy][iz - 1]);
                    } else if (k_is_odd) {
                        childPhi[i][j][k] = 0.5 * (parentPhi[ix][iy][iz] + parentPhi[ix][iy - 1][iz]);
                    } else {
                        childPhi[i][j][k] = 0.25 * (parentPhi[ix][iy][iz] + parentPhi[ix][iy - 1][iz] + parentPhi[ix][iy][iz - 1] + parentPhi[ix][iy - 1][iz - 1]);
                    }
                }
            }
        }
        //! y面束縛
        for(int j = 1; j < child_ny + 1; j += child_ny - 1) {
            int iy = ((j - 1) / 2) + child_from_iy;
            for(int i = 1; i < child_nx + 1; ++i) {
                int ix = (i / 2) + child_from_ix;
                for(int k = 1; k < child_nz + 1; ++k) {
                    int iz = (k / 2) + child_from_iz;
                    const bool i_is_odd = (i % 2 == 1);
                    const bool k_is_odd = (k % 2 == 1);

                    if (i_is_odd && k_is_odd) {
                        childPhi[i][j][k] = parentPhi[ix][iy][iz];
                    } else if (i_is_odd) {
                        childPhi[i][j][k] = 0.5 * (parentPhi[ix][iy][iz] + parentPhi[ix][iy][iz - 1]);
                    } else if (k_is_odd) {
                        childPhi[i][j][k] = 0.5 * (parentPhi[ix][iy][iz] + parentPhi[ix - 1][iy][iz]);
                    } else {
                        childPhi[i][j][k] = 0.25 * (parentPhi[ix][iy][iz] + parentPhi[ix - 1][iy][iz] + parentPhi[ix][iy][iz - 1] + parentPhi[ix - 1][iy][iz - 1]);
                    }
                }
            }
        }
        //! z面束縛
        for(int k = 1; k < child_nz + 1; k += child_nz - 1) {
            int iz = ((k - 1) / 2) + child_from_iz;
            for(int i = 1; i < child_nx + 1; ++i) {
                int ix = (i / 2) + child_from_ix;
                for(int j = 1; j < child_ny + 1; ++j) {
                    int iy = (j / 2) + child_from_iy;
                    const bool i_is_odd = (i % 2 == 1);
                    const bool j_is_odd = (j % 2 == 1);

                    if (i_is_odd && j_is_odd) {
                        childPhi[i][j][k] = parentPhi[ix][iy][iz];
                    } else if (i_is_odd) {
                        childPhi[i][j][k] = 0.5 * (parentPhi[ix][iy][iz] + parentPhi[ix][iy - 1][iz]);
                    } else if (j_is_odd) {
                        childPhi[i][j][k] = 0.5 * (parentPhi[ix][iy][iz] + parentPhi[ix - 1][iy][iz]);
                    } else {
                        childPhi[i][j][k] = 0.25 * (parentPhi[ix][iy][iz] + parentPhi[ix][iy - 1][iz] + parentPhi[ix - 1][iy][iz] + parentPhi[ix - 1][iy - 1][iz]);
                    }
                }
            }
        }
    }
}

void Grid::makeGrandChildRecursive() {
    this->makeChild(3, 3, 3, 8, 8, 8);

    if (level < 4) {
        for(auto& child : children) {
            child->makeGrandChildRecursive();
        }
    }
}