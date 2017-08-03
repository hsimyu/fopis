#include "spacecraft.hpp"

//! static 変数の実体
unsigned int Spacecraft::num_of_spacecraft = 0;

void Spacecraft::construct(const size_t nx, const size_t ny, const size_t nz) {
    ++num_of_spacecraft;

    // Node ベース, glue cell は要らない?
    objectArray::extent_gen objectExtents;
    object_map.resize(objectExtents[nx + 2][ny + 2][nz + 2]);

    // 初期化
    for(size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t k = 0; k < nz; ++k) {
                object_map[i][j][k] = false;
            }
        }
    }

    // テスト定義
    for(size_t i = nx/4; i < 3 * nx/4; ++i) {
        for (size_t j = ny / 4; j < 3 * ny / 4; ++j) {
            for (size_t k = nz / 4; k < 3 * nz / 4; ++k) {
                object_map[i][j][k] = true;
            }
        }
    }

    //! 電荷配列初期化
    for(size_t i = 0; i < nx + 2; ++i) {
        for (size_t j = 0; j < ny + 2; ++j) {
            for (size_t k = 0; k < nz + 2; ++k) {
                charge_map[i][j][k] = 0.0;
            }
        }
    }
}

bool Spacecraft::isIncluded(const Position& pos) {
    return  object_map[pos.raw_i    ][pos.raw_j    ][pos.raw_k    ] &&
            object_map[pos.raw_i + 1][pos.raw_j    ][pos.raw_k    ] &&
            object_map[pos.raw_i    ][pos.raw_j + 1][pos.raw_k    ] &&
            object_map[pos.raw_i + 1][pos.raw_j + 1][pos.raw_k    ] &&
            object_map[pos.raw_i    ][pos.raw_j    ][pos.raw_k + 1] &&
            object_map[pos.raw_i + 1][pos.raw_j    ][pos.raw_k + 1] &&
            object_map[pos.raw_i    ][pos.raw_j + 1][pos.raw_k + 1] &&
            object_map[pos.raw_i + 1][pos.raw_j + 1][pos.raw_k + 1];
}

std::ostream& operator<<(std::ostream& ost, const Spacecraft& spc) {
    ost << "name:" << spc.getName() << endl;

    return ost;
}
