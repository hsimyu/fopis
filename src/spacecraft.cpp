#include "spacecraft.hpp"

//! static 変数の実体
unsigned int Spacecraft::num_of_spacecraft = 0;

void Spacecraft::construct(const size_t nx, const size_t ny, const size_t nz) {
    ++num_of_spacecraft;

    objectArray::extent_gen extents;
    // Node ベース, glue cell は要らない?
    object_map.resize(extents[nx][ny][nz]);

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
}

std::ostream& operator<<(std::ostream& ost, const Spacecraft& spc) {
    ost << "name:" << spc.getName() << endl;

    return ost;
}
