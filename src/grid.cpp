#include "global.hpp"
#include "grid.hpp"
#include "field.hpp"
#include "particle.hpp"
#include "environment.hpp"
#include "mpiw.hpp"
#include <silo.h>
#include <random>

// Unique ID の実体
unsigned int Grid::nextID = 0;

void Grid::makeChild(const int _from_ix, const int _from_iy, const int _from_iz, const int _to_ix, const int _to_iy, const int _to_iz) {
    Grid* child = new Grid(this, _from_ix, _from_iy, _from_iz, _to_ix, _to_iy, _to_iz);
    this->addChild(child);
    incrementSumOfChild();
}

void Grid::addChild(Grid* child) {
    children.push_back(child);
}

void Grid::decrementSumOfChild() {
    --sumTotalNumOfChildGrids;

    if(level > 0) {
        parent->decrementSumOfChild();
    }
}

void Grid::incrementSumOfChild() {
    ++sumTotalNumOfChildGrids;

    if(level > 0) {
        parent->incrementSumOfChild();
    }
}

void Grid::initializeField(void){
    Field* field = new Field;
    tdArray::extent_gen tdExtents;

    const int cx = nx + 2;
    const int cy = ny + 2;
    const int cz = nz + 2;

    field->getPhi().resize(tdExtents[cx][cy][cz]);
    field->getRho().resize(tdExtents[cx][cy][cz]);

    // @note: Dirichlet Boundary
    Utils::clearBoundaryValues(field->getPhi(), cx, cy, cz);

    field->getEx().resize(tdExtents[cx-1][cy][cz]);
    field->getEy().resize(tdExtents[cx][cy-1][cz]);
    field->getEz().resize(tdExtents[cx][cy][cz-1]);

    // reference fields have the same size as nodal size
    field->getExRef().resize(tdExtents[cx][cy][cz]);
    field->getEyRef().resize(tdExtents[cx][cy][cz]);
    field->getEzRef().resize(tdExtents[cx][cy][cz]);

    field->getBx().resize(tdExtents[cx][cy-1][cz-1]);
    field->getBy().resize(tdExtents[cx-1][cy][cz-1]);
    field->getBz().resize(tdExtents[cx-1][cy-1][cz]);

    this->setField(field);
}

// root grid constructor
Grid::Grid(void){
    //! - コンストラクタにEnvironmentクラスが渡された場合、
    //! レベル0のGridを作成します.
    level = 0;
    sumTotalNumOfChildGrids = 0;

    //! UniqueなIDをセット
    id = this->getNextID();

    nx = Environment::cell_x;
    ny = Environment::cell_y;
    nz = Environment::cell_z;
    dx = Utils::Normalizer::normalizeLength(Environment::dx);

    //! @{
    //! Root Gridの場合の親グリッドは、計算空間を全て統合した空間として、
    //! その上にプロセス分割されたグリッドが乗っていると考える
    from_ix = MPIw::Environment::xrank * Environment::cell_x;
    from_iy = MPIw::Environment::yrank * Environment::cell_y;
    from_iz = MPIw::Environment::zrank * Environment::cell_z;
    to_ix = from_ix + nx;
    to_iy = from_iy + ny;
    to_iz = from_iz + nz;
    //! @note: base_x, base_y, base_zは正規化された長さ
    base_x = dx * static_cast<double>(from_ix);
    base_y = dx * static_cast<double>(from_iy);
    base_z = dx * static_cast<double>(from_iz);
    //! @}

    // Field初期化
    this->initializeField();

    //! 粒子位置の上限を設定
    //! [0, max_x)になるよう1e-20を引いておく
    const double max_x = static_cast<double>(Environment::cell_x) - 1e-20;
    const double max_y = static_cast<double>(Environment::cell_y) - 1e-20;
    const double max_z = static_cast<double>(Environment::cell_z) - 1e-20;

    // std::random_device rnd;
    const int random_src_x = 10684930 + MPIw::Environment::rank;
    const int random_src_y = 99881 + MPIw::Environment::rank;
    const int random_src_z = 861200045 + MPIw::Environment::rank;
    const int random_src_vx = 930 + MPIw::Environment::rank;
    const int random_src_vy = 98076621 + MPIw::Environment::rank;
    const int random_src_vz = 7662566 + MPIw::Environment::rank;
    std::mt19937 mt_x(random_src_x);
    std::mt19937 mt_y(random_src_y);
    std::mt19937 mt_z(random_src_z);
    std::mt19937 mt_vx(random_src_vx);
    std::mt19937 mt_vy(random_src_vy);
    std::mt19937 mt_vz(random_src_vz);

    std::uniform_real_distribution<> dist_x(0.0, max_x);
    std::uniform_real_distribution<> dist_y(0.0, max_y);
    std::uniform_real_distribution<> dist_z(0.0, max_z);

    // particlesは空のstd::vector< std::vector<Particle> >として宣言されている
    // particle types 分だけresize
    particles.resize(Environment::num_of_particle_types);

    for(int id = 0; id < Environment::num_of_particle_types; ++id){
        int pnum = Environment::ptype[id].getTotalNumber();

        //! 各粒子分のメモリをreserve
        particles[id].reserve(Environment::max_particle_num);

        //! particle_number分のコンストラクタが呼ばれる
        particles[id].resize(pnum);

        const double deviation = Environment::ptype[id].calcDeviation();
        std::normal_distribution<> dist_vx(0.0, deviation);
        std::normal_distribution<> dist_vy(0.0, deviation);
        std::normal_distribution<> dist_vz(0.0, deviation);

        //! - 粒子はレベル0グリッドにのみ所属します
        for(int i = 0; i < pnum; ++i){
            particles[id][i].setPosition(dist_x(mt_x), dist_y(mt_y), dist_z(mt_z));
            particles[id][i].setVelocity(dist_vx(mt_vx), dist_vy(mt_vy), dist_vz(mt_vz));
            particles[id][i].typeId = id;
        }
    }
}

//! child grid constructor
//! GridコンストラクタにGridが渡された場合、
//! そのGridを親とした子グリッドを生成します
Grid::Grid(Grid* g,
        const int _from_ix, const int _from_iy, const int _from_iz,
        const int _to_ix,   const int _to_iy,   const int _to_iz)
{
    const double refineRatio = 2.0;

    //! UniqueなIDをセット
    id = this->getNextID();

    parent = g;
    level = g->getLevel() + 1;
    sumTotalNumOfChildGrids = 0;

    //! @{
    //! 子グリッドの場合, from_ix, to_ix変数は純粋に親グリッドの何番目に乗っているかを表す
    //! Glue cell 分も考慮した方に乗る
    from_ix = _from_ix;
    from_iy = _from_iy;
    from_iz = _from_iz;
    to_ix = _to_ix;
    to_iy = _to_iy;
    to_iz = _to_iz;
    base_x = g->getBaseX() + g->getDX() * (from_ix - 1);
    base_y = g->getBaseY() + g->getDX() * (from_iy - 1);
    base_z = g->getBaseZ() + g->getDX() * (from_iz - 1);
    //! @}

    // patchの大きさを計算
    nx = (_to_ix - _from_ix) * 2 + 1;
    ny = (_to_iy - _from_iy) * 2 + 1;
    nz = (_to_iz - _from_iz) * 2 + 1;

    // refineRatioは2で固定
    dx = g->getDX() / refineRatio;

    checkGridValidness();

    // Field初期化
    this->initializeField();
}

void Grid::checkGridValidness() {
    const double refineRatio = 2.0;

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

    if(!isValid) MPIw::Environment::exitWithFinalize(1);
}

//! 粒子の位置から電荷を空間電荷にする
//! 基本的にはroot_gridに対してのみ呼ぶ
void Grid::updateRho() {
    tdArray& rho = field->getRho();

    for(int pid = 0; pid < Environment::num_of_particle_types; ++pid){
        int max_pnum = Environment::ptype[pid].getTotalNumber();
        double q = Environment::ptype[pid].getCharge();

        for(int pnum = 0; pnum < max_pnum; ++pnum){
            Particle& p = particles[pid][pnum];
            if(p.isValid) {
                Position pos(p);
                int i = pos.i, j = pos.j, k = pos.k;

#ifdef DEBUG
                if( (i >= nx + 2) || (j >= ny + 2) || (k >= nz + 2)) {
                    cout << Environment::rankStr() << format("[Particle Position Error]: xyz: %5f %5f %5f") % pos.x % pos.y % pos.z << endl;
                    cout << Environment::rankStr() << format("[Particle Position Error]: ijk: %d %d %d") % i % j % k << endl;
                }
#endif

                rho[i  ][j  ][k] += pos.dx2 * pos.dy2 * pos.dz2 * q;
                rho[i+1][j  ][k] += pos.dx1 * pos.dy2 * pos.dz2 * q;
                rho[i  ][j+1][k] += pos.dx2 * pos.dy1 * pos.dz2 * q;
                rho[i+1][j+1][k] += pos.dx1 * pos.dy1 * pos.dz2 * q;

                rho[i  ][j  ][k+1] += pos.dx2 * pos.dy2 * pos.dz1 * q;
                rho[i+1][j  ][k+1] += pos.dx1 * pos.dy2 * pos.dz1 * q;
                rho[i  ][j+1][k+1] += pos.dx2 * pos.dy1 * pos.dz1 * q;
                rho[i+1][j+1][k+1] += pos.dx1 * pos.dy1 * pos.dz1 * q;
            }
        }
    }
}

void Grid::updateParticleVelocity(void) {
    tdArray& exref = field->getExRef();
    tdArray& eyref = field->getEyRef();
    tdArray& ezref = field->getEzRef();

    for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
        double qm = 0.5 * (Environment::ptype[pid].getCharge()) * (Environment::ptype[pid].getMass());

        for(int pnum = 0; pnum < particles[pid].size(); ++pnum){
            Particle& p = particles[pid][pnum];
            if(p.isValid) {
                //! 毎回生成するよりコピーのが早い？
                Position pos(p);

                int i = pos.i, j = pos.j, k = pos.k;
                double v1 = qm * pos.dx2 * pos.dy2 * pos.dz2;
                double v2 = qm * pos.dx1 * pos.dy2 * pos.dz2;
                double v3 = qm * pos.dx2 * pos.dy1 * pos.dz2;
                double v4 = qm * pos.dx1 * pos.dy1 * pos.dz2;
                double v5 = qm * pos.dx2 * pos.dy2 * pos.dz1;
                double v6 = qm * pos.dx1 * pos.dy2 * pos.dz1;
                double v7 = qm * pos.dx2 * pos.dy1 * pos.dz1;
                double v8 = qm * pos.dx1 * pos.dy1 * pos.dz1;

                double ex =  v1*exref[i][j][k]
                    + v2*exref[i+1][j][k]
                    + v3*exref[i][j+1][k]
                    + v4*exref[i+1][j+1][k]
                    + v5*exref[i][j][k+1]
                    + v6*exref[i+1][j][k+1]
                    + v7*exref[i][j+1][k+1]
                    + v8*exref[i+1][j+1][k+1];
                double ey =  v1*eyref[i][j][k]
                    + v2*eyref[i+1][j][k]
                    + v3*eyref[i][j+1][k]
                    + v4*eyref[i+1][j+1][k]
                    + v5*eyref[i][j][k+1]
                    + v6*eyref[i+1][j][k+1]
                    + v7*eyref[i][j+1][k+1]
                    + v8*eyref[i+1][j+1][k+1];
                double ez =  v1*ezref[i][j][k]
                    + v2*ezref[i+1][j][k]
                    + v3*ezref[i][j+1][k]
                    + v4*ezref[i+1][j+1][k]
                    + v5*ezref[i][j][k+1]
                    + v6*ezref[i+1][j][k+1]
                    + v7*ezref[i][j+1][k+1]
                    + v8*ezref[i+1][j+1][k+1];
                double bx = 0.0;
                double by = 0.0;
                double bz = 0.0;
                double boris = 2.0/(1.0 + (bx*bx+by*by+bz*bz));

                double vx1 = p.vx + ex;
                double vy1 = p.vy + ey;
                double vz1 = p.vz + ez;

                double vxt = vx1 + vy1*bz - vz1 * by;
                double vyt = vy1 + vz1*bx - vx1 * bz;
                double vzt = vz1 + vx1*by - vy1 * bx;

                p.vx = vx1 + ex + boris*(vyt*bz - vzt*by);
                p.vy = vy1 + ey + boris*(vzt*bx - vxt*bz);
                p.vz = vz1 + ez + boris*(vxt*by - vyt*bx);
            }
        }
    }
}

void Grid::updateParticlePosition(void) {
    std::vector< std::vector<Particle> > pbuff(6);
    std::vector< std::vector<Particle> > pbuffRecv(6);

    const double slx = dx * nx;
    const double sly = dx * ny;
    const double slz = dx * nz;

    for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
        for(int i = 0; i < particles[pid].size(); ++i){
            Particle& p = particles[pid][i];
            if(p.isValid) {
                p.updatePosition();
                checkXBoundary(pbuff, p, slx);
                checkYBoundary(pbuff, p, sly);
                checkZBoundary(pbuff, p, slz);
            }
        }
    }

    MPIw::Environment::sendRecvParticlesX(pbuff, pbuffRecv);

    for(int axis = 0; axis < 2; ++axis) {
        for(int i = 0; i < pbuffRecv[axis].size(); ++i){
            Particle& p = pbuffRecv[axis][i];

            // 0方向(下側)からやってきた時、領域長分だけ座標をずらす
            if(axis == 0) {
                p.x -= slx;
            } else {
                p.x += slx;
            }

            checkYBoundary(pbuff, p, sly);
            checkZBoundary(pbuff, p, slz);
        }
    }

    MPIw::Environment::sendRecvParticlesY(pbuff, pbuffRecv);

    for(int axis = 2; axis < 4; ++axis) {
        for(int i = 0; i < pbuffRecv[axis].size(); ++i){
            Particle& p = pbuffRecv[axis][i];

            // 0方向(下側)からやってきた時、領域長分だけ座標をずらす
            if(axis == 2) {
                p.x -= sly;
            } else {
                p.x += sly;
            }
            checkZBoundary(pbuff, p, slz);
        }
    }

    MPIw::Environment::sendRecvParticlesZ(pbuff, pbuffRecv);

    for(int axis = 4; axis < 6; ++axis) {
        for(int i = 0; i < pbuffRecv[axis].size(); ++i){
            Particle& p = pbuffRecv[axis][i];

            // 0方向(下側)からやってきた時、領域長分だけ座標をずらす
            if(axis == 4) {
                p.x -= slz;
            } else {
                p.x += slz;
            }
        }
    }

    for(int i = 0; i < 6; ++i) {
        for(auto& p : pbuffRecv[i]){
#ifdef DEBUG
            if( particles[ p.typeId ].capacity() == particles[ p.typeId ].size() ) {
                cout << format("[WARNING] The size of %s array is full.: capacity = %d, size = %d") % Environment::ptype[ p.typeId ].getName() % particles[ p.typeId ].capacity() % particles[ p.typeId ].capacity()<< endl;
            }
#endif
            particles[ p.typeId ].push_back(p);
        }
    }
}

double Grid::getParticleEnergy(void) {
    double res = 0.0;

    for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
        double eachEnergy = 0.0;
        for(int i = 0; i < particles[pid].size(); ++i){
            if(particles[pid][i].isValid) eachEnergy += particles[pid][i].getSquaredMagnitudeOfVelocity();
        }
        res += (0.5 * Environment::ptype[id].getMass() * eachEnergy);
    }

    return res;
}

// 子グリッドへ場の値をコピーする
// この実装はノード to ノードの場合
void Grid::copyScalarToChildren(std::string varname){
    tdArray& tdValue = field->getScalar(varname);

    // @note: OpenMP
    for(int chidx = 0; chidx < children.size(); ++chidx) {
        tdArray& childValue = children[chidx]->getField()->getScalar(varname);

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

                    // cout << format("i, j, k = %d, %d, %d") % i % j % k << endl;
                    // cout << format("ix, iy, iz = %d, %d, %d") % ix % iy % iz << endl;
                    childValue[i][j][k] = tdValue[ix][iy][iz];

                    if(iz != child_to_iz) {
                        childValue[i][j][k + 1] = 0.5 * (tdValue[ix][iy][iz] + tdValue[ix][iy][iz + 1]);
                    }

                    if(iy != child_to_iy) {
                        childValue[i][j + 1][k] = 0.5 * (tdValue[ix][iy][iz] + tdValue[ix][iy + 1][iz]);

                        if(iz != child_to_iz) {
                            childValue[i][j + 1][k + 1] = 0.25 * (tdValue[ix][iy][iz] + tdValue[ix][iy][iz + 1] + tdValue[ix][iy + 1][iz] + tdValue[ix][iy + 1][iz + 1]);
                        }
                    }

                    if(ix != child_to_ix) {
                        childValue[i + 1][j][k] = 0.5 * (tdValue[ix][iy][iz] + tdValue[ix + 1][iy][iz]);

                        if(iz != child_to_iz) {
                            childValue[i + 1][j][k + 1] = 0.25 * (tdValue[ix][iy][iz] + tdValue[ix][iy][iz + 1] + tdValue[ix + 1][iy][iz] + tdValue[ix + 1][iy][iz + 1]);
                        }

                        if(iy != child_to_iy) {
                            childValue[i + 1][j + 1][k] = 0.25 * (tdValue[ix][iy][iz] + tdValue[ix + 1][iy][iz] + tdValue[ix][iy + 1][iz] + tdValue[ix + 1][iy + 1][iz]);

                            if(iz != child_to_iz) {
                                childValue[i + 1][j + 1][k + 1] = 0.125 *
                                    ( tdValue[ix][iy][iz] + tdValue[ix][iy][iz + 1] + tdValue[ix][iy + 1][iz] + tdValue[ix][iy + 1][iz + 1]
                                    + tdValue[ix + 1][iy][iz] + tdValue[ix + 1][iy][iz + 1] + tdValue[ix + 1][iy + 1][iz] + tdValue[ix + 1][iy + 1][iz + 1]);
                            }
                        }
                    }
                }
            }
        }
    }
}

void Grid::copyScalarToParent(std::string varname){
    tdArray& tdValue = field->getScalar(varname);

    // @note: OpenMP
    tdArray& parentValue = parent->getField()->getScalar(varname);

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

// -- AMR utility methods --
int Grid::getMaxLevel() {
    int maxLevel = level; //! 初期値は自分のレベル

    for(int i = 0; i < getChildrenLength(); ++i){
        //! 子がいる時、再帰的に子の最大レベルを取ってくる
        //! @note: パッチ数が増えた時に無駄なチェック時間が大きくなりすぎないか
        int lv = children[i]->getMaxLevel();
        maxLevel = (lv > maxLevel) ? lv : maxLevel;
    }

    return maxLevel;
}

//! 各レベルの子パッチの数を再帰的に返す
//! numOfPatches = {1, 3, 3, 1}
int* Grid::getNumOfPatches() {
    int maxChildLevel = this->getMaxChildLevel();
    int* numOfPatches = new int[maxChildLevel + 1];

    // 初期化
    for(int j = 1; j < maxChildLevel + 1; ++j){
        numOfPatches[j] = 0;
    }

    // 自分を追加
    numOfPatches[0] = 1;

    for(int i = 0; i < getChildrenLength(); ++i){
        int* tempNumOfPatches = children[i]->getNumOfPatches();
        int tempMaxChildLevel = children[i]->getMaxChildLevel() + 1;

        for(int j = 1; j < tempMaxChildLevel + 1; ++j){
            numOfPatches[j] += tempNumOfPatches[j - 1];
        }

        delete [] tempNumOfPatches;
    }

    return numOfPatches;
}

//! 自分以下のGridの子の数を順に格納していく
std::map<int, std::vector<int> > Grid::getChildMapOnRoot(void) {
    // 初期化
    std::map<int, std::vector<int> > childMap;

    this->addChildrenIDToMap(childMap);
    return childMap;
}

void Grid::addChildrenIDToMap(std::map<int, std::vector<int> >& childMap){
    childMap[this->getID()] = this->getChildrenIDs();

    for(int i = 0; i < getChildrenLength(); ++i){
        children[i]->addChildrenIDToMap(childMap);
    }
}

std::vector<int> Grid::getChildrenIDs(void) {
    std::vector<int> childIDs(this->getChildrenLength());

    for(int i = 0; i < childIDs.size(); ++i){
        childIDs[i] = children[i]->getID();
    }

    return childIDs;
}

//! 自分以下のGridのIDを順に格納していく
std::vector< std::vector<int> > Grid::getIDMapOnRoot(void) {
    // 初期化
    std::vector< std::vector<int> > idMap(this->getMaxLevel() + 1);
    this->addIDToVector(idMap);

    return idMap;
}

void Grid::addIDToVector(std::vector< std::vector<int> >& idMap){
    idMap[level].push_back( this->getID() );

    for(int i = 0; i < getChildrenLength(); ++i){
        children[i]->addIDToVector(idMap);
    }
}

// mesh nodesの座標配列を生成
float** Grid::getMeshNodes(int dim) {
    // the array of coordinate arrays
    // @note: メモリリーク防止のため必ずdeleteする
    const float real_dx = Utils::Normalizer::unnormalizeLength(dx);
    const float real_base_x = Utils::Normalizer::unnormalizeLength(base_x);
    const float real_base_y = Utils::Normalizer::unnormalizeLength(base_y);
    const float real_base_z = Utils::Normalizer::unnormalizeLength(base_z);

    float** coordinates = new float*[dim];
    coordinates[0] = new float[nx];
    for(int i = 0; i < nx; ++i) {
	coordinates[0][i] = real_base_x + real_dx * i;
    }
    coordinates[1] = new float[ny];
    for(int i = 0; i < ny; ++i) {
	coordinates[1][i] = real_base_y + real_dx * i;
    }
    coordinates[2] = new float[nz];
    for(int i = 0; i < nz; ++i) {
	coordinates[2][i] = real_base_z + real_dx * i;
    }
    return coordinates;
}


// 渡されたポインタにExtentを入力する
void Grid::addExtent(int* data[6], float* sdata[6], float* rdata[1]){
    if(level == 0) {
        data[0][id] = from_ix;
        data[1][id] = data[0][id] + (nx - 1);
        data[2][id] = from_iy;
        data[3][id] = data[0][id] + (ny - 1);
        data[4][id] = from_iz;
        data[5][id] = data[0][id] + (nz - 1);
    } else {
        data[0][id] = 2 * (from_ix - 1);
        data[1][id] = data[0][id] + (nx - 1);
        data[2][id] = 2 * (from_iy - 1);
        data[3][id] = data[0][id] + (ny - 1);
        data[4][id] = 2 * (from_iz - 1);
        data[5][id] = data[0][id] + (nz - 1);
    }

    sdata[0][id] = base_x;
    sdata[1][id] = base_x + (nx - 1) * dx;
    sdata[2][id] = base_y;
    sdata[3][id] = base_y + (ny - 1) * dx;
    sdata[4][id] = base_z;
    sdata[5][id] = base_z + (nz - 1) * dx;

    rdata[0][id] = dx;

    cout << "[ID " << id << "]" << endl;
    cout << "LogicalExtent imin, imax = " << data[0][id] << " to " << data[1][id] << endl;
    cout << "LogicalExtent jmin, jmax = " << data[2][id] << " to " << data[3][id] << endl;
    cout << "LogicalExtent kmin, kmax = " << data[4][id] << " to " << data[5][id] << endl;

    cout << "SpatialExtent xmin, xmax = " << sdata[0][id] << " to " << sdata[1][id] << endl;
    cout << "SpatialExtent ymin, ymax = " << sdata[2][id] << " to " << sdata[3][id] << endl;
    cout << "SpatialExtent zmin, zmax = " << sdata[4][id] << " to " << sdata[5][id] << endl;

    for(int i = 0; i < getChildrenLength(); ++i){
        children[i]->addExtent(data, sdata, rdata);
    }
}

// -- DATA IO methods --
void Grid::putQuadMesh(DBfile* file, std::string dataTypeName, const char* coordnames[3], int rankInGroup, DBoptlist* optListMesh, DBoptlist* optListVar){
    const int dim = 3;

    // dimensions
    // glue cellは出力しない
    int dimensions[3];
    dimensions[0] = nx;
    dimensions[1] = ny;
    dimensions[2] = nz;

    // the array of coordinate arrays
    float** coordinates = this->getMeshNodes(dim);

    // quadmesh and quadvar name
    std::string meshname = (format("/block%04d/mesh%04d%04d") % rankInGroup % MPIw::Environment::rank % id).str();
    std::string varname = (format("/block%04d/%s%04d%04d") % rankInGroup % dataTypeName % MPIw::Environment::rank % id).str();

    const char* m = meshname.c_str();
    const char* v = varname.c_str();
    DBPutQuadmesh(file, m, coordnames, coordinates, dimensions, dim, DB_FLOAT, DB_COLLINEAR, optListMesh);

    if(dataTypeName == "potential") {
        const char* varnames[1] = {dataTypeName.c_str()};
        float* tdArray = Utils::getTrueCells(this->getField()->getPhi());
        DBPutQuadvar1(file, v, m, tdArray, dimensions, dim, NULL, 0, DB_FLOAT, DB_NODECENT, optListVar);
        delete [] tdArray;
    } else if(dataTypeName == "rho") {
        float* tdArray = Utils::getTrueCells(this->getField()->getRho());
        DBPutQuadvar1(file, v, m, tdArray, dimensions, dim, NULL, 0, DB_FLOAT, DB_NODECENT, optListVar);
        delete [] tdArray;
    } else if(dataTypeName == "efield") {
        const char* varnames[3] = {"ex", "ey", "ez"};
        float* vars[3];
        vars[0] = Utils::getTrueEdges(this->getField()->getEx(), 0);
        vars[1] = Utils::getTrueEdges(this->getField()->getEy(), 1);
        vars[2] = Utils::getTrueEdges(this->getField()->getEz(), 2);
        DBPutQuadvar(file, v, m, 3, varnames, vars, dimensions, dim, NULL, 0, DB_FLOAT, DB_EDGECENT, NULL);
        delete [] vars[0];
        delete [] vars[1];
        delete [] vars[2];
    } else if(dataTypeName == "bfield") {
        const char* varnames[3] = {"bx", "by", "bz"};
        float* vars[3];
        vars[0] = Utils::getTrueFaces(this->getField()->getBx(), 0);
        vars[1] = Utils::getTrueFaces(this->getField()->getBy(), 1);
        vars[2] = Utils::getTrueFaces(this->getField()->getBz(), 2);
        DBPutQuadvar(file, v, m, 3, varnames, vars, dimensions, dim, NULL, 0, DB_FLOAT, DB_FACECENT, NULL);
        delete [] vars[0];
        delete [] vars[1];
        delete [] vars[2];
    } else {
        throw std::invalid_argument("[ERROR] Invalid argument was passed to putQuadMesh().");
    }

    for(int i = 0; i < children.size(); ++i) {
        children[i]->putQuadMesh(file, dataTypeName, coordnames, rankInGroup, optListMesh, optListVar);
    }

    delete [] coordinates[0];
    delete [] coordinates;
}

Grid::~Grid(){
    //! delete all particles
    //! vector内のparticleは自動でデストラクタが呼ばれる
    particles.erase(particles.begin(), particles.end());

    // reserveしてあった分を削除する
    particles.shrink_to_fit();

    delete field;

    //! delete all children
    std::vector<Grid*>::iterator it = children.begin();
    while(it != children.end()) {
        it = children.erase(it);
    }

    cout << "Grid Destructor Called!" << endl;
}

// --- stdout friend function ----
void printGridInfo(std::ostream& ost, Grid* g, int childnum) {
    std::string tab = "";
    for(int i = 0; i < g->getLevel(); ++i) tab += "    ";

    if(g->getLevel() > 0) ost << tab << "--- child [" << childnum << "] ---" << endl;
    ost << tab << "id: " << g->getID() << endl;
    ost << tab << "level: " << g->getLevel() << endl;
    ost << tab << "dx: " << format("%10.5e") % Utils::Normalizer::unnormalizeLength(g->getDX()) << "m" << endl;
    ost << tab << "nx, ny, nz: " << format("%1%x%2%x%3%") % g->getNX() % g->getNY() % g->getNZ() << " grids [total]" << endl;
    ost << tab << "nx,ny,nz(+): " << format("%1%x%2%x%3%") % (g->getNX() + 2) % (g->getNY() + 2) % (g->getNZ() + 2) << " grids [with glue cells]" << endl;
    ost << tab << "parent from: " << format("%1%,%2%,%3%") % g->getFromIX() % g->getFromIY() % g->getFromIZ() << endl;
    ost << tab << "parent to  : " << format("%1%,%2%,%3%") % g->getToIX() % g->getToIY() % g->getToIZ() << endl;
    ost << tab << "base positions (normalized): " << format("%1%,%2%,%3%") % g->getBaseX() % g->getBaseY() % g->getBaseZ() << endl;
    ost << tab << "sumNumOfChild: " << g->getSumOfChild() << endl;
    ost << tab << "numOfChild: " << g->getChildrenLength() << endl;

    if(g->getChildrenLength() > 0) {
        std::vector<Grid*>& children = g->getChildren();
        for(unsigned int i = 0; i < children.size(); ++i) {
            printGridInfo(ost, children[i], i);
        }
    }
}

std::ostream& operator<<(std::ostream& ost, Grid* g){
    ost << "[Grid Info]" << std::endl;
    printGridInfo(ost, g, 0);
    return ost;
}

