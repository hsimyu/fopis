#include "global.hpp"
#include "position.hpp"
#include "grid.hpp"
#include "field.hpp"
#include "particle.hpp"
#include "environment.hpp"
#include "mpiw.hpp"
#include "utils.hpp"
#include "normalizer.hpp"
#include <silo.h>
#include <random>
#include <algorithm>

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

//! 場の resize を行う
void Grid::initializeField(void){
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

void Grid::initializeObject(void) {
    //! TODO: オブジェクト数を数える
    if (Environment::isRootNode) cout << "-- Defining Objects -- " << endl;

    std::map<std::string, ObjectNodes> defined_objects;
    std::map<std::string, unsigned int> num_cmat_map;

    for (const auto& object_info : Environment::objects_info) {
        std::string obj_name = object_info.name;
        //! 物体関連の設定を関連付けされた obj 形式ファイルから読み込む
        defined_objects[obj_name] = Utils::getObjectNodesFromObjFile(object_info.file_name);
        num_cmat_map[obj_name] = defined_objects[obj_name].size();
    }

    for(const auto& object_nodes : defined_objects) {
        const auto obj_name = object_nodes.first;
        const auto& node_array = object_nodes.second;

        //! innerと判定されたやつだけ渡す
        ObjectNodes inner_node_array;
        ObjectNodes glue_node_array;
        bool is_object_in_this_node = false;

        for(const auto& node_pair : node_array) {
            const auto cmat_itr = node_pair.first;
            const auto& node_pos = node_pair.second;

            const auto i = node_pos[0];
            const auto j = node_pos[1];
            const auto k = node_pos[2];

            if (isInnerNode(i, j, k)) {
                is_object_in_this_node = true;
                inner_node_array[cmat_itr] = getRelativePosition<unsigned int>(i, j, k);
            } else if (isGlueNode(i, j, k)) {
                glue_node_array[cmat_itr] = getRelativePosition<unsigned int>(i, j, k);
            }
        }

        //! Comm作成 (物体が入っていないならnullになる)
        MPIw::Environment::makeNewComm(obj_name, is_object_in_this_node);
        // if (is_object_in_this_node) {
        //     cout << Environment::rankStr() << "Object " << obj_name << " was defined in me." << endl;
        // }

        //! 物体定義点がゼロでも Spacecraft オブジェクトだけは作成しておいた方がよい
        Spacecraft spc(nx, ny, nz, num_cmat_map[obj_name], obj_name, inner_node_array, glue_node_array);
        if (Environment::isRootNode) cout << spc << endl;
        objects.emplace_back( spc );
    }
}

void Grid::initializeObjectsCmatrix(void) {
    if (Environment::isRootNode) cout << "-- Initializing Objects Capacity Matrix --" << endl;
    tdArray& rho = field->getRho();
    tdArray& phi = field->getPhi();

    for(auto& obj : objects) {
        const auto num_cmat = obj.getCmatSize();

        { //! Progress Manager のライフタイムを区切る
            Utils::ProgressManager pm(num_cmat, "cmat_solve");

            for(unsigned int cmat_col_itr = 0; cmat_col_itr < num_cmat; ++cmat_col_itr ) {
                if (Environment::isRootNode) pm.update(cmat_col_itr);

                // rhoを初期化
                Utils::initializeTdarray(rho);
                Utils::initializeTdarray(phi);

                //! 各頂点に単位電荷を付与
                if (obj.isMyCmat(cmat_col_itr)) {
                    const auto& cmat_pos = obj.getCmatPos(cmat_col_itr);
                    rho[cmat_pos.i][cmat_pos.j][cmat_pos.k] = 1.0;
                }

                solvePoisson();

                for(unsigned int cmat_row_itr = 0; cmat_row_itr < num_cmat; ++cmat_row_itr ) {
                    double value = 0.0;
                    if (obj.isMyCmat(cmat_row_itr)) {
                        //! phiの値がB_{ij}の値になっている
                        const auto& target_pos = obj.getCmatPos(cmat_row_itr);
                        value = phi[target_pos.i][target_pos.j][target_pos.k];
                    }

                    if (obj.isDefined()) {
                        //! bcastの代わりにsumしてしまう
                        value = MPIw::Environment::Comms[obj.getName()].sum(value);
                        obj.setCmatValue(cmat_col_itr, cmat_row_itr, value);
                    }
                }
            }

            //! 物体が有効でないなら解く必要なし
            if (obj.isDefined()) obj.makeCmatrixInvert();
        }
    }
}

// root grid constructor
Grid::Grid(void) : field(std::make_unique<Field>()) {
    //! - コンストラクタにEnvironmentクラスが渡された場合、
    //! レベル0のGridを作成します.
    level = 0;
    sumTotalNumOfChildGrids = 0;

    //! UniqueなIDをセット
    id = this->getNextID();

    nx = Environment::cell_x;
    ny = Environment::cell_y;
    nz = Environment::cell_z;
    dx = Normalizer::normalizeLength(Environment::dx);
    dt = Normalizer::normalizeTime(Environment::dt);

    //! @{
    //! Root Gridの場合の親グリッドは、計算空間を全て統合した空間として、
    //! その上にプロセス分割されたグリッドが乗っていると考える
    from_ix = Environment::getAssignedXBegin();
    from_iy = Environment::getAssignedYBegin();
    from_iz = Environment::getAssignedZBegin();
    to_ix = Environment::getAssignedXEnd();
    to_iy = Environment::getAssignedYEnd();
    to_iz = Environment::getAssignedZEnd();
    //! @note: base_x, base_y, base_zは正規化された長さ
    base_x = dx * static_cast<double>(from_ix);
    base_y = dx * static_cast<double>(from_iy);
    base_z = dx * static_cast<double>(from_iz);
    //! @}

    // Field初期化
    this->initializeField();

    // 物体初期化
    this->initializeObject();
    this->initializeObjectsCmatrix();

    //! - 粒子位置の上限を設定
    double max_x = static_cast<double>(Environment::cell_x);
    double max_y = static_cast<double>(Environment::cell_y);
    double max_z = static_cast<double>(Environment::cell_z);

    //! - 上側境界にいる場合は外側にはみ出した粒子を生成しないようにする
    if(!Environment::isPeriodic(AXIS::x, AXIS_SIDE::up)) max_x -= 1.0;
    if(!Environment::isPeriodic(AXIS::y, AXIS_SIDE::up)) max_y -= 1.0;
    if(!Environment::isPeriodic(AXIS::z, AXIS_SIDE::up)) max_z -= 1.0;

    //! - particlesは空のstd::vector< std::vector<Particle> >として宣言されている
    //! - particle types 分だけresize
    particles.resize(Environment::num_of_particle_types);

    for(int id = 0; id < Environment::num_of_particle_types; ++id){
        int pnum = Environment::ptype[id].getTotalNumber();

        //! 各粒子分のメモリをreserveしておく
        particles[id].reserve(Environment::max_particle_num);

        for(int i = 0; i < pnum; ++i){
            Particle p(id);
            p.generateNewPosition(0.0, max_x, 0.0, max_y, 0.0, max_z);
            p.generateNewVelocity();

            //! 物体がある場合は生成時にチェックする
            for(const auto& obj : objects) {
                if (obj.isDefined()) obj.removeInnerParticle(p);
            }

            if (p.isValid) particles[id].emplace_back(p);
        }
    }
}

//! child grid constructor
//! GridコンストラクタにGridが渡された場合、
//! そのGridを親とした子グリッドを生成します
Grid::Grid(Grid* g,
        const int _from_ix, const int _from_iy, const int _from_iz,
        const int _to_ix,   const int _to_iy,   const int _to_iz) : field(std::make_unique<Field>())
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

    if(!isValid) MPIw::Environment::abort(1);
}

//! 粒子の位置から密度を計算する
float* Grid::getDensity(const int pid) {
    // ZONECENTなので-1する
    const int xsize = this->getXNodeSize() - 1;
    const int ysize = this->getYNodeSize() - 1;
    const int zsize = this->getZNodeSize() - 1;
    float* zone_density = new float[xsize * ysize * zsize];
    const int maxitr = xsize*ysize*zsize;

    //! initialize
    for(int i = 0; i < maxitr; ++i){
        zone_density[i] = 0.0f;
    }

    const auto size = static_cast<float>(Normalizer::unnormalizeDensity(Environment::ptype[pid].getSize()));

    for(int pnum = 0; pnum < particles[pid].size(); ++pnum){
        Particle& p = particles[pid][pnum];

        if(p.isValid) {
            Position pos(p);
            int i = pos.i - 1; // 対応する zone 番号へ変換
            int j = pos.j - 1;
            int k = pos.k - 1;

            int itr = i + xsize*j + xsize*ysize*k;
            zone_density[itr] += size;
        }
    }

    return zone_density;
}

//! @note: childrenのエネルギーも取る?
double Grid::getEFieldEnergy(void) const {
    return field->getEfieldEnergy();
}

double Grid::getBFieldEnergy(void) const {
    return field->getBfieldEnergy();
}

// 子グリッドへ場の値をコピーする
// この実装はノード to ノードの場合
void Grid::copyScalarToChildren(std::string varname){
    tdArray& tdValue = field->getScalar(varname);

    // @note: OpenMP
    for(int chidx = 0; chidx < children.size(); ++chidx) {
        tdArray& childValue = children[chidx]->getScalar(varname);

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

int Grid::getXNodeSize(void) const {
    //! 周期境界の場合は上側境界と下側境界の間の空間も有効な空間となるので、
    //! 上側のノードを1つ増やす
    return (level == 0 && Environment::isPeriodic(AXIS::x, AXIS_SIDE::up)) ? nx + 1 : nx;
}

int Grid::getYNodeSize(void) const {
    return (level == 0 && Environment::isPeriodic(AXIS::y, AXIS_SIDE::up)) ? ny + 1 : ny;
}

int Grid::getZNodeSize(void) const {
    return (level == 0 && Environment::isPeriodic(AXIS::z, AXIS_SIDE::up)) ? nz + 1 : nz;
}

// mesh nodesの座標配列を生成
float** Grid::getMeshNodes(int dim) {
    // the array of coordinate arrays
    // @note: メモリリーク防止のため必ずdeleteする
    const float real_dx = Normalizer::unnormalizeLength(dx);
    const float real_base_x = Normalizer::unnormalizeLength(base_x);
    const float real_base_y = Normalizer::unnormalizeLength(base_y);
    const float real_base_z = Normalizer::unnormalizeLength(base_z);

    float** coordinates = new float*[dim];
    // root にいて正方向の端でない場合は+1分出力する
    int xsize = this->getXNodeSize();
    coordinates[0] = new float[xsize];

    for(int i = 0; i < xsize; ++i) {
        coordinates[0][i] = real_base_x + real_dx * i;
    }

    int ysize = this->getYNodeSize();
    coordinates[1] = new float[ysize];

    for(int i = 0; i < ysize; ++i) {
        coordinates[1][i] = real_base_y + real_dx * i;
    }

    int zsize = this->getZNodeSize();
    coordinates[2] = new float[zsize];

    for(int i = 0; i < zsize; ++i) {
        coordinates[2][i] = real_base_z + real_dx * i;
    }
    return coordinates;
}

//! for DATA IO
float* Grid::getTrueNodes(const tdArray& x3D){
    int xsize = this->getXNodeSize();
    int ysize = this->getYNodeSize();
    int zsize = this->getZNodeSize();
    float* x1D = new float[xsize*ysize*zsize];

    //! C-based indicing
    //! 上側境界にいない時 or 周期境界である時、上側のglue cellの値も出力する
    for(int k = 1; k < zsize + 1; ++k){
        for(int j = 1; j < ysize + 1; ++j){
            for(int i = 1; i < xsize + 1; ++i){
                x1D[(i-1) + (j-1)*xsize + (k-1)*xsize*ysize] = static_cast<float>(x3D[i][j][k]);
            }
        }
    }

    return x1D;
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
    dimensions[0] = this->getXNodeSize();
    dimensions[1] = this->getYNodeSize();
    dimensions[2] = this->getZNodeSize();

    // the array of coordinate arrays
    float** coordinates = this->getMeshNodes(dim);

    // quadmesh and quadvar name
    std::string meshname = (format("/block%04d/mesh%04d%04d") % rankInGroup % MPIw::Environment::rank % id).str();
    std::string varname = (format("/block%04d/%s%04d%04d") % rankInGroup % dataTypeName % MPIw::Environment::rank % id).str();

    const char* m = meshname.c_str();
    const char* v = varname.c_str();
    DBPutQuadmesh(file, m, coordnames, coordinates, dimensions, dim, DB_FLOAT, DB_COLLINEAR, optListMesh);

    delete [] coordinates[0];
    delete [] coordinates[1];
    delete [] coordinates[2];
    delete [] coordinates;

    if(dataTypeName == "potential") {
        float* tdArray = this->getTrueNodes(field->getPhi());
        DBPutQuadvar1(file, v, m, tdArray, dimensions, dim, NULL, 0, DB_FLOAT, DB_NODECENT, optListVar);
        delete [] tdArray;
    } else if(dataTypeName == "rho") {
        float* tdArray = this->getTrueNodes(field->getRho());
        DBPutQuadvar1(file, v, m, tdArray, dimensions, dim, NULL, 0, DB_FLOAT, DB_NODECENT, optListVar);
        delete [] tdArray;
    } else if(dataTypeName == "efield") {
        // const char* varnames[3] = {"ex", "ey", "ez"};
        // float* vars[3];
        // vars[0] = Utils::getTrueEdges(this->getField()->getEx(), 0);
        // vars[1] = Utils::getTrueEdges(this->getField()->getEy(), 1);
        // vars[2] = Utils::getTrueEdges(this->getField()->getEz(), 2);
        // DBPutQuadvar(file, v, m, 3, varnames, vars, dimensions, dim, NULL, 0, DB_FLOAT, DB_EDGECENT, optListVar);
        // delete [] vars[0];
        // delete [] vars[1];
        // delete [] vars[2];
        const char* varnames[3] = {"ex", "ey", "ez"};
        float* vars[3];
        vars[0] = getTrueNodes(field->getExRef());
        vars[1] = getTrueNodes(field->getEyRef());
        vars[2] = getTrueNodes(field->getEzRef());
        DBPutQuadvar(file, v, m, 3, varnames, vars, dimensions, dim, NULL, 0, DB_FLOAT, DB_NODECENT, optListVar);
        delete [] vars[0];
        delete [] vars[1];
        delete [] vars[2];
    } else if(dataTypeName == "bfield") {
        const char* varnames[3] = {"bx", "by", "bz"};
        float* vars[3];
        vars[0] = getTrueNodes(field->getBxRef());
        vars[1] = getTrueNodes(field->getByRef());
        vars[2] = getTrueNodes(field->getBzRef());
        DBPutQuadvar(file, v, m, 3, varnames, vars, dimensions, dim, NULL, 0, DB_FLOAT, DB_NODECENT, optListVar);
        delete [] vars[0];
        delete [] vars[1];
        delete [] vars[2];
    } else if(dataTypeName == "density") {
        // zone centに変更
        dimensions[0] -= 1;
        dimensions[1] -= 1;
        dimensions[2] -= 1;
        
        const auto num_of_ptypes = Environment::num_of_particle_types;

        if(num_of_ptypes == 1) {
            float* tdArray = this->getDensity(0);
            const char* vname = Environment::ptype[0].getName().c_str();
            DBPutQuadvar1(file, vname, m, tdArray, dimensions, dim, NULL, 0, DB_FLOAT, DB_ZONECENT, optListVar);
            delete [] tdArray;
        } else {
            float** vars = new float*[num_of_ptypes];
            char** vnames = new char*[num_of_ptypes];
            for(int pid = 0; pid < num_of_ptypes; ++pid){
                std::string pname = Environment::ptype[pid].getName();
                vnames[pid] = new char[pname.size() + 1];
                std::strcpy(vnames[pid], pname.c_str());

                // 密度配列の取得
                vars[pid] = this->getDensity(pid);
            }

            DBPutQuadvar(file, v, m, num_of_ptypes, vnames, vars, dimensions, dim, NULL, 0, DB_FLOAT, DB_ZONECENT, optListVar);

            for(int pid = 0; pid < num_of_ptypes; ++pid){
                delete [] vars[pid];
                delete [] vnames[pid];
            }
            delete [] vnames;
            delete [] vars;
        }
    } else {
        throw std::invalid_argument("[ERROR] Invalid argument was passed to putQuadMesh().");
    }
    
    for(int i = 0; i < children.size(); ++i) {
        children[i]->putQuadMesh(file, dataTypeName, coordnames, rankInGroup, optListMesh, optListVar);
    }

}

Grid::~Grid(){
    //! delete all particles
    //! vector内のparticleは自動でデストラクタが呼ばれる
    particles.erase(particles.begin(), particles.end());

    // reserveしてあった分を削除する
    particles.shrink_to_fit();

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
    ost << tab << "dx: " << format("%10.5e") % Normalizer::unnormalizeLength(g->getDx()) << "m" << endl;
    ost << tab << "dt: " << format("%10.5e") % Normalizer::unnormalizeTime(g->getDt()) << "m" << endl;
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

