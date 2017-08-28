#include "global.hpp"
#include "position.hpp"
#include "grid.hpp"
#include "field.hpp"
#include "particle.hpp"
#include "environment.hpp"
#include "mpiw.hpp"
#include "utils.hpp"
#include "normalizer.hpp"
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

    auto& rho = field->getRho();
    for(int i = 0; i < Environment::num_of_particle_types + 1; ++i) {
        //! 総和のtdArray + 粒子種毎のtdArrayを直接配置で生成する
        rho.emplace_back(tdExtents[cx][cy][cz], boost::fortran_storage_order());
    }

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

    for (const auto& object_info : Environment::objects_info) {
        std::string obj_name = object_info.name;
        //! 物体関連の設定を関連付けされた obj 形式ファイルから読み込む
        ObjectDataFromFile object_data = ObjectUtils::getObjectNodesFromObjFile(object_info.file_name);

        unsigned int num_cmat = object_data.nodes.size();
        const auto& node_array = object_data.nodes;
        const auto& face_array = object_data.faces;

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
                inner_node_array[cmat_itr] = getRelativePosition<int>(i, j, k);
            } else if (isGlueNode(i, j, k)) {
                glue_node_array[cmat_itr] = getRelativePosition<int>(i, j, k);
            }
        }

        //! 面判定
        ObjectFaces inner_face_array;
        for(const auto& face_values : face_array) {
            const auto face_type = face_values[0];
            const auto i = face_values[1];
            const auto j = face_values[2];
            const auto k = face_values[3];

            if (isInnerFaceWithGlue(face_type, i, j, k)) {
                const auto& rel = getRelativePosition<int>(i, j, k);
                inner_face_array.push_back({{face_type, rel[0], rel[1], rel[2]}});
            }
        }

        //! 物体定義点がゼロでも Spacecraft オブジェクトだけは作成しておいた方がよい
        //! emplace_back で Spacecraft object を直接構築
        objects.emplace_back(nx, ny, nz, num_cmat, object_info, inner_node_array, glue_node_array, inner_face_array);

        //! Comm作成 (物体が入っていないならnullになる)
        MPIw::Environment::makeNewComm(obj_name, is_object_in_this_node);
        if (MPIw::Environment::isRootNode(obj_name)) {
            cout << Environment::rankStr() << "is set to Root Node for " << obj_name << "." << endl;
            cout << objects[ objects.size() - 1 ] << endl;
        }
    }
}

void Grid::initializeObjectsCmatrix(void) {
    if (Environment::isRootNode) cout << "-- Initializing Objects Capacity Matrix --" << endl;
    RhoArray& rho = field->getRho();
    tdArray& phi = field->getPhi();

    for(auto& obj : objects) {
        const auto num_cmat = obj.getCmatSize();

        { //! Progress Manager のライフタイムを区切る
            Utils::ProgressManager pm(num_cmat, "cmat_solve");

            for(unsigned int cmat_col_itr = 0; cmat_col_itr < num_cmat; ++cmat_col_itr ) {
                if (Environment::isRootNode) pm.update(cmat_col_itr);

                // rhoを初期化
                Utils::initializeRhoArray(rho);
                Utils::initialize3DArray(phi);

                //! 各頂点に単位電荷を付与
                if (obj.isMyCmat(cmat_col_itr)) {
                    const auto& cmat_pos = obj.getCmatPos(cmat_col_itr);
                    rho[0][cmat_pos.i][cmat_pos.j][cmat_pos.k] = 1.0;
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
    if(!Environment::isNotBoundary(AXIS::x, AXIS_SIDE::up)) max_x -= 1.0;
    if(!Environment::isNotBoundary(AXIS::y, AXIS_SIDE::up)) max_y -= 1.0;
    if(!Environment::isNotBoundary(AXIS::z, AXIS_SIDE::up)) max_z -= 1.0;

    //! - particlesは空のstd::vector< std::vector<Particle> >として宣言されている
    //! - particle types 分だけresize
    particles.resize(Environment::num_of_particle_types);

    for(int id = 0; id < Environment::num_of_particle_types; ++id){
        //! 各粒子分のメモリをreserveしておく
        particles[id].reserve(Environment::max_particle_num);

        //! 初期化時は背景粒子のみ生成
        if (Environment::getParticleType(id)->getType() == "ambient") {
            auto ambient_particle_ptr = Environment::getAmbientParticleType(id);
            int pnum = ambient_particle_ptr->getTotalNumber();

            for(int i = 0; i < pnum; ++i){
                Particle p = ambient_particle_ptr->generateNewParticle(0.0, max_x, 0.0, max_y, 0.0, max_z);

                //! 物体がある場合は生成時にチェックする
                for(const auto& obj : objects) {
                    if (obj.isDefined()) obj.removeInnerParticle(p);
                }

                if (p.isValid) particles[id].push_back( std::move(p) );
            }
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
    return (level == 0 && Environment::isNotBoundary(AXIS::x, AXIS_SIDE::up)) ? nx + 1 : nx;
}

int Grid::getYNodeSize(void) const {
    return (level == 0 && Environment::isNotBoundary(AXIS::y, AXIS_SIDE::up)) ? ny + 1 : ny;
}

int Grid::getZNodeSize(void) const {
    return (level == 0 && Environment::isNotBoundary(AXIS::z, AXIS_SIDE::up)) ? nz + 1 : nz;
}

//! for DATA IO
//! 粒子の位置から密度を計算する
boost::multi_array<float, 3> Grid::getDensity(const int pid) const {
    // ZONECENTなので-1する
    const int xsize = this->getXNodeSize() - 1;
    const int ysize = this->getYNodeSize() - 1;
    const int zsize = this->getZNodeSize() - 1;

    boost::multi_array<float, 3> zones(boost::extents[xsize][ysize][zsize], boost::fortran_storage_order());
    const auto size = static_cast<float>(Normalizer::unnormalizeDensity(Environment::getParticleType(pid)->getSize()));

    for(int pnum = 0; pnum < particles[pid].size(); ++pnum){
        const Particle& p = particles[pid][pnum];

        if(p.isValid) {
            Position pos(p);
            zones[pos.i - 1][pos.j - 1][pos.k - 1] += size;
        }
    }

    // RVO
    return zones;
}

//! Glueノードを含まないデータを生成
boost::multi_array<float, 3> Grid::getTrueNodes(const tdArray& x3D, const double unnorm) const {
    int xsize = this->getXNodeSize();
    int ysize = this->getYNodeSize();
    int zsize = this->getZNodeSize();

    boost::multi_array<float, 3> true_nodes(boost::extents[xsize][ysize][zsize], boost::fortran_storage_order());

    for(int k = 1; k < zsize + 1; ++k){
        for(int j = 1; j < ysize + 1; ++j){
            for(int i = 1; i < xsize + 1; ++i){
                true_nodes[i - 1][j - 1][k - 1] = static_cast<float>(x3D[i][j][k] * unnorm);
            }
        }
    }

    // RVO
    return true_nodes;
}

//! RhoArray用
boost::multi_array<float, 3> Grid::getTrueNodes(const RhoArray& rho, const int pid, const double unnorm) const {
    int xsize = this->getXNodeSize();
    int ysize = this->getYNodeSize();
    int zsize = this->getZNodeSize();

    boost::multi_array<float, 3> true_nodes(boost::extents[xsize][ysize][zsize], boost::fortran_storage_order());

    for(int k = 1; k < zsize + 1; ++k){
        for(int j = 1; j < ysize + 1; ++j){
            for(int i = 1; i < xsize + 1; ++i){
                true_nodes[i - 1][j - 1][k - 1] = static_cast<float>(rho[pid][i][j][k] * unnorm);
            }
        }
    }

    // RVO
    return true_nodes;
}

// -- DATA IO methods --
void Grid::putFieldData(HighFive::Group& group, const std::string& data_type_name, const std::string& i_timestamp) {
    auto getGroup = [](auto& g, const std::string& group_name) {
        if (g.exist(group_name)) {
            return g.getGroup(group_name);
        } else {
            return g.createGroup(group_name);
        }
    };

    const std::string& level_str = "level" + std::to_string(level);
    HighFive::Group local_group = getGroup(group, level_str);

    if (data_type_name == "potential") {
        auto values = this->getTrueNodes(field->getPhi(), Normalizer::unnormalizePotential(1.0));
        auto dataset = local_group.createDataSet<float>(data_type_name, HighFive::DataSpace::From(values));
        dataset.write(values);
    } else if(data_type_name == "rho") {
        auto values = this->getTrueNodes(field->getRho(), 0, Normalizer::unnormalizeRho(1.0));
        auto dataset = local_group.createDataSet<float>(data_type_name, HighFive::DataSpace::From(values));
        dataset.write(values);
    } else if(data_type_name == "efield") {
        HighFive::Group data_type_group = getGroup(local_group, data_type_name);

        const std::array<std::string, 3> axis{{"ex", "ey", "ez"}};
        for(const auto& group_name : axis) {
            if (group_name == "ex") {
                auto values = this->getTrueNodes(field->getExRef(), Normalizer::unnormalizeEfield(1.0));
                auto dataset = local_group.createDataSet<float>(group_name, HighFive::DataSpace::From(values));
                dataset.write(values);
            } else if (group_name == "ey") {
                auto values = this->getTrueNodes(field->getEyRef(), Normalizer::unnormalizeEfield(1.0));
                auto dataset = local_group.createDataSet<float>(group_name, HighFive::DataSpace::From(values));
                dataset.write(values);
            } else if (group_name == "ez") {
                auto values = this->getTrueNodes(field->getEzRef(), Normalizer::unnormalizeEfield(1.0));
                auto dataset = local_group.createDataSet<float>(group_name, HighFive::DataSpace::From(values));
                dataset.write(values);
            }
        }
    } else if(data_type_name == "bfield") {
        HighFive::Group data_type_group = getGroup(local_group, data_type_name);
        const std::array<std::string, 3> axis{{"bx", "by", "bz"}};

        for(const auto& group_name : axis) {
            if (group_name == "bx") {
                auto values = this->getTrueNodes(field->getBxRef(), Normalizer::unnormalizeBfield(1.0));
                auto dataset = local_group.createDataSet<float>(group_name, HighFive::DataSpace::From(values));
                dataset.write(values);
            } else if (group_name == "by") {
                auto values = this->getTrueNodes(field->getByRef(), Normalizer::unnormalizeBfield(1.0));
                auto dataset = local_group.createDataSet<float>(group_name, HighFive::DataSpace::From(values));
                dataset.write(values);
            } else if (group_name == "bz") {
                auto values = this->getTrueNodes(field->getBzRef(), Normalizer::unnormalizeBfield(1.0));
                auto dataset = local_group.createDataSet<float>(group_name, HighFive::DataSpace::From(values));
                dataset.write(values);
            }
        }
    } else if(data_type_name == "density") {
        HighFive::Group data_type_group = getGroup(local_group, data_type_name);
        for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
            const std::string& pname = Environment::getParticleType(pid)->getName();
            auto values = this->getDensity(pid);
            auto dataset = data_type_group.createDataSet<float>(pname, HighFive::DataSpace::From(values));
            dataset.write(values);
        }
    }

    for(int i = 0; i < children.size(); ++i) {
        children[i]->putFieldData(group, data_type_name, i_timestamp);
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

