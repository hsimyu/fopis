#include "global.hpp"
#include "position.hpp"
#include "grid.hpp"
#include "field.hpp"
#include "particle.hpp"
#include "environment.hpp"
#include "mpiw.hpp"
#include "utils.hpp"
#include "normalizer.hpp"
#include "dataio.hpp"
#include <random>
#include <algorithm>

// Unique ID の実体
unsigned int Grid::nextID = 0;

// Grid 基底クラス用のコンストラクタ
Grid::Grid(void) : field(std::make_unique<Field>()) {
    sumTotalNumOfChildGrids = 0;

    //! UniqueなIDをセット
    id = this->getNextID();
}

void Grid::makeChild(const int _from_ix, const int _from_iy, const int _from_iz, const int _to_ix, const int _to_iy, const int _to_iz) {
    this->addChild(
        std::make_unique<ChildGrid>(this, _from_ix, _from_iy, _from_iz, _to_ix, _to_iy, _to_iz)
    );
    incrementSumOfChild();
}

void Grid::addChild(std::unique_ptr<ChildGrid>&& child) {
    children.push_back( std::move(child) );
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
void Grid::putFieldData(HighFive::Group& group, const std::string& data_type_name, const std::string& i_timestamp) const {
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
    auto it = children.begin();
    while(it != children.end()) {
        it = children.erase(it);
    }

    cout << "Grid Destructor Called!" << endl;
}

// --- stdout friend function ----
void printGridInfo(std::ostream& ost, std::shared_ptr<Grid> g, int childnum) {
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

    /*if(g->getChildrenLength() > 0) {
        std::vector<Grid*>& children = g->getChildren();
        for(unsigned int i = 0; i < children.size(); ++i) {
            printGridInfo(ost, children[i], i);
        }
    }*/
}

std::ostream& operator<<(std::ostream& ost, std::shared_ptr<Grid> g){
    ost << "[Grid Info]" << std::endl;
    printGridInfo(ost, g, 0);
    return ost;
}

