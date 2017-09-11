#include "grid.hpp"
#include "field.hpp"
#include "normalizer.hpp"

// Unique ID の実体
unsigned int Grid::nextID = 0;

// Grid 基底クラス用のコンストラクタ
Grid::Grid(void) : child_map(boost::extents[0][0][0]), field(std::make_unique<Field>()) {
    sumTotalNumOfChildGrids = 0;

    //! UniqueなIDをセット
    id = this->getNextID();
}

//! @note: childrenのエネルギーも取る?
double Grid::getEFieldEnergy(void) const {
    return field->getEfieldEnergy();
}

double Grid::getBFieldEnergy(void) const {
    return field->getBfieldEnergy();
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
    field->getPoissonResidual().resize(tdExtents[cx][cy][cz]);
    field->getPoissonError().resize(tdExtents[cx][cy][cz]);

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

//! ChildMap の resize と 初期化
void Grid::initializeChildMap(void) {
    ChildDefinedMapInt::extent_gen mapExtentGen;
    child_map.resize(mapExtentGen[nx + 2][ny + 2][nz + 2]);

    for(int i = 0; i < nx + 2; ++i) {
        for(int j = 0; j < ny + 2; ++j) {
            for(int k = 0; k < nz + 2; ++k) {
                child_map[i][j][k] = CHILD_MAP_TAG::NOT_EXIST;
            }
        }
    }
}

Grid::~Grid(){
    //! delete all particles
    if (particles.size() > 0) {
        particles.erase(particles.begin(), particles.end());

        // reserveしてあった分を削除する
        particles.shrink_to_fit();
    }

    cout << "Grid Destructor Called! " << id << endl;
}

void Grid::printInfo() const {
    std::string tab = "";
    for(int i = 0; i < level; ++i) tab += "    ";
    if (level == 0) {
        cout << "[Grid Info]" << endl;
    } else {
        cout << tab << "--- child [" << id << "] ---" << endl;
    }
    cout << tab << "id: " << id << endl;
    cout << tab << "level: " << level << endl;
    cout << tab << "dx: " << format("%10.5e") % Normalizer::unnormalizeLength(dx) << "m" << endl;
    cout << tab << "dt: " << format("%10.5e") % Normalizer::unnormalizeTime(dt) << "m" << endl;
    cout << tab << "nx, ny, nz: " << format("%1%x%2%x%3%") % nx % ny % nz << " grids [total]" << endl;
    cout << tab << "nx,ny,nz(+): " << format("%1%x%2%x%3%") % (nx + 2) % (ny + 2) % (nz + 2) << " grids [with glue cells]" << endl;
    cout << tab << "parent from: " << format("%1%,%2%,%3%") % from_ix % from_iy % from_iz << endl;
    cout << tab << "parent to  : " << format("%1%,%2%,%3%") % to_ix % to_iy % to_iz << endl;
    cout << tab << "base positions (normalized): " << format("%1%,%2%,%3%") % base_x % base_y % base_z << endl;
    cout << tab << "sumNumOfChild: " << this->getSumOfChild() << endl;
    cout << tab << "numOfChild: " << this->getChildrenLength() << endl;

    for(const auto& child : children) {
        child->printInfo();
    }
}
