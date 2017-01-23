#include <tdpic.h>
#include <random>

// Unique ID
unsigned int Grid::nextID = 0;
unsigned int Grid::getNextID(void) {
    return Grid::nextID++;
}

void Grid::resetNextID(void) {
    nextID = 0;
}

// accessors
void   Grid::setFromIX(int _ix) { from_ix = _ix; }
void   Grid::setFromIY(int _iy) { from_iy = _iy; }
void   Grid::setFromIZ(int _iz) { from_iz = _iz; }
int    Grid::getFromIX(void) const { return from_ix; }
int    Grid::getFromIY(void) const { return from_iy; }
int    Grid::getFromIZ(void) const { return from_iz; }

void   Grid::setToIX(int _ix) { to_ix = _ix; }
void   Grid::setToIY(int _iy) { to_iy = _iy; }
void   Grid::setToIZ(int _iz) { to_iz = _iz; }
int    Grid::getToIX(void) const { return to_ix; }
int    Grid::getToIY(void) const { return to_iy; }
int    Grid::getToIZ(void) const { return to_iz; }

void   Grid::setBaseX(double _x){ base_x = _x; }
void   Grid::setBaseY(double _y){ base_y = _y; }
void   Grid::setBaseZ(double _z){ base_z = _z; }
double Grid::getBaseX(void)  const { return base_x; }
double Grid::getBaseY(void)  const { return base_y; }
double Grid::getBaseZ(void)  const { return base_z; }

unsigned int Grid::getID(void) const { return id; }

void Grid::setNX(int _x){ nx = _x; }
void Grid::setNY(int _y){ ny = _y; }
void Grid::setNZ(int _z){ nz = _z; }
int  Grid::getNX(void) const { return nx; }
int  Grid::getNY(void) const { return ny; }
int  Grid::getNZ(void) const { return nz; }

void Grid::setLevel(int l){ level = l; }
int  Grid::getLevel(void) const { return level; }
void   Grid::setDX(double _dx){ dx = _dx; }
double Grid::getDX(void) const { return dx; }

void  Grid::setParent(Grid* g){ parent = g; }
Grid* Grid::getParent(void){ return parent; }

void Grid::makeChild(
        const int from_ix, const int from_iy, const int from_iz,
        const int to_ix,   const int to_iy,   const int to_iz)
{
    Grid* child = new Grid(this, from_ix, from_iy, from_iz, to_ix, to_iy, to_iz);

    this->addChild(child);
    incrementSumOfChild();
}

void Grid::addChild(Grid* child) { children.push_back(child); }
std::vector<Grid*>& Grid::getChildren(void) {
    // 参照にしないと新しいポインタが生まれてしまう？
    return children;
}

int Grid::getChildrenLength(void) const {
    return children.size();
}

int Grid::getSumOfChild(void) const {
    return sumTotalNumOfChildGrids;
}

void Grid::setSumOfChild(const int s) {
    sumTotalNumOfChildGrids = s;
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

    field->getEx().resize(tdExtents[cx-1][cy][cz]);
    field->getEy().resize(tdExtents[cx][cy-1][cz]);
    field->getEz().resize(tdExtents[cx][cy][cz-1]);

    field->getBx().resize(tdExtents[cx][cy-1][cz-1]);
    field->getBy().resize(tdExtents[cx-1][cy][cz-1]);
    field->getBz().resize(tdExtents[cx-1][cy-1][cz]);

    this->setField(field);
}

// root grid constructor
Grid::Grid(const Environment* env){
    //! - コンストラクタにEnvironmentクラスが渡された場合、
    //! レベル0のGridを作成します.
    level = 0;
    sumTotalNumOfChildGrids = 0;

    //! UniqueなIDをセット
    id = this->getNextID();

    nx = env->cell_x;
    ny = env->cell_y;
    nz = env->cell_z;
    dx = env->dx;

    //! @{
    //! Root Gridの場合の親グリッドは、計算空間を全て統合した空間として、
    //! その上にプロセス分割されたグリッドが乗っていると考える
    from_ix = env->xrank * env->cell_x;
    from_iy = env->yrank * env->cell_y;
    from_iz = env->zrank * env->cell_z;
    to_ix = from_ix + nx;
    to_iy = from_iy + ny;
    to_iz = from_iz + nz;
    base_x = dx * static_cast<double>(from_ix);
    base_y = dx * static_cast<double>(from_iy);
    base_z = dx * static_cast<double>(from_iz);
    //! @}

    // Field初期化
    this->initializeField();

    //! 粒子位置の上限を設定
    //! [0, max_x)になるよう1e-20を引いておく
    const double max_x = static_cast<double>(env->cell_x) - 1e-20;
    const double max_y = static_cast<double>(env->cell_y) - 1e-20;
    const double max_z = static_cast<double>(env->cell_z) - 1e-20;

    // std::random_device rnd;
    const int random_src_x = 10684930;
    const int random_src_y = 99881;
    const int random_src_z = 861200045;
    const int random_src_vx = 930;
    const int random_src_vy = 98076621;
    const int random_src_vz = 7662566;
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
    particles.resize(env->num_of_particle_types);

    for(int id = 0; id < env->num_of_particle_types; ++id){
        int pnum = env->ptype[id].getTotalNumber();
        //! particle_number分のコンストラクタが呼ばれる
        particles[id].resize(pnum);

        const double deviation = Utils::Normalizer::normalizeVelocity( env->ptype[id].calcDeviation() );
        std::normal_distribution<> dist_vx(0.0, deviation);
        std::normal_distribution<> dist_vy(0.0, deviation);
        std::normal_distribution<> dist_vz(0.0, deviation);

        //! - 粒子はレベル0グリッドにのみ所属します
        for(int i = 0; i < pnum; ++i){
            particles[id][i].setPosition(dist_x(mt_x), dist_y(mt_y), dist_z(mt_z));
            particles[id][i].setVelocity(dist_vx(mt_vx), dist_vy(mt_vy), dist_vz(mt_vz));
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
    if( to_ix > parent->getNX() ){
        std::cerr << "[ERROR] A child patch's x-extent exceeds the parent's extent. : " << to_ix << " > " << parent->getNX() << endl;
        isValid = false;
    }

    // y extent
    if( to_iy > parent->getNY() ){
        std::cerr << "[ERROR] A child patch's y-extent exceeds the parent's extent. : " << to_iy << " > " << parent->getNY() << endl;
        isValid = false;
    }

    // z extent
    if( to_iz > parent->getNZ() ){
        std::cerr << "[ERROR] A child patch's z-extent exceeds the parent's extent. : " << to_iz << " > " << parent->getNZ() << endl;
        isValid = false;
    }

    if(!isValid) MPIw::Environment::exitWithFinalize(1);
}

void Grid::setField(Field* f){ field = f; }
Field* Grid::getField(void){ return field; }

//! 粒子の位置から電荷を空間電荷にする
void Grid::updateRho(const Environment* env) {
    tdArray& rho = field->getRho();

    ParticleType* ptype = env->ptype;
    for(int id = 0; id < env->num_of_particle_types; ++id){
        int pnum = ptype[id].getTotalNumber();

        for(int i = 0; i < pnum; ++i){
            double x = particles[id][i].getX();
            double y = particles[id][i].getY();
            double z = particles[id][i].getZ();

            int gx_lower = floor(x);
            double delta_gx = x - gx_lower;

            int gy_lower = floor(y);
            double delta_gy = y - gy_lower;

            int gz_lower = floor(z);
            double delta_gz = z - gz_lower;

            // glue cell分を考慮
            gx_lower += 1; gy_lower += 1; gz_lower += 1;

            double q = ptype[id].getCharge();

#ifdef DEBUG
            if(gx_lower + 1 >= env->cell_x + 2 || gy_lower + 1 >= env->cell_y + 2 || gz_lower + 1 >= env->cell_z + 2) {
                cout << env->rankStr() << format("[Particle]: %5f %5f %5f") % x % y % z << endl;
                cout << env->rankStr() << format("[Particle]: int + 1: %d %d %d") % (gx_lower+1) % (gy_lower+1) % (gz_lower+1) << endl;
            }
#endif

            rho[gx_lower    ][gy_lower    ][gz_lower    ] += (1.0 - delta_gx) * (1.0 - delta_gy) * (1.0 - delta_gz) * q;
            rho[gx_lower + 1][gy_lower    ][gz_lower    ] += delta_gx * (1.0 - delta_gy) * (1.0 - delta_gz) * q;
            rho[gx_lower    ][gy_lower + 1][gz_lower    ] += (1.0 - delta_gx) * delta_gy * (1.0 - delta_gz) * q;
            rho[gx_lower + 1][gy_lower + 1][gz_lower    ] += delta_gx * delta_gy * (1.0 - delta_gz) * q;

            rho[gx_lower    ][gy_lower    ][gz_lower + 1] += (1.0 - delta_gx) * (1.0 - delta_gy) * delta_gz * q;
            rho[gx_lower + 1][gy_lower    ][gz_lower + 1] += delta_gx * (1.0 - delta_gy) * delta_gz * q;
            rho[gx_lower    ][gy_lower + 1][gz_lower + 1] += (1.0 - delta_gx) * delta_gy * delta_gz * q;
            rho[gx_lower + 1][gy_lower + 1][gz_lower + 1] += delta_gx * delta_gy * delta_gz * q;
        }
    }

    // clear values on glue cell
    Utils::clearBoundaryValues(rho, env->cell_x + 2, env->cell_y + 2, env->cell_z + 2);
}

float** Grid::getMeshNodes(int dim) {
    // the array of coordinate arrays
    // @note: メモリリーク防止のため必ずdeleteする
    float** coordinates = new float*[dim];
    coordinates[0] = new float[nx];
    for(int i = 0; i < nx; ++i) {
	coordinates[0][i] = base_x + dx * i;
    }
    coordinates[1] = new float[ny];
    for(int i = 0; i < ny; ++i) {
	coordinates[1][i] = base_y + dx * i;
    }
    coordinates[2] = new float[nz];
    for(int i = 0; i < nz; ++i) {
	coordinates[2][i] = base_z + dx * i;
    }
    return coordinates;
}

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

int Grid::getMaxChildLevel() {
    return this->getMaxLevel() - level;
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

void Grid::putQuadMesh(DBfile* file, char* coordnames[3], char* varnames[1], int rankInGroup, DBoptlist* optListMesh, DBoptlist* optListVar){
    const int dim = 3;

    // the array of coordinate arrays
    float** coordinates = this->getMeshNodes(dim);

    //! @TODO: varnamesによってtdArrayを変える
    double* tdArray = Utils::getTrueCells(this->getField()->getPhi());

    // dimensions
    // glue cellも出力
    int dimensions[3];
    dimensions[0] = nx;
    dimensions[1] = ny;
    dimensions[2] = nz;
    double* vars[] = {tdArray};

    std::string meshname = (format("/block%04d/mesh%04d%04d") % rankInGroup % MPIw::Environment::rank % id).str();
    std::string varname = (format("/block%04d/%s%04d%04d") % rankInGroup % varnames[0] % MPIw::Environment::rank % id).str();

    char* m = const_cast<char*>(meshname.c_str());
    char* v = const_cast<char*>(varname.c_str());

    DBPutQuadmesh(file, m, coordnames, coordinates, dimensions, dim, DB_FLOAT, DB_COLLINEAR, optListMesh);
    // DBPutQuadvar(file, v, m, 1, varnames, vars, dimensions, dim, NULL, 0, DB_DOUBLE, DB_NODECENT, NULL);
    DBPutQuadvar1(file, v, m, tdArray, dimensions, dim, NULL, 0, DB_DOUBLE, DB_NODECENT, optListVar);

    for(int i = 0; i < children.size(); ++i) {
        children[i]->putQuadMesh(file, coordnames, varnames, rankInGroup, optListMesh, optListVar);
    }

    delete [] vars[0];
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
    for(int i = 0; i < g->getLevel(); ++i) tab += "  ";

    if(g->getLevel() > 0) ost << tab << "--- child [" << childnum << "] ---" << endl;
    ost << tab << "id: " << g->getID() << endl;
    ost << tab << "level: " << g->getLevel() << endl;
    ost << tab << "dx: " << format("%10.5e") % g->getDX() << "m" << endl;
    ost << tab << "nx, ny, nz: " << format("%1%x%2%x%3%") % g->getNX() % g->getNY() % g->getNZ() << " grids [total]" << endl;
    ost << tab << "nx,ny,nz(+): " << format("%1%x%2%x%3%") % (g->getNX() + 2) % (g->getNY() + 2) % (g->getNZ() + 2) << " grids [with glue cells]" << endl;
    ost << tab << "parent from: " << format("%1%,%2%,%3%") % g->getFromIX() % g->getFromIY() % g->getFromIZ() << endl;
    ost << tab << "parent to  : " << format("%1%,%2%,%3%") % g->getToIX() % g->getToIY() % g->getToIZ() << endl;
    ost << tab << "base positions: " << format("%1%,%2%,%3%") % g->getBaseX() % g->getBaseY() % g->getBaseZ() << endl;
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

