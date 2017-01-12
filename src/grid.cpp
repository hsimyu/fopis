#include <tdpic.h>
#include <random>

// accessors
void Grid::setBaseIX(int _ix){ base_ix = _ix; }
void Grid::setBaseIY(int _iy){ base_iy = _iy; }
void Grid::setBaseIZ(int _iz){ base_iz = _iz; }
int  Grid::getBaseIX(void) const { return base_ix; }
int  Grid::getBaseIY(void) const { return base_iy; }
int  Grid::getBaseIZ(void) const { return base_iz; }

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

void Grid::makeChild(const int _base_ix, const int _base_iy, const int _base_iz, const int _nx, const int _ny, const int _nz) {
    Grid* child = new Grid(this, _base_ix, _base_iy, _base_iz, _nx, _ny, _nz);

    this->addChild(child);
}

void Grid::addChild(Grid* child) { children.push_back(child); }
std::vector<Grid*>& Grid::getChildren(void) {
    // 参照にしないと新しいポインタが生まれてしまう？
    return children;
}

unsigned int Grid::getChildrenLength(void) const {
    return children.size();
}

// root grid constructor
Grid::Grid(const Environment* env){
    //! - コンストラクタにEnvironmentクラスが渡された場合、
    //! レベル0のGridを作成します.
    level = 0;

    nx = env->cell_x;
    ny = env->cell_y;
    nz = env->cell_z;
    dx = env->dx;

    //! @{
    //! Root Gridの場合の親グリッドは、計算空間を全て統合した空間として、
    //! その上にプロセス分割されたグリッドが乗っていると考える
    base_ix = env->xrank * env->cell_x;
    base_iy = env->yrank * env->cell_y;
    base_iz = env->zrank * env->cell_z;
    //! @}

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
Grid::Grid(Grid* g, const int _base_ix, const int _base_iy, const int _base_iz, const int _nx, const int _ny, const int _nz){
    const double refineRatio = 2.0;

    parent = g;
    level = g->getLevel() + 1;

    // patchの大きさを指定
    nx = _nx;
    ny = _ny;
    nz = _nz;

    // refineRatioは2で固定
    dx = g->getDX() / refineRatio;

    //! @{
    //! 子グリッドの場合, base_ix変数は純粋に親グリッドの何番目に乗っているかを表す
    base_ix = _base_ix;
    base_iy = _base_ix;
    base_iz = _base_ix;
    //! @}

    checkGridValidness();
}

void Grid::checkGridValidness() {
    const double refineRatio = 2.0;

    bool isValid = true;

    if( nx % 2 == 0 ) {
        std::cerr << "[ERROR] x-extent is not odd number. : " << nx << endl;
        isValid = false;
    }

    if( ny % 2 == 0 ) {
        std::cerr << "[ERROR] y-extent is not odd number. : " << ny << endl;
        isValid = false;
    }

    if( nz % 2 == 0 ) {
        std::cerr << "[ERROR] z-extent is not odd number. : " << nz << endl;
        isValid = false;
    }

    // x extent
    if( (base_ix + (nx - 1)/ refineRatio) > parent->getNX() ){
        std::cerr << "[ERROR] A child patch's x-extent exceeds the parent's extent. : " << (base_ix + (nx - 1)/ refineRatio) << " > " << parent->getNX() << endl;
        isValid = false;
    }

    // y extent
    if( (base_iy + (ny - 1)/ refineRatio) > parent->getNY() ){
        std::cerr << "[ERROR] A child patch's y-extent exceeds the parent's extent. : " << (base_iy + (ny - 1)/ refineRatio) << " > " << parent->getNY() << endl;
        isValid = false;
    }

    // z extent
    if( (base_iz + (nz - 1)/ refineRatio) > parent->getNZ() ){
        std::cerr << "[ERROR] A child patch's z-extent exceeds the parent's extent. : " << (base_iz + (nz - 1)/ refineRatio) << " > " << parent->getNZ() << endl;
        isValid = false;
    }

    if(!isValid) MPI::Environment::exitWithFinalize(1);
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
	coordinates[0][i] = dx * (base_ix + i);
    }
    coordinates[1] = new float[ny];
    for(int i = 0; i < ny; ++i) {
	coordinates[1][i] = dx * (base_iy + i);
    }
    coordinates[2] = new float[nz];
    for(int i = 0; i < nz; ++i) {
	coordinates[2][i] = dx * (base_iz + i);
    }
    return coordinates;
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
    ost << tab << "level: " << g->getLevel() << endl;
    ost << tab << "dx: " << format("%10.5e") % g->getDX() << "m" << endl;
    ost << tab << "nx, ny, nz: " << format("%1%x%2%x%3%") % g->getNX() % g->getNY() % g->getNZ() << " grids [total]" << endl;
    ost << tab << "base grids: " << format("%1%,%2%,%3%") % g->getBaseIX() % g->getBaseIY() % g->getBaseIZ() << endl;
    ost << tab << "numOfChild: " << g->getChildrenLength() << endl;

    if(g->getChildrenLength() > 0) {
        std::vector<Grid*>& children = g->getChildren();
        for(int i = 0; i < children.size(); ++i) {
            printGridInfo(ost, children[i], i);
        }
    }
}

std::ostream& operator<<(std::ostream& ost, Grid* g){
    ost << "[Grid]" << std::endl;
    printGridInfo(ost, g, 0);
    return ost;
}

