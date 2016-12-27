#include <tdpic.h>
#include <random>

Grid::Grid(const Environment* env){
    //! - コンストラクタにEnvironmentクラスが渡された場合、
    //! レベル0のGridを作成します.
    level = 0;

    //! @{
    //! Root Gridの場合の親グリッドは、計算空間を全て統合した空間として、
    //! その上にプロセス分割されたグリッドが乗っていると考える
    base_x = static_cast<double>(env->xid * env->cell_x);
    base_y = static_cast<double>(env->yid * env->cell_y);
    base_z = static_cast<double>(env->zid * env->cell_z);
    //! @}

    nx = env->cell_x;
    ny = env->cell_y;
    nz = env->cell_z;
    dx = env->dx;

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

void Grid::setBaseX(int _x){ base_x = _x; }
void Grid::setBaseY(int _y){ base_y = _y; }
void Grid::setBaseZ(int _z){ base_z = _z; }
int Grid::getBaseX(void){ return base_x; }
int Grid::getBaseY(void){ return base_y; }
int Grid::getBaseZ(void){ return base_z; }
void Grid::setNX(int _x){ nx = _x; }
void Grid::setNY(int _y){ ny = _y; }
void Grid::setNZ(int _z){ nz = _z; }
int Grid::getNX(void){ return nx; }
int Grid::getNY(void){ return ny; }
int Grid::getNZ(void){ return nz; }

void Grid::setLevel(int l){ level = l; }
int Grid::getLevel(void){ return level; }

void Grid::setDX(double _dx){ dx = _dx; }
double Grid::getDX(void){ return dx; }

void Grid::setParent(Grid* g){ parent = g; }
Grid* Grid::getParent(void){ return parent; }

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
                cout << format("[P%d Particle]: %5f %5f %5f") % env->myid % x % y % z << endl;
                cout << format("[P%d Particle]: int + 1: %d %d %d") % env->myid % (gx_lower+1) % (gy_lower+1) % (gz_lower+1) << endl;
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
    // @note: メモリリークしそう
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


Grid::~Grid(){
    //! delete all particles
    //! vector内のparticleは自動でデストラクタが呼ばれる
    particles.erase(particles.begin(), particles.end());

    // reserveしてあった分を削除する
    particles.shrink_to_fit();

    delete field;
}
