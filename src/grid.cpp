#include <tdpic.h>
#include <random>
using boost::format;

Grid::Grid(const Environment* env){
    //! - コンストラクタにEnvironmentクラスが渡された場合、
    //! レベル0のGridを作成します.
    level = 0;
    base_x = 0.0;
    base_y = 0.0;
    base_z = 0.0;

    nx = env->cell_x;
    ny = env->cell_y;
    nz = env->cell_z;
    dx = env->dx;

    // 粒子位置の上限を設定
    // 下限はbase_xになる
    const float max_x = env->cell_x;
    const float max_y = env->cell_y;
    const float max_z = env->cell_z;

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

    std::uniform_real_distribution<> dist_x(base_x, max_x);
    std::uniform_real_distribution<> dist_y(base_y, max_y);
    std::uniform_real_distribution<> dist_z(base_z, max_z);

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
    threeDArray& rho = field->getRho();

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
