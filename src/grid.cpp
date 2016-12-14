#include <tdpic.h>
#include <random>

Grid::Grid(const Environment* env){
    //! - コンストラクタにEnvironmentクラスが渡された場合、
    //! レベル0のGridを作成します.
    level = 0;
    base_x = 0.0;
    base_y = 0.0;
    base_z = 0.0;

    // 粒子数は自動算出する
    const float max_x = env->nx * env->dx, max_y = env->ny * env->dx, max_z = env->nz * env->dx;

    // std::random_device rnd;
    const int random_src_x = 10684930;
    const int random_src_y = 99881;
    const int random_src_z = 861200045;
    std::mt19937 mt_x(random_src_x);
    std::mt19937 mt_y(random_src_y);
    std::mt19937 mt_z(random_src_z);

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

        //! - 粒子はレベル0グリッドにのみ所属します
        for(int i = 0; i < pnum; ++i){
            particles[id][i].setPosition(dist_x(mt_x), dist_y(mt_y), dist_z(mt_z));
        }
    }
}

Grid::~Grid(){
    //! delete all particles
    //! vector内のparticleは自動でデストラクタが呼ばれる
    particles.erase(particles.begin(), particles.end());

    // reserveしてあった分を削除する
    particles.shrink_to_fit();
}

void Grid::setBaseX(int _x){ base_x = _x; }
void Grid::setBaseY(int _y){ base_y = _y; }
void Grid::setBaseZ(int _z){ base_z = _z; }
int Grid::getBaseX(void){ return base_x; }
int Grid::getBaseY(void){ return base_y; }
int Grid::getBaseZ(void){ return base_z; }

void Grid::setLevel(int l){ level = l; }
int Grid::getLevel(void){ return level; }

void Grid::setParent(Grid* g){ parent = g; }
Grid* Grid::getParent(void){ return parent; }
