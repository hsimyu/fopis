#include <tdpic.h>
#include <random>

Grid::Grid(const Environment* env, const ParticleType* ptype){
    //! - コンストラクタにEnvironmentクラスが渡された場合、
    //! レベル0のGridを作成します.
    level = 0;
    base_x = 0.0;
    base_y = 0.0;
    base_z = 0.0;

    // 粒子数は自動算出する
    const int particle_number = 40;
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

    // std::vector< Particle[] > particle_arrays;
    for(int id = 0; id < env->particle_types; ++id){
        //! 各粒子種の配列を保持するためのvector
        //! ここで一度に確保しておけば連続になる？
        // std::unique_ptr<Particle[]> local_particles(new Particle[ptype[id]->getTotalNumber()]);
        std::unique_ptr<Particle[]> local_particles(new Particle[particle_number]);

        //! - 粒子はレベル0グリッドにのみ所属します
        // particles[id] = new Particle[particle_number]();

        for(int i = 0; i < particle_number; ++i){
            local_particles[i].setPosition(dist_x(mt_x), dist_y(mt_y), dist_z(mt_z));
        }

        //! std::moveするとcopy-constructorが呼ばれない
        particles.push_back(std::move(local_particles));
    }
}

Grid::~Grid(){
    //! delete all particles
    //! vector内のparticleはunique_ptrなので自動削除される
    particles.erase(particles.begin(), particles.end());
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
