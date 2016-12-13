#include <tdpic.h>

Grid::Grid(const Environment* env){
    //! - コンストラクタにEnvironmentクラスが渡された場合、
    //! レベル0のGridを作成します.
    level = 0;
    base_x = 0.0;
    base_y = 0.0;
    base_z = 0.0;

    /*
    // particles = new Particle[env->particle_types];
    particles.length = env->particle_types;

    //! - 粒子はレベル0グリッドにのみ所属します
    for(int id = 0; id < env->particle_types; ++id){
        particles[id] = new Particle[4]();

        particles[id][0]->setX(10.0);
    }*/
}

Grid::~Grid(){
    // delete [] children;
    // delete [] particles;
}
