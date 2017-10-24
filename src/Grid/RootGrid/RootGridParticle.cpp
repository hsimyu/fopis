#include "grid.hpp"

void RootGrid::injectParticlesFromBoundary(void) {
    static std::vector< std::vector<double> > residual;
    static bool isFirstCall = true;

    //! staticな残余変数の初期化
    if(isFirstCall) {
        residual.resize(Environment::getNumOfAmbientParticles());

        for(int i = 0; i < residual.size(); ++i) {
            residual[i].resize(6);
            for(int j = 0; j < 6; ++j) {
                residual[i][j] = 0.0;
            }
        }

        isFirstCall = false;
    }

    //! - 粒子位置の上限を設定
    double max_x = static_cast<double>(Environment::cell_x);
    double max_y = static_cast<double>(Environment::cell_y);
    double max_z = static_cast<double>(Environment::cell_z);

    //! - 上側境界にいる場合は外側にはみ出した粒子を生成しないようにする
    if(Environment::isBoundary(AXIS::x, AXIS_SIDE::up)) max_x -= 1.0;
    if(Environment::isBoundary(AXIS::y, AXIS_SIDE::up)) max_y -= 1.0;
    if(Environment::isBoundary(AXIS::z, AXIS_SIDE::up)) max_z -= 1.0;

	auto ambient_ptype_ptr_list = Environment::getAmbientParticleTypes();
    for(int itr = 0; itr < Environment::getNumOfAmbientParticles(); ++itr) {
        auto ambient_particle_ptr = ambient_ptype_ptr_list[itr];

        std::vector<double> flux = ambient_particle_ptr->calcFlux();
        const auto pid = ambient_particle_ptr->getId();

        if(Environment::isBoundary(AXIS::x, AXIS_SIDE::low)) {
            const double injecting_area_size = dx * dx * (this->getZNodeSize() - 1) * (this->getYNodeSize() - 1);

            const int index = 0;
            const double inject_num_double = dt * injecting_area_size * flux[index];
            const int inject_num = static_cast<int>(floor(inject_num_double + residual[itr][index]));
            residual[itr][index] += inject_num_double - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Velocity vel = ambient_particle_ptr->generateNewInjectionVelocity(AXIS::x, AXIS_SIDE::low);
                Particle p = ambient_particle_ptr->generateNewParticle(0.0, vel.vx * dt, 0.0, max_y, 0.0, max_z, vel);
                particles[pid].push_back( std::move(p) );
            }
        }

        if(Environment::isBoundary(AXIS::x, AXIS_SIDE::up)) {
            const double injecting_area_size = dx * dx * (this->getZNodeSize() - 1) * (this->getYNodeSize() - 1);

            const int index = 1;
            const double inject_num_double = dt * injecting_area_size * flux[index];
            const int inject_num = static_cast<int>(floor(inject_num_double + residual[itr][index]));
            residual[itr][index] += inject_num_double - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Velocity vel = ambient_particle_ptr->generateNewInjectionVelocity(AXIS::x, AXIS_SIDE::up);

                //! 負方向速度をxの最大値から引いた点までがありうる範囲
                Particle p = ambient_particle_ptr->generateNewParticle(max_x + vel.vx * dt, max_x, 0.0, max_y, 0.0, max_z, vel);
                particles[pid].push_back( std::move(p) );
            }
        }

        if(Environment::isBoundary(AXIS::y, AXIS_SIDE::low)) {
            const double injecting_area_size = dx * dx * (this->getZNodeSize() - 1) * (this->getXNodeSize() - 1);

            const int index = 2;
            const double inject_num_double = dt * injecting_area_size * flux[index];
            const int inject_num = static_cast<int>(floor(inject_num_double + residual[itr][index]));
            residual[itr][index] += inject_num_double - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Velocity vel = ambient_particle_ptr->generateNewInjectionVelocity(AXIS::y, AXIS_SIDE::low);
                Particle p = ambient_particle_ptr->generateNewParticle(0.0, max_x, 0.0, vel.vy * dt, 0.0, max_z, vel);
                particles[pid].push_back( std::move(p) );
            }
        }

        if(Environment::isBoundary(AXIS::y, AXIS_SIDE::up)) {
            const double injecting_area_size = dx * dx * (this->getZNodeSize() - 1) * (this->getXNodeSize() - 1);

            const int index = 3;
            const double inject_num_double = dt * injecting_area_size * flux[index];
            const int inject_num = static_cast<int>(floor(inject_num_double + residual[itr][index]));
            residual[itr][index] += inject_num_double - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Velocity vel = ambient_particle_ptr->generateNewInjectionVelocity(AXIS::y, AXIS_SIDE::up);
                Particle p = ambient_particle_ptr->generateNewParticle(0.0, max_x, max_y + vel.vy * dt, max_y, 0.0, max_z, vel);
                particles[pid].push_back( std::move(p) );
            }
        }

        if(Environment::isBoundary(AXIS::z, AXIS_SIDE::low)) {
            const double injecting_area_size = dx * dx * (this->getXNodeSize() - 1) * (this->getYNodeSize() - 1);

            const int index = 4;
            const double inject_num_double = dt * injecting_area_size * flux[index];
            const int inject_num = static_cast<int>(floor(inject_num_double + residual[itr][index]));
            residual[itr][index] += inject_num_double - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Velocity vel = ambient_particle_ptr->generateNewInjectionVelocity(AXIS::z, AXIS_SIDE::low);
                Particle p = ambient_particle_ptr->generateNewParticle(0.0, max_x, 0.0, max_y, 0.0, vel.vz * dt, vel);
                particles[pid].push_back( std::move(p) );
            }
        }

        if(Environment::isBoundary(AXIS::z, AXIS_SIDE::up)) {
            const double injecting_area_size = dx * dx * (this->getXNodeSize() - 1) * (this->getYNodeSize() - 1);

            const int index = 5;
            const double inject_num_double = dt * injecting_area_size * flux[index];
            const int inject_num = static_cast<int>(floor(inject_num_double + residual[itr][index]));
            residual[itr][index] += inject_num_double - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Velocity vel = ambient_particle_ptr->generateNewInjectionVelocity(AXIS::z, AXIS_SIDE::up);
                Particle p = ambient_particle_ptr->generateNewParticle(0.0, max_x, 0.0, max_y, max_z + vel.vz * dt, max_z, vel);
                particles[pid].push_back( std::move(p) );
            }
        }
    }
}
