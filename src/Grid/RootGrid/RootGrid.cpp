#include "grid.hpp"
#include "normalizer.hpp"
#include "field.hpp"
#include "utils.hpp"
#include "dataio.hpp"

RootGrid::RootGrid() : Grid() {
    level = 0;

    //! UniqueなIDをセット
    constexpr int minimum_id_offset = 10;
    id = minimum_id_offset * MPIw::Environment::rank + this->getNextID();

    nx = Environment::cell_x;
    ny = Environment::cell_y;
    nz = Environment::cell_z;
    dx = Normalizer::normalizeLength(Environment::dx);
    dt = Normalizer::normalizeTime(Environment::dt);

    //! @{
    //! Root Gridの場合の親グリッドは、計算空間を全て統合した空間として、
    //! その上にプロセス分割されたグリッドが乗っていると考える
    from_ix = Environment::getAssignedXBegin();
    from_iy = Environment::getAssignedYBegin();
    from_iz = Environment::getAssignedZBegin();
    to_ix = Environment::getAssignedXEnd();
    to_iy = Environment::getAssignedYEnd();
    to_iz = Environment::getAssignedZEnd();
    //! @note: base_x, base_y, base_zは正規化された長さ
    base_x = dx * static_cast<double>(from_ix);
    base_y = dx * static_cast<double>(from_iy);
    base_z = dx * static_cast<double>(from_iz);
    //! @}

    //! ChildMap初期化
    this->initializeChildMap();

    // Field初期化
    this->initializeField();

    // 物体初期化
    this->initializeObjects();
    this->initializeObjectsCmatrix();

    //! - 粒子位置の上限を設定
    double max_x = static_cast<double>(Environment::cell_x);
    double max_y = static_cast<double>(Environment::cell_y);
    double max_z = static_cast<double>(Environment::cell_z);

    //! - 上側境界にいる場合は外側にはみ出した粒子を生成しないようにする
    if (Environment::isBoundary(AXIS::x, AXIS_SIDE::up)) max_x -= 1.0;
    if (Environment::isBoundary(AXIS::y, AXIS_SIDE::up)) max_y -= 1.0;
    if (Environment::isBoundary(AXIS::z, AXIS_SIDE::up)) max_z -= 1.0;

    //! - particlesは空のstd::vector< std::vector<Particle> >として宣言されている
    //! - particle types 分だけresize
    particles.resize(Environment::num_of_particle_types);

    for(int id = 0; id < Environment::num_of_particle_types; ++id){
        //! 各粒子分のメモリをreserveしておく
        particles[id].reserve(Environment::getParticleType(id)->getTotalNumber() * 2);

        //! 初期化時は背景粒子のみ生成
        if (Environment::getParticleType(id)->getType() == "ambient") {
            auto ambient_particle_ptr = Environment::getAmbientParticleType(id);
            int pnum = ambient_particle_ptr->getTotalNumber();

            for(int i = 0; i < pnum; ++i){
                Particle p = ambient_particle_ptr->generateNewParticle(0.0, max_x, 0.0, max_y, 0.0, max_z);

                //! 物体がある場合は生成時にチェックする
                for(const auto& obj : objects) {
                    if (obj.isDefined()) obj.removeInnerParticle(p);
                }

                if (p.isValid) particles[id].push_back( std::move(p) );
            }
        }
    }
}

void RootGrid::mainLoopES() {
    auto time_counter = Utils::TimeCounter::getInstance();

    time_counter->begin("resetObjects");
    this->resetObjects();

    // -- timing: 0 --

    //! 子グリッド上のループを1回分先に呼び出し、
    //! 各関数においてももう一度呼び出すことで時間感覚を合わせる
    for (auto& child : children) {
        // -- timing: 0.5 dt --
        child->mainLoop();
    }

    // -- timing: dt --
    // 速度更新
    time_counter->switchTo("updateParticleVelocity");
    this->updateParticleVelocity();

    // 粒子放出
    time_counter->switchTo("emitParticlesFromObjects");
    this->emitParticlesFromObjects();

    // 粒子位置更新
    time_counter->switchTo("updateParticlePosition");
    this->updateParticlePosition();

    // 粒子注入
    time_counter->switchTo("injectParticleFromBoundary");
    this->injectParticlesFromBoundary();

    // 密度更新
    time_counter->switchTo("updateDensity");
    this->updateDensity();

    // 静電計算の場合
    // 新しい位置に対応する電荷密度算出
    time_counter->switchTo("updateRho");
    this->updateRho();

    // Poisson を解く
    time_counter->switchTo("solvePoisson");
    this->solvePoisson();

    // 電場更新
    time_counter->switchTo("updateEfield");
    this->updateEfield();

    // 物体情報更新
    time_counter->switchTo("updateObjects");
    this->updateObjects();

    //! 不要な粒子を削除
    constexpr unsigned int remove_invalid_particle_width = 50;
    if (Environment::timestep % remove_invalid_particle_width == 0) {
        time_counter->switchTo("removeInvalidParticles");
        this->removeInvalidParticles();
    }

    time_counter->end();
}

void RootGrid::mainLoopEM() {
    auto time_counter = Utils::TimeCounter::getInstance();

    time_counter->begin("resetObjects");
    this->resetObjects();

    // -- timing: t + 0.5 dt --
    // 速度更新
    time_counter->switchTo("updateParticleVelocity");
    this->updateParticleVelocity();

    // 粒子放出
    time_counter->switchTo("emitParticlesFromObjects");
    this->emitParticlesFromObjects();

    // 位置更新
    time_counter->switchTo("updateParticlePosition");
    this->updateParticlePosition(); // jx, jy, jz もここで update される

    // 粒子注入
    time_counter->switchTo("injectParticleFromBoundary");
    this->injectParticlesFromBoundary();

    // 新しい位置に対応する電荷密度算出
    time_counter->switchTo("updateRho");
    this->updateRho();

    time_counter->switchTo("updateEfieldFDTD");
    this->updateEfieldFDTD();

    // time_counter->switchTo("solvePoissonCorrection");
    // this->solvePoissonCorrection();

    // 磁場更新
    time_counter->switchTo("updateBfield");
    this->updateBfield();

    //! 不要な粒子を削除
    constexpr unsigned int remove_invalid_particle_width = 50;
    if (Environment::timestep % remove_invalid_particle_width == 0) {
        time_counter->switchTo("removeInvalidParticles");
        this->removeInvalidParticles();
    }

    time_counter->end();
}

int RootGrid::getXNodeSize(void) const {
    //! 周期境界の場合は上側境界と下側境界の間の空間も有効な空間となるので、
    //! 上側のノードを1つ増やす
    return (level == 0 && Environment::isNotBoundary(AXIS::x, AXIS_SIDE::up)) ? nx + 1 : nx;
}

int RootGrid::getYNodeSize(void) const {
    return (level == 0 && Environment::isNotBoundary(AXIS::y, AXIS_SIDE::up)) ? ny + 1 : ny;
}

int RootGrid::getZNodeSize(void) const {
    return (level == 0 && Environment::isNotBoundary(AXIS::z, AXIS_SIDE::up)) ? nz + 1 : nz;
}

inline void RootGrid::decrementSumOfChild() { --sumTotalNumOfChildGrids; }
inline void RootGrid::incrementSumOfChild() { ++sumTotalNumOfChildGrids; }