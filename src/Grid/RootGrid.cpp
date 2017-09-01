#include "grid.hpp"
#include "normalizer.hpp"
#include "field.hpp"
#include "utils.hpp"
#include "dataio.hpp"

RootGrid::RootGrid() : Grid() {
    level = 0;

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

    // Field初期化
    this->initializeField();

    // 物体初期化
    this->initializeObject();
    this->initializeObjectsCmatrix();

    //! - 粒子位置の上限を設定
    double max_x = static_cast<double>(Environment::cell_x);
    double max_y = static_cast<double>(Environment::cell_y);
    double max_z = static_cast<double>(Environment::cell_z);

    //! - 上側境界にいる場合は外側にはみ出した粒子を生成しないようにする
    if(!Environment::isNotBoundary(AXIS::x, AXIS_SIDE::up)) max_x -= 1.0;
    if(!Environment::isNotBoundary(AXIS::y, AXIS_SIDE::up)) max_y -= 1.0;
    if(!Environment::isNotBoundary(AXIS::z, AXIS_SIDE::up)) max_z -= 1.0;

    //! - particlesは空のstd::vector< std::vector<Particle> >として宣言されている
    //! - particle types 分だけresize
    particles.resize(Environment::num_of_particle_types);

    for(int id = 0; id < Environment::num_of_particle_types; ++id){
        //! 各粒子分のメモリをreserveしておく
        particles[id].reserve(Environment::max_particle_num);

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

void RootGrid::solvePoisson(void) {
    const int DEFAULT_ITERATION_LOOP = 500;
    field->solvePoisson(DEFAULT_ITERATION_LOOP, dx);
}

void RootGrid::updateEfield(void) {
    field->updateEfield(dx);
}

void RootGrid::updateEfieldFDTD(void) {
    field->updateEfieldFDTD(dx, dt);
}

void RootGrid::updateBfield(void) {
    field->updateBfield(dx, nx, ny, nz, dt);
}

void RootGrid::initializeObject(void) {
    if (Environment::isRootNode) cout << "-- Defining Objects -- " << endl;

    for (const auto& object_info : Environment::objects_info) {
        std::string obj_name = object_info.name;
        //! 物体関連の設定を関連付けされた obj 形式ファイルから読み込む
        ObjectDataFromFile object_data = ObjectUtils::getObjectNodesFromObjFile(object_info.file_name);

        unsigned int num_cmat = object_data.nodes.size();
        const auto& node_array = object_data.nodes;

        //! innerと判定されたやつだけ渡す
        ObjectNodes inner_node_array;
        bool is_object_in_this_node = false;

        for(const auto& node_pair : node_array) {
            const auto cmat_itr = node_pair.first;
            const auto& node_pos = node_pair.second;

            const auto i = node_pos[0];
            const auto j = node_pos[1];
            const auto k = node_pos[2];

            if (isInnerNode(i, j, k)) {
                is_object_in_this_node = true;
                inner_node_array[cmat_itr] = Environment::getRelativePositionOnRootWithGlue(i, j, k);
            }
        }

        const auto& cell_array = object_data.cells;
        ObjectCells inner_cell_array;
        for(const auto& cell_pos : cell_array) {
            if (isInnerCellWithGlue(cell_pos[0], cell_pos[1], cell_pos[2])) {
                auto rel_pos = Environment::getRelativePositionOnRootWithGlue(cell_pos[0], cell_pos[1], cell_pos[2]);
                inner_cell_array.push_back( {{rel_pos[0], rel_pos[1], rel_pos[2], cell_pos[3]}} );
            }
        }

        //! 物体定義点がゼロでも Spacecraft オブジェクトだけは作成しておいた方がよい
        //! emplace_back で Spacecraft object を直接構築
        objects.emplace_back(nx, ny, nz, num_cmat, object_info, inner_node_array, inner_cell_array);

        //! Comm作成 (物体が入っていないならnullになる)
        MPIw::Environment::makeNewComm(obj_name, is_object_in_this_node);
        if (MPIw::Environment::isRootNode(obj_name)) {
            cout << Environment::rankStr() << "is set to Root Node for " << obj_name << "." << endl;
            cout << objects[ objects.size() - 1 ] << endl;
        }
    }
}

void RootGrid::initializeObjectsCmatrix(void) {
    if (Environment::isRootNode) cout << "-- Initializing Objects Capacity Matrix --" << endl;
    RhoArray& rho = field->getRho();
    tdArray& phi = field->getPhi();

    for(auto& obj : objects) {
        // rhoを初期化
        Utils::initializeRhoArray(rho);
        Utils::initialize3DArray(phi);

        //! データを読み込み
        if (Environment::useExistingCapacityMatrix) {
            const auto is_valid_load = IO::loadCmatrixData(obj);
            if (is_valid_load) continue;
        }

        //! データを使わずに初期化
        const auto num_cmat = obj.getCmatSize();

        std::unique_ptr<Utils::ProgressManager> pm;
        if (MPIw::Environment::isRootNode(obj.getName())) {
            pm = std::make_unique<Utils::ProgressManager>(num_cmat, "cmat_solve");
        }

        for(unsigned int cmat_col_itr = 0; cmat_col_itr < num_cmat; ++cmat_col_itr ) {
            if (MPIw::Environment::isRootNode(obj.getName())) pm->update(cmat_col_itr);

            //! 該当する頂点に単位電荷を付与
            if (obj.isMyCmat(cmat_col_itr)) {
                const auto& cmat_pos = obj.getCmatPos(cmat_col_itr);
                rho[0][cmat_pos.i][cmat_pos.j][cmat_pos.k] = 1.0;
            }

            solvePoisson();

            for(unsigned int cmat_row_itr = 0; cmat_row_itr < num_cmat; ++cmat_row_itr ) {
                double value = 0.0;
                if (obj.isMyCmat(cmat_row_itr)) {
                    //! phiの値がB_{ij}の値になっている
                    const auto& target_pos = obj.getCmatPos(cmat_row_itr);
                    value = phi[target_pos.i][target_pos.j][target_pos.k];
                }

                if (obj.isDefined()) {
                    //! bcastの代わりにsumしてしまう
                    value = MPIw::Environment::Comms[obj.getName()].sum(value);
                    obj.setCmatValue(cmat_col_itr, cmat_row_itr, value);
                }
            }

            //! 付与した単位電荷を消去する
            if (obj.isMyCmat(cmat_col_itr)) {
                const auto& cmat_pos = obj.getCmatPos(cmat_col_itr);
                rho[0][cmat_pos.i][cmat_pos.j][cmat_pos.k] = 0.0;
            }
        }

        //! 物体が有効でないなら解く必要なし
        if (obj.isDefined()) obj.makeCmatrixInvert();

        if (MPIw::Environment::isRootNode(obj.getName())) {
            IO::writeCmatrixData(obj);
        }
    }
}

void RootGrid::resetObjects() {
    //! 物体上の一時的な情報を初期化
    for(auto& obj : objects) {
        if(obj.isDefined()) {
            obj.resetCurrent();
        }
    }
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
    if(!Environment::isNotBoundary(AXIS::x, AXIS_SIDE::up)) max_x -= 1.0;
    if(!Environment::isNotBoundary(AXIS::y, AXIS_SIDE::up)) max_y -= 1.0;
    if(!Environment::isNotBoundary(AXIS::z, AXIS_SIDE::up)) max_z -= 1.0;

	auto ambient_ptype_ptr_list = Environment::getAmbientParticleTypes();
    for(int itr = 0; itr < Environment::getNumOfAmbientParticles(); ++itr) {
        auto ambient_particle_ptr = ambient_ptype_ptr_list[itr];

        std::vector<double> flux = ambient_particle_ptr->calcFlux(*this);
        const auto pid = ambient_particle_ptr->getId();

        if(!Environment::isNotBoundary(AXIS::x, AXIS_SIDE::low)) {
            const int index = 0;
            const int inject_num = floor(dt * flux[index] + residual[itr][index]);
            residual[itr][index] += dt * flux[index] - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Velocity vel = ambient_particle_ptr->generateNewVelocity();

                //! 流入方向速度に変換
                //! 実際はフラックスを積分して割合を求める必要がある?
                while (vel.vx <= 0.0) {
                    vel = ambient_particle_ptr->generateNewVelocity();
                }

                Particle p = ambient_particle_ptr->generateNewParticle(0.0, vel.vx * dt, 0.0, max_y, 0.0, max_z, vel);
                particles[pid].push_back( std::move(p) );
            }
        }

        if(!Environment::isNotBoundary(AXIS::x, AXIS_SIDE::up)) {
            const int index = 1;
            const int inject_num = floor(dt * flux[index] + residual[itr][index]);
            residual[itr][index] += dt * flux[index] - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Velocity vel = ambient_particle_ptr->generateNewVelocity();

                while (vel.vx >= 0.0) {
                    vel = ambient_particle_ptr->generateNewVelocity();
                }

                //! 負方向速度をxの最大値から引いた点までがありうる範囲
                Particle p = ambient_particle_ptr->generateNewParticle(max_x + vel.vx * dt, max_x, 0.0, max_y, 0.0, max_z, vel);
                particles[pid].push_back( std::move(p) );
            }
        }

        if(!Environment::isNotBoundary(AXIS::y, AXIS_SIDE::low)) {
            const int index = 2;
            const int inject_num = floor(dt * flux[index] + residual[itr][index]);
            residual[itr][index] += dt * flux[index] - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Velocity vel = ambient_particle_ptr->generateNewVelocity();

                while (vel.vy <= 0.0) {
                    vel = ambient_particle_ptr->generateNewVelocity();
                }

                Particle p = ambient_particle_ptr->generateNewParticle(0.0, max_x, 0.0, vel.vy * dt, 0.0, max_z, vel);
                particles[pid].push_back( std::move(p) );
            }
        }

        if(!Environment::isNotBoundary(AXIS::y, AXIS_SIDE::up)) {
            const int index = 3;
            const int inject_num = floor(dt * flux[index] + residual[itr][index]);
            residual[itr][index] += dt * flux[index] - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Velocity vel = ambient_particle_ptr->generateNewVelocity();

                while (vel.vy >= 0.0) {
                    vel = ambient_particle_ptr->generateNewVelocity();
                }

                Particle p = ambient_particle_ptr->generateNewParticle(0.0, max_x, max_y + vel.vy * dt, max_y, 0.0, max_z, vel);
                particles[pid].push_back( std::move(p) );
            }
        }

        if(!Environment::isNotBoundary(AXIS::z, AXIS_SIDE::low)) {
            const int index = 4;
            const int inject_num = floor(dt * flux[index] + residual[itr][index]);
            residual[itr][index] += dt * flux[index] - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Velocity vel = ambient_particle_ptr->generateNewVelocity();

                while (vel.vz <= 0.0) {
                    vel = ambient_particle_ptr->generateNewVelocity();
                }

                Particle p = ambient_particle_ptr->generateNewParticle(0.0, max_x, 0.0, max_y, 0.0, vel.vz * dt, vel);
                particles[pid].push_back( std::move(p) );
            }
        }

        if(!Environment::isNotBoundary(AXIS::z, AXIS_SIDE::up)) {
            const int index = 5;
            const int inject_num = floor(dt * flux[index] + residual[itr][index]);
            residual[itr][index] += dt * flux[index] - inject_num;

            for(int i = 0; i < inject_num; ++i) {
                Velocity vel = ambient_particle_ptr->generateNewVelocity();

                while (vel.vz >= 0.0) {
                    vel = ambient_particle_ptr->generateNewVelocity();
                }

                Particle p = ambient_particle_ptr->generateNewParticle(0.0, max_x, 0.0, max_y, max_z + vel.vz * dt, max_z, vel);
                particles[pid].push_back( std::move(p) );
            }
        }
    }
}

void RootGrid::emitParticlesFromObjects(void) {
    for(auto& obj : objects) {
        if (obj.isDefined() && obj.hasEmitParticles()) {
            obj.emitParticles(particles);
        }
    }
}

//! 粒子の位置から電荷を空間電荷にする
void RootGrid::updateRho() {
    RhoArray& rho = field->getRho();

#ifdef CHARGE_CONSERVATION
    // 電荷保存則をcheckするため、古いrhoを保持する
    auto old_rho = rho;
#endif

    //! rhoを初期化
    Utils::initializeRhoArray(rho);

    for(int pid = 0; pid < Environment::num_of_particle_types; ++pid){
        double q = Environment::getParticleType(pid)->getChargeOfSuperParticle();
        const auto rho_idx = pid + 1;

        for(auto& p : particles[pid]) {
            if(p.isValid) {
                for(auto& obj : objects) {
                    //! 物体中にいた場合には自動的に invalid になる
                    obj.distributeInnerParticleCharge(p);
                }

                //! もし物体内でなければ
                if (p.isValid) {
                    const auto pos = p.getPosition();
                    const int i = pos.i, j = pos.j, k = pos.k;

                    rho[rho_idx][i  ][j  ][k] += pos.dx2 * pos.dy2 * pos.dz2 * q;
                    rho[rho_idx][i+1][j  ][k] += pos.dx1 * pos.dy2 * pos.dz2 * q;
                    rho[rho_idx][i  ][j+1][k] += pos.dx2 * pos.dy1 * pos.dz2 * q;
                    rho[rho_idx][i+1][j+1][k] += pos.dx1 * pos.dy1 * pos.dz2 * q;
                    rho[rho_idx][i  ][j  ][k+1] += pos.dx2 * pos.dy2 * pos.dz1 * q;
                    rho[rho_idx][i+1][j  ][k+1] += pos.dx1 * pos.dy2 * pos.dz1 * q;
                    rho[rho_idx][i  ][j+1][k+1] += pos.dx2 * pos.dy1 * pos.dz1 * q;
                    rho[rho_idx][i+1][j+1][k+1] += pos.dx1 * pos.dy1 * pos.dz1 * q;
                }
            }
        }
    }

    for(auto& obj : objects) {
        if (obj.isDefined()) {
            //! 物体に配分された電荷を現在の rho に印加する
            obj.applyCharge(rho);
        }
    }

    for(int pid = 1; pid < Environment::num_of_particle_types + 1; ++pid) {
        for(int i = 1; i < nx + 2; ++i) {
            for(int j = 1; j < ny + 2; ++j) {
                for(int k = 1; k < nz + 2; ++k) {
                    rho[0][i][j][k] += rho[pid][i][j][k];
                }
            }
        }
    }

    //! rho を隣に送る
    for(int pid = 0; pid < Environment::num_of_particle_types + 1; ++pid) {
        MPIw::Environment::sendRecvField(rho[pid]);
    }

    //! 物体が存在する場合、電荷再配分が必要になる
    if (objects.size() > 0) {
        //! 一度 Poisson を解いて phi を更新
        solvePoisson();

        for(auto& obj : objects) {
            if (obj.isDefined()) obj.redistributeCharge(rho, field->getPhi());
        }

        MPIw::Environment::sendRecvField(rho[0]);
    }

#ifdef CHARGE_CONSERVATION
    if (Environment::solver_type == "EM") {
        field->checkChargeConservation(old_rho, 1.0, dx);
    }
#endif
}

inline void RootGrid::decrementSumOfChild() { --sumTotalNumOfChildGrids; }
inline void RootGrid::incrementSumOfChild() { ++sumTotalNumOfChildGrids; }
