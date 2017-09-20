#include "grid.hpp"
#include "normalizer.hpp"
#include "field.hpp"
#include "utils.hpp"
#include "dataio.hpp"

#define USE_BOOST
#include "simple_vtk.hpp"

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

    //! ChildMap初期化
    this->initializeChildMap();

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
    constexpr int PRE_LOOP_NUM = 10;
    constexpr int POST_LOOP_NUM = 250;

    if (this->getChildrenLength() > 0) {
        this->solvePoissonPSOR(PRE_LOOP_NUM);
        this->restrictPhiToChildrenBoundary();

        for(auto& child : children) {
            child->solvePoisson();
        }
    }

    this->solvePoissonPSOR(POST_LOOP_NUM);

    if (this->getChildrenLength() > 0) {
        this->correctChildrenPhi();
    }
}

void RootGrid::solvePoissonPSOR(const int loopnum) {
    auto& phi = field->getPhi();
    auto& poisson_error = field->getPoissonError();
    auto& rho = field->getRho();

    const double omega = 2.0/(1.0 + sin(M_PI/(phi.shape()[0] - 2))); // spectral radius
    const double rho_coeff = pow(dx, 2) / Normalizer::eps0;

    const int cx_with_glue = phi.shape()[0];
    const int cy_with_glue = phi.shape()[1];
    const int cz_with_glue = phi.shape()[2];

    constexpr double required_error = 1.0e-7;

    field->setBoundaryConditionPhi();

    const bool is_not_boundary[6] = {
        Environment::isNotBoundary(AXIS::x, AXIS_SIDE::low),
        Environment::isNotBoundary(AXIS::x, AXIS_SIDE::up),
        Environment::isNotBoundary(AXIS::y, AXIS_SIDE::low),
        Environment::isNotBoundary(AXIS::y, AXIS_SIDE::up),
        Environment::isNotBoundary(AXIS::z, AXIS_SIDE::low),
        Environment::isNotBoundary(AXIS::z, AXIS_SIDE::up),
    };

    poisson_error = phi;

    for(int loop = 1; loop <= loopnum; ++loop) {
        //! is_not_boundaryは各方向の境界について
        //! 「その境界は計算空間の境界でない」か「その境界が計算空間の境界であり、周期境界である」場合にtrueとなるため
        //! 各方向の端要素の計算をIterationで計算する場合にチェックする必要がある
        for(int k = 1; k < cz_with_glue - 1; ++k){
            if((k != 1 || is_not_boundary[4]) && (k != cz_with_glue - 2 || is_not_boundary[5])) {
                for(int j = 1; j < cy_with_glue - 1; ++j){
                    if((j != 1 || is_not_boundary[2]) && (j != cy_with_glue - 2 || is_not_boundary[3])) {
                        for(int i = 1; i < cx_with_glue - 1; ++i){
                            if((i != 1 || is_not_boundary[0]) && (i != cx_with_glue - 2 || is_not_boundary[1])) {
                                phi[i][j][k] = (1.0 - omega) * phi[i][j][k] + omega*(phi[i+1][j][k] + phi[i-1][j][k] + phi[i][j+1][k] + phi[i][j-1][k] + phi[i][j][k+1] + phi[i][j][k-1] + rho_coeff * rho[0][i][j][k])/6.0;
                            }
                        }
                    }
                }
            }
        }

        //! @note: 実際には部分部分をSORで計算して送信というのを繰り返す方が収束効率がよい
        MPIw::Environment::sendRecvField(phi);

        if ( (loop % 10 == 0) && (this->checkPhiResidual() < required_error) ) {
            cout << "performed " << loop << " iterations." << endl;
            break;
        }
    }

    field->setBoundaryConditionPhi();

    //! 全グリッド上のエラーを更新
    for(int k = 0; k < cz_with_glue; ++k){
        for(int j = 0; j < cy_with_glue; ++j){
            for(int i = 0; i < cx_with_glue; ++i){
                poisson_error[i][j][k] = phi[i][j][k] - poisson_error[i][j][k];
            }
        }
    }
}

//! 電位分布の残差の最大ノルムを返す
double RootGrid::checkPhiResidual() {
    double residual = 0.0;
    const double normalized_eps = Normalizer::eps0;

    const bool is_not_boundary[6] = {
        Environment::isNotBoundary(AXIS::x, AXIS_SIDE::low),
        Environment::isNotBoundary(AXIS::x, AXIS_SIDE::up),
        Environment::isNotBoundary(AXIS::y, AXIS_SIDE::low),
        Environment::isNotBoundary(AXIS::y, AXIS_SIDE::up),
        Environment::isNotBoundary(AXIS::z, AXIS_SIDE::low),
        Environment::isNotBoundary(AXIS::z, AXIS_SIDE::up),
    };

    auto& phi = field->getPhi();
    auto& rho = field->getRho();
    auto& poisson_residual = field->getPoissonResidual();

    const int cx_with_glue = phi.shape()[0];
    const int cy_with_glue = phi.shape()[1];
    const int cz_with_glue = phi.shape()[2];

    for(int k = 1; k < cz_with_glue - 1; ++k){
        if((k != 1 || is_not_boundary[4]) && (k != cz_with_glue - 2 || is_not_boundary[5])) {
            for(int j = 1; j < cy_with_glue - 1; ++j){
                if((j != 1 || is_not_boundary[2]) && (j != cy_with_glue - 2 || is_not_boundary[3])) {
                    for(int i = 1; i < cx_with_glue - 1; ++i){
                        if((i != 1 || is_not_boundary[0]) && (i != cx_with_glue - 2 || is_not_boundary[1])) {
                            double source_value = rho[0][i][j][k]/normalized_eps;
                            double tmp_res = field->poissonOperator(phi, i, j, k) + source_value;
                            poisson_residual[i][j][k] = tmp_res;

                            if (fabs(tmp_res / source_value) > residual) {
                                residual = fabs(tmp_res / source_value);
                            }
                        }
                    }
                }
            }
        }
    }

    if(MPIw::Environment::numprocs > 1) {
        residual = MPIw::Environment::Comms["world"].max(residual);
    }
    return residual;
}

//! @brief 差分法で電場を更新する
//! e = - (p_+1 - p_+0)/dx
void RootGrid::updateEfield() {
    auto& ex = field->getEx();
    auto& ey = field->getEy();
    auto& ez = field->getEz();
    auto& phi = field->getPhi();

    const int cx_with_glue = ex.shape()[0] + 1; // nx + 2
    const int cy_with_glue = ey.shape()[1] + 1;
    const int cz_with_glue = ez.shape()[2] + 1;
    const double per_dx = 1.0 / dx;

    //! @note:隣と通信しなくてもいい
    //! phiが通信してあるため、端の要素を通信なしで計算可能
    for(int i = 0; i < cx_with_glue; ++i){
        for(int j = 0; j < cy_with_glue; ++j){
            for(int k = 0; k < cz_with_glue; ++k){
                //! 各方向には1つ少ないのでcx-1まで
                if(i < cx_with_glue - 1) ex[i][j][k] = (phi[i][j][k] - phi[i + 1][j][k]) * per_dx;
                if(j < cy_with_glue - 1) ey[i][j][k] = (phi[i][j][k] - phi[i][j + 1][k]) * per_dx;
                if(k < cz_with_glue - 1) ez[i][j][k] = (phi[i][j][k] - phi[i][j][k + 1]) * per_dx;
            }
        }
    }

    this->updateReferenceEfield();

    for(auto& child : children) {
        child->updateEfield();
    }
}

//! Reference Efield (ノード上で定義される電場) を更新する
void RootGrid::updateReferenceEfield() {
    auto& ex = field->getEx();
    auto& ey = field->getEy();
    auto& ez = field->getEz();
    auto& phi = field->getPhi();
    auto& exref = field->getExRef();
    auto& eyref = field->getEyRef();
    auto& ezref = field->getEzRef();
    const int cx_with_glue = ex.shape()[0] + 1; // nx + 2
    const int cy_with_glue = ey.shape()[1] + 1;
    const int cz_with_glue = ez.shape()[2] + 1;

    //! reference 更新
    for(int i = 1; i < cx_with_glue - 1; ++i){
        for(int j = 1; j < cy_with_glue - 1; ++j){
            for(int k = 1; k < cz_with_glue - 1; ++k){
                exref[i][j][k] = 0.5 * (ex[i-1][j][k] + ex[i][j][k]);
                eyref[i][j][k] = 0.5 * (ey[i][j-1][k] + ey[i][j][k]);
                ezref[i][j][k] = 0.5 * (ez[i][j][k-1] + ez[i][j][k]);
            }
        }
    }

    //! 外側境界の条件設定
    if (Environment::isBoundary(AXIS::x, AXIS_SIDE::low)) {
        for(int j = 1; j < cy_with_glue - 1; ++j){
            for(int k = 1; k < cz_with_glue - 1; ++k){
                exref[0][j][k] = 0.0;

                // Boundary である場合、更新 (これ正しい??)
                exref[1][j][k] = 0.5 * (phi[3][j][k] - 4.0 * phi[2][j][k] + 3.0 * phi[1][j][k]);
            }
        }
    }

    if (Environment::isBoundary(AXIS::x, AXIS_SIDE::up)) {
        for(int j = 1; j < cy_with_glue - 1; ++j){
            for(int k = 1; k < cz_with_glue - 1; ++k){
                exref[cx_with_glue - 1][j][k] = 0.0;
                exref[cx_with_glue - 2][j][k] = 0.5 * (-phi[cx_with_glue - 4][j][k] + 4.0 * phi[cx_with_glue - 3][j][k] - 3.0 * phi[cx_with_glue - 2][j][k]);
            }
        }
    }

    if (Environment::isBoundary(AXIS::y, AXIS_SIDE::low)) {
        for(int i = 1; i < cx_with_glue - 1; ++i){
            for(int k = 1; k < cz_with_glue - 1; ++k){
                eyref[i][0][k] = 0.0;
                eyref[i][1][k] = 0.5 * (phi[i][3][k] - 4.0 * phi[i][2][k] + 3.0 * phi[i][1][k]);
            }
        }
    }

    if (Environment::isBoundary(AXIS::y, AXIS_SIDE::up)) {
        for(int i = 1; i < cx_with_glue - 1; ++i){
            for(int k = 1; k < cz_with_glue - 1; ++k){
                eyref[i][cy_with_glue - 1][k] = 0.0;
                eyref[i][cy_with_glue - 2][k] = 0.5 * (-phi[i][cy_with_glue - 4][k] + 4.0 * phi[i][cy_with_glue - 3][k] - 3.0 * phi[i][cy_with_glue - 2][k]);
            }
        }
    }

    if (Environment::isBoundary(AXIS::z, AXIS_SIDE::low)) {
        for(int i = 1; i < cx_with_glue - 1; ++i){
            for(int j = 1; j < cy_with_glue - 1; ++j){
                ezref[i][j][0] = 0.0;
                ezref[i][j][1] = 0.5 * (phi[i][j][3] - 4.0 * phi[i][j][2] + 3.0 * phi[i][j][1]);
            }
        }
    }

    if (Environment::isBoundary(AXIS::z, AXIS_SIDE::up)) {
        for(int i = 1; i < cx_with_glue - 1; ++i){
            for(int j = 1; j < cy_with_glue - 1; ++j){
                ezref[i][j][cz_with_glue - 1] = 0.0;
                ezref[i][j][cz_with_glue - 2] = 0.5 * (-phi[i][j][cz_with_glue - 4] + 4.0 * phi[i][j][cz_with_glue - 3] - 3.0 * phi[i][j][cz_with_glue - 2]);
            }
        }
    }

    MPIw::Environment::sendRecvField(exref);
    MPIw::Environment::sendRecvField(eyref);
    MPIw::Environment::sendRecvField(ezref);
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

    //! 子の電荷更新
    for(auto& child : children) {
        child->updateRho();
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

void RootGrid::insertAMRBlockInfo(SimpleVTK& vtk_gen, const std::string& data_type_name, const std::string& i_timestamp) const {
    vtk_gen.beginBlock(0);
        vtk_gen.beginDataSet(id);
        vtk_gen.setAMRBoxNode(from_ix, to_ix, from_iy, to_iy, from_iz, to_iz);
        vtk_gen.setFile(data_type_name + "_id_" + std::to_string(id) + "_" + i_timestamp + ".vti");
        vtk_gen.endDataSet();
    vtk_gen.endBlock();

    for(const auto& child : children) {
        child->insertAMRBlockInfo(vtk_gen, data_type_name, i_timestamp);
    }
}