#include "grid.hpp"
#include "field.hpp"
#include "utils.hpp"
#include "dataio.hpp"

void RootGrid::initializeObject(void) {
    if (Environment::isRootNode) cout << "-- Defining Objects -- " << endl;

    for (const auto& object_info : Environment::objects_info) {
        std::string obj_name = object_info.name;
        //! 物体関連の設定を関連付けされた obj 形式ファイルから読み込む
        ObjectDataFromFile object_data = ObjectUtils::getObjectNodesFromObjFile(object_info.file_name);

        size_t num_cmat = object_data.nodes.size();
        const auto& whole_node_array = object_data.nodes;

        //! innerと判定されたやつだけ渡す
        ObjectNodes inner_node_array;
        ObjectNodes glue_node_array;
        bool is_object_in_this_node = false;

        for(const auto& node_pair : whole_node_array) {
            const auto cmat_itr = node_pair.first;
            const auto& node_pos = node_pair.second;

            const auto i = node_pos[0];
            const auto j = node_pos[1];
            const auto k = node_pos[2];

            if (isInnerNode(i, j, k)) {
                is_object_in_this_node = true;
                inner_node_array[cmat_itr] = Environment::getRelativePositionOnRootWithGlue(i, j, k);
            } else if (isGlueNode(i, j, k)) {
                glue_node_array[cmat_itr] = Environment::getRelativePositionOnRootWithGlue(i, j, k);
            }
        }

        const auto& cell_array = object_data.cells;
        ObjectCells inner_cell_array;
        for(const auto& cell_pos : cell_array) {
            if (isInnerCellWithGlue(cell_pos[0], cell_pos[1], cell_pos[2])) {
                auto rel_pos = Environment::getRelativePositionOnRootWithGlue(cell_pos[0], cell_pos[1], cell_pos[2]);
                inner_cell_array.push_back( {{rel_pos[0], rel_pos[1], rel_pos[2]}} );
            }
        }

        //! 物体定義点がゼロでも Spacecraft オブジェクトだけは作成しておいた方がよい
        //! emplace_back で Spacecraft object を直接構築
        objects.emplace_back(
            nx, ny, nz, num_cmat, object_info,
            inner_node_array, glue_node_array, whole_node_array,
            inner_cell_array, object_data.textures, object_data.connected_list
        );

        //! Comm作成 (物体が入っていないならnullになる)
        MPIw::Environment::makeNewComm(obj_name, is_object_in_this_node);

        auto& obj = objects[ objects.size() - 1 ];
        if (obj.isDefined()) {
            obj.initializeAfterMakeComm();

            if (MPIw::Environment::isRootNode(obj_name)) {
                cout << Environment::rankStr() << "is set to Root Node for " << obj_name << "." << endl;
                cout << obj << endl;
            }
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
        if (Environment::getOptions().useExistingCapacityMatrix()) {
            const auto is_valid_load = IO::loadCmatrixData(obj);
            if (is_valid_load) continue;
        }

        //! データを使わずに初期化
        const auto num_cmat = obj.getCmatSize();

        std::unique_ptr<Utils::ProgressManager> pm;
        if (MPIw::Environment::isRootNode(obj.getName())) {
            std::unique_ptr<Utils::ProgressManager> tmp_pm(new Utils::ProgressManager(num_cmat, "cmat_solve"));
            pm = std::move(tmp_pm);
        }

        for(unsigned int cmat_col_itr = 0; cmat_col_itr < num_cmat; ++cmat_col_itr ) {
            if (MPIw::Environment::isRootNode(obj.getName())) pm->update(cmat_col_itr);

            //! 該当する頂点に単位電荷を付与
            if (obj.isMyCmat(cmat_col_itr)) {
                const auto& cmat_pos = obj.getCmatPos(cmat_col_itr);
                rho[0][cmat_pos.i][cmat_pos.j][cmat_pos.k] = 1.0;
            }

            solvePoisson();

            //#pragma omp parallel for ordered shared(obj)
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

    for(auto& obj : objects) {
        //! 初期電荷オフセットを付与する
        obj.initializeChargeMapOffset(phi);
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

void RootGrid::emitParticlesFromObjects(void) {
    for(auto& obj : objects) {
        if (obj.isDefined() && obj.hasEmitParticles()) {
            obj.emitParticles(particles);
        }
    }
}
