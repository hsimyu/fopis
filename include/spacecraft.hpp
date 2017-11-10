#ifndef __TDPIC_SPACECRAFT_H_INCLUDED__
#define __TDPIC_SPACECRAFT_H_INCLUDED__

#include <array>
#include "global.hpp"
#include "field.hpp"
#include "position.hpp"
#include <Spacecraft/Cmatrix.hpp>
#include <Spacecraft/Material.hpp>
#include <Spacecraft/Incident.hpp>

using ObjectDefinedMapBool = boost::multi_array<bool, 3>;
using ObjectDefinedMapInt = boost::multi_array<int, 3>;
using ObjectNodes = std::map< unsigned int, std::array<int, 3> >;
using ObjectNodeTextures = std::map< unsigned int, std::vector<int> >;
using ObjectChargeMap = boost::multi_array<double, 2>;
using ObjectConnectivityList = std::vector< std::vector<unsigned int> >;
using ObjectCells = std::vector<std::array<int, 3>>;
using PropertyPair = std::map<std::string, double>;

struct ObjectDataFromFile {
    ObjectNodes nodes;
    ObjectNodeTextures textures;
    ObjectCells cells;
    ObjectConnectivityList connected_list;
};

class ParticleType;
class Particle;
class SimpleVTK;
class Grid;

//! @class: Spacecraft
class Spacecraft {
private:
    Grid* parent_grid;

    static const std::map<std::string, PropertyPair> material_property_list;
    static unsigned int num_of_spacecraft;

    //! 計算中での物体名
    std::string name;

    //! 実際に読み込むobjファイル名
    std::string file_name;

    std::string surface_type;
    size_t num_cmat;
    bool is_defined_in_this_process;
    unsigned int plot_potential_mapping_width;
    double potential;
    double total_charge;
    std::vector<double> current;

    bool is_potential_fixed;
    double fixed_potential;

    double initial_potential_offset;

    struct LocalParticleEmissionInfo {
        Position relative_emission_position;
        std::array<double, 3> emission_vector;
        double emission_process_number = -1.0;
    };

    std::map<int, LocalParticleEmissionInfo> emit_particle_info;
    std::vector<IncidentInfo_t> incident_events;
    std::map<int, std::string> material_names;
    std::map<int, double> material_capacitances;

    //! セルベースのオブジェクト定義マップ
    ObjectDefinedMapBool object_cell_map;

    //! 電荷定義マップ([particle_id, cmat_index])
    ObjectChargeMap charge_map;
    ObjectChargeMap temporary_charge_map;

    //! キャパシタンス行列
    Cmatrix capacity_matrix;

    //! 自分が担当するキャパシタンス行列の番号と対応する物体の位置を格納する
    std::map<size_t, Position> capacity_matrix_relation;

    //! Glueキャパシタンス行列の番号と対応する物体の位置を格納する
    std::map<size_t, Position> glue_capacity_matrix_relation;

    //! 全体のキャパシタンス行列の番号と対応する物体の位置を格納する
    std::map<size_t, Position> whole_capacity_matrix_relation;

    //! キャパシタンス行列の番号と対応する静電容量の値を格納する
    std::map<size_t, double> capacitance_map;

    //! 再計算しないために初期化後に保持する
    double total_cmat_value;

    //! コンストラクタ内部処理共通化用
    void construct(
        const size_t, const size_t, const size_t, const ObjectInfo_t&,
        const ObjectNodes&, const ObjectNodes&, const ObjectNodes&,
        const ObjectCells&, const ObjectNodeTextures&, const ObjectConnectivityList&
    );

    //! 表面電荷の総量
    auto getTotalCharge() const;
    double getMaxPotential() const;
    double getMinPotential() const;
    double getNodeCharge(const unsigned int cmat_itr) const;

    //! 実際の電荷配分関数
    void distributeInnerParticleChargeToXsurface(const Position& pos, const int id, const double charge);
    void distributeInnerParticleChargeToYsurface(const Position& pos, const int id, const double charge);
    void distributeInnerParticleChargeToZsurface(const Position& pos, const int id, const double charge);

    //! 実際の電荷再配分関数
    void redistributeChargeForPerfectConductor(RhoArray& rho, const tdArray& phi);
    void redistributeChargeForDielectric(RhoArray& rho, const tdArray& phi);

    //! 粒子放出時の電荷計算用
    void subtractChargeOfParticleFromXsurface(const Position& pos, const int id, const double charge);
    void subtractChargeOfParticleFromYsurface(const Position& pos, const int id, const double charge);
    void subtractChargeOfParticleFromZsurface(const Position& pos, const int id, const double charge);

    //! 粒子放出に参加するプロセスの数を計算
    void initializeEmissionParticleInfo();

    //! 接続リストは物体が定義されたプロセスが持つことにしてしまう
    ObjectConnectivityList connected_list;

    //! VTKに表面マップなどを追加するための関数
    void insertPoints(SimpleVTK& gen) const;
    void insertConnectivity(SimpleVTK& gen) const;
    void insertOffsets(SimpleVTK& gen) const;
    void insertTypes(SimpleVTK& gen) const;

    //! 表面にマップするデータを物体定義済みプロセスから取得する
    template<typename T>
    std::vector<T> getWholePotentialMap(const tdArray& phi) const;
    void saveWholeNodePositions(const ObjectNodes& whole_nodes);

    //! 表面上の点かどうかを確認する
    bool isXsurfaceMinus(const int i, const int j, const int k) const;
    bool isXsurfacePlus(const int i, const int j, const int k) const;
    bool isYsurfaceMinus(const int i, const int j, const int k) const;
    bool isYsurfacePlus(const int i, const int j, const int k) const;
    bool isZsurfaceMinus(const int i, const int j, const int k) const;
    bool isZsurfacePlus(const int i, const int j, const int k) const;

    bool isXsurfacePoint(const Position& pos, const int sign) const;
    bool isXsurfaceCmatNode(const Position& pos, const int sign) const;
    bool isXsurfaceCmatNode(const Position& pos) const;

    bool isYsurfacePoint(const Position& pos, const int sign) const;
    bool isZsurfacePoint(const Position& pos, const int sign) const;

public:
    Spacecraft(
        const size_t nx, const size_t ny, const size_t nz,
        const unsigned int _num_cmat,
        const ObjectInfo_t& obj_info,
        const ObjectNodes& nodes,
        const ObjectNodes& glue_nodes,
        const ObjectNodes& whole_nodes,
        const ObjectCells& cells,
        const ObjectNodeTextures& textures,
        const ObjectConnectivityList& clist) :
        name(obj_info.name),
        surface_type(obj_info.surface_type),
        num_cmat(_num_cmat),
        current{},
        emit_particle_info{},
        material_names{obj_info.materials},
        material_capacitances{},
        object_cell_map(boost::extents[0][0][0]),
        charge_map{},
        temporary_charge_map{},
        capacity_matrix{0},
        capacity_matrix_relation{},
        glue_capacity_matrix_relation{},
        whole_capacity_matrix_relation{},
        capacitance_map{} {
        construct(nx, ny, nz, obj_info, nodes, glue_nodes, whole_nodes, cells, textures, clist);
    }

    //! MPI Comm 生成後に呼び出す必要がある初期化
    void initializeAfterMakeComm();

    //! 初期電荷オフセットを付与
    void initializeChargeMapOffset(const tdArray& phi);

    // アクセサ
    void setName(const std::string _name) { name = _name; }
    std::string getName() const { return name; }
    std::string getFileName() const { return file_name; }
    void setParent(Grid* g) { parent_grid = g; }

    auto getPotential(void) const { return potential; }
    auto getFixedPotential(void) const { return fixed_potential; }
    void setFixedPotential(const double val) { fixed_potential = val; }

    auto getCmatSize(void) const { return num_cmat; }
    Position getCmatPos(const unsigned int);
    unsigned int getCmatNumber(const int i, const int j, const int k) const;
    unsigned int getCmatNumber(const Position& pos) const;
    bool isCmatNode(const int i, const int j, const int k) const;

    auto getCmatValue(const unsigned int col, const unsigned int row) const {
        return capacity_matrix(col, row);
    }

    void setCmatValue(const unsigned int col, const unsigned int row, const double value) {
        capacity_matrix(col, row) = value;
    };
    void updateTotalCmatValue();

    auto getChargeMap() const {
        return charge_map;
    }

    auto& getChargeMap() {
        return charge_map;
    }

    // 判定用関数
    template<typename T>
    bool isMyCmat(const T cmat_number) const { return (capacity_matrix_relation.count(cmat_number) > 0); }

    bool isDefined(void) const { return is_defined_in_this_process; }
    bool isContaining(const Particle&) const;
    bool isContaining(const Position&) const;
    bool isDielectricSurface() const { return (surface_type == "dielectric"); }

    bool isPlotTiming(const unsigned int timestep) const { return (plot_potential_mapping_width > 0 && (timestep % plot_potential_mapping_width == 0)); }

    //! 粒子放出用
    void emitParticles(ParticleArray& parray, const double unit_length);
    bool hasEmitParticles() const {return (emit_particle_info.size() > 0);}
    bool hasSecondaryParticles() const;

    //! 二次電子管理系
    void addIncidentEvent(const std::shared_ptr<ParticleType> ptype_ptr, const double remaining_time, const Position& incident_pos, const Velocity& incident_vel, AXIS axis);
    void clearIncidentEvents();

    bool isValidEmission(Particle& p) const;

    // その他ユーティリティ関数
    void makeCmatrixInvert(void);
    void removeInnerParticle(Particle&) const;
    void distributeInnerParticleCharge(Particle&);
    void distributeInnerParticleChargeForSecondary(Particle& p, AXIS axis);
    void sumWholeCharge();

    //! 電荷再配分
    void redistributeCharge(RhoArray& rho, const tdArray& phi);

    void resetCurrent();
    void sumCurrent();

    // IO関数
    void plotPotentialMapping(const int timestep, const tdArray& phi) const;
    std::string getLogHeader() const;
    std::string getLogEntry() const;
    friend std::ostream& operator<<(std::ostream&, const Spacecraft&);
};

namespace ObjectUtils {
    ObjectDataFromFile getObjectNodesFromObjFile(const std::string& obj_file_name);

    template<typename T>
    int getTextureIndex(const T&& counts);

    template<typename T>
    T determineOneValueFromFourElems(T v1, T v2, T v3, T v4);
}

#endif
