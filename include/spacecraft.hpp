#ifndef __TDPIC_SPACECRAFT_H_INCLUDED__
#define __TDPIC_SPACECRAFT_H_INCLUDED__
#include "global.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include "position.hpp"
#include "field.hpp"

using ObjectDefinedMapBool = boost::multi_array<bool, 3>;
using ObjectDefinedMapInt = boost::multi_array<int, 3>;
using ObjectNodes = std::map< unsigned int, std::array<int, 3> >;
using ObjectNodeTextures = std::map< unsigned int, std::vector<int> >;
using ObjectCells = std::vector<std::array<int, 3>>;
using PropertyPair = std::map<std::string, double>;

struct ObjectDataFromFile {
    ObjectNodes nodes;
    ObjectNodeTextures textures;
    ObjectCells cells;
};

class Particle;

//! @class: Spacecraft
class Spacecraft {
private:
    static const std::map<std::string, PropertyPair> material_property_list;
    static unsigned int num_of_spacecraft;

    //! 計算中での物体名
    std::string name;

    //! 実際に読み込むobjファイル名
    std::string file_name;

    std::string surface_type;
    size_t num_cmat;
    bool is_defined_in_this_process;
    double potential;
    double potential_fix;
    double total_charge;
    std::vector<double> current;

    struct LocalParticleEmissionInfo {
        Position relative_emission_position;
        std::array<double, 3> emission_vector;
    };

    std::map<int, LocalParticleEmissionInfo> emit_particle_info;
    std::map<int, std::string> material_names;
    std::map<int, double> material_capacitances;

    //! セルベースのオブジェクト定義マップ
    ObjectDefinedMapBool object_cell_map;

    //! ノードベースの電荷定義マップ
    RhoArray charge_map;

    //! キャパシタンス行列
    using Cmatrix = boost::numeric::ublas::matrix<double>;
    Cmatrix capacity_matrix;

    //! キャパシタンス行列の番号と対応する物体の位置を格納する
    std::map<size_t, Position> capacity_matrix_relation;
    //! キャパシタンス行列の番号と対応する静電容量の値を格納する
    std::map<size_t, double> capacitance_map;

    //! 再計算しないために初期化後に保持する
    double total_cmat_value;

    //! コンストラクタ内部処理共通化用
    void construct(const size_t, const size_t, const size_t, const ObjectInfo_t&, const ObjectNodes&, const ObjectCells&, const ObjectNodeTextures&);

    //! 電荷の総量が変化していないかの check 用
    auto getTotalCharge(const RhoArray&) const;

    //! 実際の電荷配分関数
    void distributeInnerParticleChargeToXsurface(const Position& pos, const int id, const double charge);
    void distributeInnerParticleChargeToYsurface(const Position& pos, const int id, const double charge);
    void distributeInnerParticleChargeToZsurface(const Position& pos, const int id, const double charge);

public:
    Spacecraft(const size_t nx, const size_t ny, const size_t nz,
        const unsigned int _num_cmat, const ObjectInfo_t& obj_info,
        const ObjectNodes& nodes, const ObjectCells& cells, const ObjectNodeTextures& textures) :
        name(obj_info.name),
        surface_type(obj_info.surface_type),
        num_cmat(_num_cmat),
        current{},
        emit_particle_info{},
        material_names{obj_info.materials},
        material_capacitances{},
        object_cell_map(boost::extents[0][0][0]),
        charge_map{},
        capacity_matrix(0, 0),
        capacity_matrix_relation{},
        capacitance_map{} {
        construct(nx, ny, nz, obj_info, nodes, cells, textures);
    }

    // アクセサ
    void setName(const std::string _name) { name = _name; }
    std::string getName() const { return name; }
    std::string getFileName() const { return file_name; }

    auto getPotential(void) const { return potential; }
    auto getPotentialFix(void) const { return potential_fix; }
    void setPotentialFix(const double val) { potential_fix = val; }

    auto getCmatSize(void) const { return num_cmat; }
    Position getCmatPos(const unsigned int);

    auto getCmatValue(const unsigned int col, const unsigned int row) const { return capacity_matrix(col, row); }
    void setCmatValue(const unsigned int col, const unsigned int row, const double value) {
        capacity_matrix(col, row) = value;
    };
    void updateTotalCmatValue();

    // 判定用関数
    bool isMyCmat(const unsigned int cmat_number) const { return (capacity_matrix_relation.count(cmat_number) > 0); }
    bool isDefined(void) const { return is_defined_in_this_process; }
    bool isContaining(const Particle&) const;
    bool isContaining(const Position&) const;
    bool isDielectricSurface() const { return (surface_type == "dielectric"); }

    //! 粒子放出用
    void emitParticles(ParticleArray& parray);
    bool hasEmitParticles() const {return (emit_particle_info.size() > 0);}
    bool isValidEmission(Particle& p) const;
    void subtractChargeOfParticle(const Particle& p);

    // その他ユーティリティ関数
    void makeCmatrixInvert(void);
    void removeInnerParticle(Particle&) const;
    void distributeInnerParticleCharge(Particle&);
    void applyCharge(RhoArray&) const;
    void redistributeCharge(RhoArray&, const tdArray&);
    void resetCurrent();
    void sumCurrent();

    // IO関数
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
