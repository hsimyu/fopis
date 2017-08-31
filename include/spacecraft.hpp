#ifndef __TDPIC_SPACECRAFT_H_INCLUDED__
#define __TDPIC_SPACECRAFT_H_INCLUDED__
#include "global.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include "position.hpp"
#include "field.hpp"

using ObjectDefinedMapBool = boost::multi_array<bool, 3>;
using ObjectDefinedMapInt = boost::multi_array<int, 3>;
using ObjectNodes = std::map< unsigned int, std::array<int, 3> >;
using ObjectCells = std::vector<std::array<int, 4>>;
using PropertyPair = std::map<std::string, double>;

struct ObjectDataFromFile {
    ObjectNodes nodes;
    ObjectCells cells;
};

class Particle;

//! @class: Spacecraft
class Spacecraft {
private:
    static const std::map<std::string, PropertyPair> material_property_list;
    static unsigned int num_of_spacecraft;
    std::string name;
    std::string surface_type;
    size_t num_cmat;
    bool is_defined_in_this_process;
    double potential;
    double potential_fix;
    double total_charge;
    std::vector<double> current;
    std::vector<int> emit_particle_ids;
    std::map<int, std::string> materials;

    //! オブジェクト定義マップとキャパシティ定義マップ
    ObjectDefinedMapBool object_node_map;
    ObjectDefinedMapInt object_cell_map;
    RhoArray charge_map;

    //! キャパシタンス行列
    using Cmatrix = boost::numeric::ublas::matrix<double>;
    Cmatrix capacity_matrix;

    //! キャパシタンス行列の番号と対応する物体の位置を格納する
    std::map<size_t, Position> capacity_matrix_relation;

    //! 再計算しないために初期化後に保持する
    double total_cmat_value;

    //! コンストラクタ内部処理共通化用
    void construct(const size_t, const size_t, const size_t, const ObjectInfo_t&, const ObjectNodes&, const ObjectNodes&, const ObjectCells&);

    //! 電荷の総量が変化していないかの check 用
    auto getTotalCharge(const RhoArray&) const;

public:
    Spacecraft(const size_t nx, const size_t ny, const size_t nz,
        const unsigned int _num_cmat, const ObjectInfo_t& obj_info,
        const ObjectNodes& nodes, const ObjectNodes& glue_nodes, const ObjectCells& cells) :
        name(obj_info.name),
        surface_type(obj_info.surface_type),
        num_cmat(_num_cmat),
        current{},
        emit_particle_ids{},
        materials{obj_info.materials},
        object_node_map(boost::extents[0][0][0]),
        object_cell_map(boost::extents[0][0][0]),
        charge_map{},
        capacity_matrix(0, 0) {
        construct(nx, ny, nz, obj_info, nodes, glue_nodes, cells);
    }

    // アクセサ
    void setName(const std::string _name) { name = _name; }
    std::string getName() const { return name; }

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
    bool hasEmitParticles() const {return (emit_particle_ids.size() > 0);}
    bool isValidEmission(Particle& p) const;
    void subtractChargeOfParticle(const Particle& p);

    // その他ユーティリティ関数
    void makeCmatrixInvert(void);
    void removeInnerParticle(Particle&) const;
    void distributeInnerParticleCharge(Particle&);
    void applyCharge(RhoArray&) const;
    void redistributeCharge(RhoArray&, const tdArray&);
    void resetCurrent();

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
