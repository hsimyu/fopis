#ifndef __TDPIC_SPACECRAFT_H_INCLUDED__
#define __TDPIC_SPACECRAFT_H_INCLUDED__
#include "global.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include "position.hpp"

using ObjectDefinedMap = boost::multi_array<bool, 3>;
using ObjectNodes = std::map< unsigned int, std::array<int, 3> >;
using ObjectFaces = std::vector< std::array<int, 4> >;

struct ObjectDataFromFile {
    ObjectNodes nodes;
    ObjectFaces faces;
};

class Particle;

//! @class: Spacecraft
class Spacecraft {
private:
    static unsigned int num_of_spacecraft;
    std::string name;
    size_t num_cmat;
    bool is_defined_in_this_process;
    double potential;
    double potential_fix;
    double total_charge;

    //! オブジェクト定義マップとキャパシティ定義マップ
    ObjectDefinedMap object_node_map;
    ObjectDefinedMap object_xface_map;
    ObjectDefinedMap object_yface_map;
    ObjectDefinedMap object_zface_map;
    tdArray charge_map;

    //! キャパシタンス行列
    using Cmatrix = boost::numeric::ublas::matrix<double>;
    Cmatrix capacity_matrix;

    //! キャパシタンス行列の番号と対応する物体の位置を格納する
    std::map<size_t, Position> capacity_matrix_relation;

    //! 再計算しないために初期化後に保持する
    double total_cmat_value;

    //! コンストラクタ内部処理共通化用
    void construct(const size_t, const size_t, const size_t, const ObjectNodes&, const ObjectNodes&, const ObjectFaces&);

    //! 電荷の総量が変化していないかの check 用
    auto getTotalCharge(const tdArray&) const;
public:
    Spacecraft(const size_t nx, const size_t ny, const size_t nz,
        const unsigned int _num_cmat, const std::string _name,
        const ObjectNodes& nodes, const ObjectNodes& glue_nodes,
        const ObjectFaces& faces) :
        name(_name),
        num_cmat(_num_cmat),
        object_node_map(boost::extents[0][0][0]),
        object_xface_map(boost::extents[0][0][0]),
        object_yface_map(boost::extents[0][0][0]),
        object_zface_map(boost::extents[0][0][0]),
        charge_map(boost::extents[0][0][0]),
        capacity_matrix(0, 0) {
        construct(nx, ny, nz, nodes, glue_nodes, faces);
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
    void setTotalCmatValue(const double val) { total_cmat_value = val; }

    // 判定用関数
    bool isMyCmat(const unsigned int cmat_number) const { return (capacity_matrix_relation.count(cmat_number) > 0); }
    bool isDefined(void) const { return is_defined_in_this_process; }
    bool isIncluded(const Particle&) const;

    // その他ユーティリティ関数
    void makeCmatrixInvert(void);
    void removeInnerParticle(Particle&) const;
    void distributeInnerParticleCharge(Particle&);
    void applyCharge(tdArray&) const;
    void redistributeCharge(tdArray&, const tdArray&);

    // IO関数
    std::string getLogHeader() const;
    std::string getLogEntry() const;
    friend std::ostream& operator<<(std::ostream&, const Spacecraft&);
};

namespace ObjectUtils {
    ObjectDataFromFile getObjectNodesFromObjFile(const std::string& obj_file_name);
}
#endif
