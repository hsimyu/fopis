#ifndef __TDPIC_SPACECRAFT_H_INCLUDED__
#define __TDPIC_SPACECRAFT_H_INCLUDED__
#include "global.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include "position.hpp"

class Particle;

//! @class: Spacecraft
class Spacecraft {
private:
    static unsigned int num_of_spacecraft;
    std::string name;
    size_t num_cmat;
    objectArray object_map;
    tdArray charge_map;

    double potential;
    double potential_bias;

    //! キャパシタンス行列
    using Cmatrix = boost::numeric::ublas::matrix<double>;
    Cmatrix capacity_matrix;

    //! キャパシタンス行列の番号と対応する物体の位置を格納する
    std::map<size_t, Position> capacity_matrix_relation;

    //! 再計算しないために初期化後に保持する
    double total_cmat_value;

    //! コンストラクタ内部処理共通化用
    using ObjNodes = std::map< unsigned int, std::array<unsigned int, 3> >;
    void construct(const size_t, const size_t, const size_t, const ObjNodes&);

    //! 電荷の総量が変化していないかの check 用
    auto getTotalCharge(const tdArray&) const;
public:
    Spacecraft(const size_t nx, const size_t ny, const size_t nz, const ObjNodes& nodes) :
        name("Spacecraft_" + std::to_string(num_of_spacecraft)),
        object_map(boost::extents[0][0][0]),
        charge_map(boost::extents[0][0][0]),
        capacity_matrix(0, 0) {
        construct(nx, ny, nz, nodes);
    }

    Spacecraft(const size_t nx, const size_t ny, const size_t nz, const std::string _name, const ObjNodes& nodes) :
        name(_name),
        object_map(boost::extents[0][0][0]),
        charge_map(boost::extents[0][0][0]),
        capacity_matrix(0, 0) {
        construct(nx, ny, nz, nodes);
    }

    void setName(const std::string _name) { name = _name; }
    std::string getName() const { return name; }

    auto getPotential(void) const { return potential; }
    auto getPotentialBias(void) const { return potential_bias; }
    void setPotentialBias(const double val) { potential_bias = val; }

    auto getCmatSize(void) const { return num_cmat; }
    Position getCmatPos(const unsigned int);
    bool isMyCmat(const unsigned int cmat_number) const {
        return (capacity_matrix_relation.count(cmat_number) > 0);
    }

    auto getCmatValue(const unsigned int col, const unsigned int row) {
        return capacity_matrix(col, row);
    }

    void setCmatValue(const unsigned int col, const unsigned int row, const double value) {
        capacity_matrix(col, row) = value;
    };

    void setTotalCmatValue(const double val) { total_cmat_value = val; }
    void makeCmatrixInvert(void);

    bool isIncluded(const Particle&) const;
    void removeInnerParticle(Particle&) const;
    void distributeInnerParticleCharge(Particle&);
    void applyCharge(tdArray&) const;
    void redistributeCharge(tdArray&, const tdArray&);
    friend std::ostream& operator<<(std::ostream&, const Spacecraft&);
};
#endif
