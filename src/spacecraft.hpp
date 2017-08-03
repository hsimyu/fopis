#ifndef __TDPIC_SPACECRAFT_H_INCLUDED__
#define __TDPIC_SPACECRAFT_H_INCLUDED__
#include "global.hpp"

class Position;
class Particle;

//! @class: Spacecraft
class Spacecraft {
private:
    static unsigned int num_of_spacecraft;
    std::string name;
    size_t num_cmat;
    objectArray object_map;
    tdArray charge_map;

    //! キャパシタンス行列
    using Cmatrix = boost::multi_array<double, 2>;
    Cmatrix capacity_matrix;

    //! キャパシタンス行列の番号と対応する物体の位置を格納する
    std::map<size_t, int[3]> capacity_matrix_relation;

    //! コンストラクタ共通化用(必要か?)
    void construct(const size_t, const size_t, const size_t);

public:
    Spacecraft(const size_t nx, const size_t ny, const size_t nz) :
        name("spacecraft" + std::to_string(num_of_spacecraft)),
        object_map(boost::extents[0][0][0]),
        charge_map(boost::extents[0][0][0]) {
        construct(nx, ny, nz);
    }

    Spacecraft(const size_t nx, const size_t ny, const size_t nz, const std::string _name) :
        name(_name),
        object_map(boost::extents[0][0][0]),
        charge_map(boost::extents[0][0][0]) {
        construct(nx, ny, nz);
    }

    void setName(const std::string _name) { name = _name; }
    std::string getName() const { return name; }

    auto getCmatSize(void) const { return num_cmat; }

    void setCmatValue(const size_t col, const size_t row, const double value) {
        capacity_matrix[col][row] = value;
    };

    bool isIncluded(const Particle&) const;
    void removeInnerParticle(Particle&) const;
    void distributeInnerParticleCharge(Particle&);
    void applyCharge(tdArray&) const;
    friend std::ostream& operator<<(std::ostream&, const Spacecraft&);
};
#endif
