#ifndef __TDPIC_SPACECRAFT_CMATRIX_H_INCLUDED__
#define __TDPIC_SPACECRAFT_CMATRIX_H_INCLUDED__

#include <Eigen/Dense>

//! キャパシタンス行列の実装
class Cmatrix {
    private:
        using CmatrixImpl = Eigen::MatrixXd;
        CmatrixImpl capacity_matrix;

    public:
        Cmatrix(const int size) : capacity_matrix{size, size} {}

        template<typename T>
        double& operator() (const T col, const T row) {
            return capacity_matrix(col, row);
        }

        // const operator
        template<typename T>
        double operator() (const T col, const T row) const {
            return capacity_matrix(col, row);
        }

        template<typename T>
        void resize(const T size) {
            capacity_matrix.resize(size, size);
        }

        Cmatrix& inverse() {
            capacity_matrix = capacity_matrix.inverse();
            return *this;
        }

        double total() const {
            double res = 0.0;
            for(size_t col = 0; col < capacity_matrix.cols(); ++col) {
                for(size_t row = 0; row < capacity_matrix.rows(); ++row) {
                    res += capacity_matrix(col, row);
                }
            }
            return res;
        }
};

#endif