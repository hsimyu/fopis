#ifndef __TDPIC_RANDOM_DISTRIBUTION_H_INCLUDED__
#define __TDPIC_RANDOM_DISTRIBUTION_H_INCLUDED__

#include <random>

#define _USE_MATH_DEFINES
#include <cmath>

namespace RandomDistribution {
    class CosineEmission {
        private:
            std::mt19937_64 mt_azimuth_angle;
            std::mt19937_64 mt_depression_angle;
            std::uniform_real_distribution<> dist_uniform{0.0, 1.0};

        public:
            CosineEmission(const int seed) : mt_azimuth_angle(586936 + seed), mt_depression_angle(8785997666 + seed) {}

            //! [-pi/2, pi/2]の cosine/2 分布 (法線からの角度)
            //! acos(2x-1) - pi/2 = -asin(2x-1)
            double genDepressionAngle() {
                return std::asin(2.0 * dist_uniform(mt_depression_angle) - 1.0);
            }

            //! [0, 2pi]の一様分布
            double genAzimuthAngle() {
                return 2.0 * M_PI * dist_uniform(mt_azimuth_angle);
            }
    };
}

#endif