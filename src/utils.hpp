#ifndef __TDPIC_UTILS_H_INCLUDED__
#define __TDPIC_UTILS_H_INCLUDED__
#include <string>
#include <picojson.h>
#include "global.hpp"

namespace Utils {
    void printTotalMemory(void);
    std::string prettyMemoryString(double);
    std::string readFile(const std::string&);
    picojson::value::object readJSONFile(const std::string&);

    float* getTrueCells(const tdArray&);
    float* getTrueEdges(const tdArray&, const int);
    float* getTrueFaces(const tdArray&, const int);
    void convert1Dto3Darray(double*, const int, const int, const int, tdArray&);
    void clearBoundaryValues(tdArray&, const int, const int, const int);

    void createDir(std::string);

    //! 物理量を無次元量として扱う際の正規化を担当するクラス
    class Normalizer {
        protected:
            Normalizer();
            ~Normalizer();

        public:
            //! @note: initializeの際に次元単位をセットする必要がある
            static double x_unit;
            static double t_unit;
            static double m_unit;
            static double e_unit;

            // Normalize Utilities
            static double normalizeLength(double raw_x) {
                return raw_x / x_unit;
            }

            static double unnormalizeLength(double normalized_x) {
                return normalized_x * x_unit;
            }

            static double normalizeVelocity(double raw_v) {
                return t_unit * raw_v / x_unit;
            }

            static double unnormalizeVelocity(double normalized_v) {
                return normalized_v * x_unit / t_unit;
            }

            static double normalizeTime(double raw_t) {
                return raw_t / t_unit;
            }

            static double unnormalizeTime(double normalized_t) {
                return normalized_t * t_unit;
            }

            static double normalizeCharge(double raw_e) {
                return raw_e / e_unit;
            }

            static double unnormalizeCharge(double normalized_e) {
                return normalized_e * e_unit;
            }

            static double normalizeMass(double raw_mass) {
                return raw_mass / m_unit;
            }

            static double unnormalizeMass(double normalized_mass) {
                return normalized_mass * m_unit;
            }

            //! kg*m^2/s^2 -> 1;
            static double normalizeEnergy(double raw_energy) {
                return raw_energy * pow(t_unit, 2) / (m_unit * pow(x_unit, 2));
            }

            //! 1 -> kg*m^2/s^2;
            static double unnormalizeEnergy(double normalized_energy) {
                return normalized_energy * m_unit * pow(x_unit, 2) / pow(t_unit, 2); //kg * m^2/s^2;
            }

            //! s^4*A^2/kg*m^3 == C^2*s^2/kg*m^3 -> 1;
            static double normalizeEpsilon(double raw_epsilon) {
                return raw_epsilon * (pow(x_unit, 3) * m_unit) / (pow(e_unit, 2) * pow(t_unit, 2));
            }

            //! 1 -> C^2*s^2/kg*m^3;
            static double unnormalizeEpsilon(double normalized_eps) {
                return normalized_eps * (pow(e_unit, 2) * pow(t_unit, 2)) / (pow(x_unit, 3) * m_unit);
            }

    };
}


#endif
