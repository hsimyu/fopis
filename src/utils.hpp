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

    // 軸方向を指定する文字列をindexに変換する
    int getAxisIndex(const std::string);
    int getLowOrUpIndex(const std::string);

    float* getTrueEdges2(tdArray const&, tdArray const&, tdArray const&);
    float* getTrueEdges(tdArray const&, const int);
    float* getTrueFaces(const tdArray&, const int);
    void convert1Dto3Darray(double*, const int, const int, const int, tdArray&);
    void clearBoundaryValues(tdArray&, const int, const int, const int);
    void initializeTdarray(tdArray&);

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

            //! V = m^2 kg / s^3 A = m^2 kg C / s^2 -> 1
            static double normalizePotential(double raw_phi) {
                return raw_phi * pow(t_unit, 2) / (pow(x_unit, 2) * m_unit * e_unit);
            }

            //! 1 -> m^2 kg / s^2 C = V
            static double unnormalizePotential(double normalized_phi) {
                return normalized_phi * (pow(x_unit, 2) * m_unit) / (pow(t_unit, 2) * e_unit);
            }

            //! V/m = m kg / s^3 A = m kg / s^2 C -> 1
            static double normalizeEfield(double raw_efield) {
                return raw_efield * pow(t_unit, 2) * e_unit / (x_unit * m_unit);
            }

            //! 1 -> m kg / C s^2 = V / m
            static double unnormalizeEfield(double normalized_efield) {
                return normalized_efield * (x_unit * m_unit) / (pow(t_unit, 2) * e_unit);
            }

            //! kg m^2 / s^2 -> 1
            static double normalizeEnergy(double raw_energy) {
                return raw_energy * pow(t_unit, 2) / (m_unit * pow(x_unit, 2));
            }

            //! 1 -> kg m^2 / s^2
            static double unnormalizeEnergy(double normalized_energy) {
                return normalized_energy * m_unit * pow(x_unit, 2) / pow(t_unit, 2);
            }

            //! s^4 A^2 / kg m^3 = C^2 s^2 / kg m^3 -> 1
            static double normalizeEpsilon(double raw_epsilon) {
                return raw_epsilon * (pow(x_unit, 3) * m_unit) / (pow(e_unit, 2) * pow(t_unit, 2));
            }

            //! 1 -> C^2 s^2 / kg m^3
            static double unnormalizeEpsilon(double normalized_eps) {
                return normalized_eps * (pow(e_unit, 2) * pow(t_unit, 2)) / (pow(x_unit, 3) * m_unit);
            }

            //! m kg / s^2 A^2 = m kg / C^2 -> 1
            static double normalizeMu(double raw_mu) {
                return raw_mu * pow(e_unit, 2) / (m_unit * x_unit);
            }

            //! 1 -> m kg / C^2
            static double unnormalizeMu(double normalized_mu) {
                return normalized_mu * m_unit * x_unit / pow(e_unit, 2);
            }

            static double normalizeDensity(double raw_density) {
                return raw_density * pow(x_unit, 3);
            }

            static double unnormalizeDensity(double normalized_density) {
                return normalized_density / pow(x_unit, 3);
            }

            static double normalizeFlux(double raw_flux) {
                return raw_flux * pow(x_unit, 2) * t_unit;
            }

            static double unnormalizeFlux(double normalized_flux) {
                return normalized_flux / (pow(x_unit, 2) * t_unit);
            }

            static double normalizeArea(double raw_area) {
                return raw_area / pow(x_unit, 2);
            }

            static double unnormalizeArea(double normalized_area) {
                return normalized_area * pow(x_unit, 2);
            }

            static double normalizeFrequency(double raw_freq) {
                return raw_freq * t_unit;
            }

            static double unnormalizeFrequency(double normalized_freq) {
                return normalized_freq / t_unit;
            }
    };
}


#endif
