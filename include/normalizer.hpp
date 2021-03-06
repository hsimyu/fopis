#ifndef __TDPIC_NORMALIZER_H_INCLUDED__
#define __TDPIC_NORMALIZER_H_INCLUDED__
#include "global.hpp"

//! 物理量を無次元量として扱う際の正規化を担当するクラス
class Normalizer {
protected:
    // インスタンス化しない
    Normalizer();
    ~Normalizer();

    //! @note: initializeの際に次元単位をセットする必要がある
    static double x_unit;
    static double t_unit;
    static double m_unit;
    static double e_unit;

    //! @note: よく使う単位は事前計算しておくと早い
    static double V_unit;
    static double A_unit;
    static double x_pow_3;

    //! 単位系が変更された際に内部の値を更新する
    static void updatePhysicalConstants(void) {
        x_pow_3 = pow(x_unit, 3);

        //! V = m^2 kg / C s^2
        V_unit = pow(x_unit, 2) * m_unit / (e_unit * pow(t_unit, 2));

        //! A = C / s
        A_unit = e_unit / t_unit;

        eps0 = normalizeEpsilon(::eps0);
        mu0 = normalizeMu(::mu0);
        c = normalizeVelocity(::c);
    }

public:
    //! 正規化された物理定数
    static double eps0;
    static double mu0;
    static double c;

    static double getEps0(const int grid_division_level) {
        return pow(pow(2, grid_division_level), 3) * eps0;
    }

    static void setLengthUnit(const double _u) {
        x_unit = _u;
        updatePhysicalConstants();
    }

    static void setTimeUnit(const double _u) {
        t_unit = _u;
        updatePhysicalConstants();
    }

    static void setMassUnit(const double _u) {
        m_unit = _u;
        updatePhysicalConstants();
    }

    static void setChargeUnit(const double _u) {
        e_unit = _u;
        updatePhysicalConstants();
    }

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

    //! 1 -> kg m / s^2 = N
    static double unnormalizeForce(double normalized_force) {
        return normalized_force * m_unit * x_unit / (t_unit * t_unit);
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

    //! V = m^2 kg / s^3 A = m^2 kg / C s^2 -> 1
    static double normalizePotential(double raw_phi) {
        return raw_phi / V_unit;
    }

    //! 1 -> m^2 kg / s^2 C = V
    static double unnormalizePotential(double normalized_phi) {
        return normalized_phi * V_unit;
    }

    //! C / s -> 1
    static double normalizeCurrent(double raw_current) {
        return raw_current / A_unit;
    }

    //! 1 -> C / s
    static double unnormalizeCurrent(double normalized_current) {
        return normalized_current * A_unit;
    }

    //! C / s m^2 -> 1
    static double normalizeCurrentDensity(double raw_current_density) {
        return raw_current_density * (x_unit * x_unit) / A_unit;
    }

    //! 1 -> C / s m^2
    static double unnormalizeCurrentDensity(double normalized_current_density) {
        return normalized_current_density * A_unit / (x_unit * x_unit);
    }

    //! 1 -> C / m^3
    static double unnormalizeRho(double normalized_rho) {
        return normalized_rho * e_unit / x_pow_3;
    }

    //! V/m = m kg / s^3 A = m kg / s^2 C -> 1
    static double normalizeEfield(double raw_efield) {
        return raw_efield * x_unit / V_unit;
    }

    //! 1 -> m kg / C s^2 = V / m
    static double unnormalizeEfield(double normalized_efield) {
        return normalized_efield * V_unit / x_unit;
    }

    //! T = kg / C s -> 1
    static double normalizeBfield(double raw_bfield) {
        return raw_bfield * (t_unit * e_unit) / m_unit;
    }

    //! 1 -> T = kg / C s
    static double unnormalizeBfield(double normalized_bfield) {
        return normalized_bfield * m_unit / (t_unit * e_unit);
    }

    //! kg m^2 / s^2 -> 1
    static double normalizeEnergy(double raw_energy) {
        return raw_energy * pow(t_unit, 2) / (m_unit * pow(x_unit, 2));
    }

    //! 1 -> kg m^2 / s^2
    static double unnormalizeEnergy(double normalized_energy) {
        return normalized_energy * m_unit * pow(x_unit, 2) / pow(t_unit, 2);
    }

    //! 1 -> kg m^2 / s^2 = J, or eV (J / C)
    static double unnormalizeEnergy(double normalized_energy, std::string unit_type) {
        if (unit_type == "eV" || unit_type == "electronvolt" || unit_type == "ev") {
            return normalized_energy * m_unit * pow(x_unit, 2) / (e_unit * pow(t_unit, 2));
        } else {
            return normalized_energy * m_unit * pow(x_unit, 2) / pow(t_unit, 2);
        }
    }

    //! s^4 A^2 / kg m^3 = C^2 s^2 / kg m^3 -> 1
    static double normalizeEpsilon(double raw_epsilon) {
        return raw_epsilon * (x_pow_3 * m_unit) / (pow(e_unit, 2) * pow(t_unit, 2));
    }

    //! 1 -> C^2 s^2 / kg m^3
    static double unnormalizeEpsilon(double normalized_eps) {
        return normalized_eps * (pow(e_unit, 2) * pow(t_unit, 2)) / (x_pow_3 * m_unit);
    }

    //! F = C / V -> 1
    static double normalizeCapacitance(double raw_capacitance) {
        return raw_capacitance * V_unit / e_unit;
    }

    //! 1 -> C / V = F
    static double unnormalizeCapacitance(double normalized_capacitance) {
        return normalized_capacitance * e_unit / V_unit;
    }

    //! 1 / Ωm = A / V m -> 1
    static double normalizeConductivity(double raw_conductivity) {
        return raw_conductivity * V_unit * x_unit / A_unit;
    }

    //! 1 -> A / V m
    static double unnormalizeConductivity(double normalized_conductivity) {
        return normalized_conductivity * A_unit / (V_unit * x_unit);
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
        return raw_density * x_pow_3;
    }

    static double unnormalizeDensity(double normalized_density) {
        return normalized_density / x_pow_3;
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

    static double convertEnergyToVelocity(const double normalized_energy) {
        return std::sqrt(2.0 * normalized_energy);
    }

    static double convertEnergyToVelocity(const double normalized_energy, const double mass) {
        return std::sqrt(2.0 * normalized_energy / mass);
    }
};

#endif
