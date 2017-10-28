#include "normalizer.hpp"

//! Normalizer の static 変数の実体
double Normalizer::x_unit = 1.0;
double Normalizer::t_unit = 1.0;
double Normalizer::e_unit = 1.0;
double Normalizer::m_unit = 1.0;

double Normalizer::V_unit = 1.0;
double Normalizer::A_unit = 1.0;
double Normalizer::x_pow_3 = 1.0;

//! Normalizer の正規化された static 物理定数 の実体
double Normalizer::eps0 = ::eps0;
double Normalizer::mu0 = ::mu0;
double Normalizer::c = ::c;