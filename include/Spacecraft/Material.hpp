#ifndef __TDPIC_SPACECRAFT_MATERIAL_H_INCLUDED__
#define __TDPIC_SPACECRAFT_MATERIAL_H_INCLUDED__

//! 物体表面素材の情報を保持する構造体
class MaterialInfo_t {
    public:
        MaterialInfo_t(const std::string _name) : name{_name} {}

        //! 電磁的な性質
        double relative_permittivity = 1.0;
        double conductivity = 1.0;

        //! 二次電子に関する係数
        double atomic_number = 1.0;
        double fermi_energy = 1.0;
        double epsi_max = 1.0;
        double delta_max = 1.0;

        std::string name;
};

#endif