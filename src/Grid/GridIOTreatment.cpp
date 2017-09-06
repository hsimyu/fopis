#include "grid.hpp"
#include "normalizer.hpp"
#include "dataio.hpp"

//! for DATA IO
//! 粒子の位置から密度を計算する
boost::multi_array<float, 3> Grid::getDensity(const int pid) const {
    // ZONECENTなので-1する
    const int xsize = this->getXNodeSize() - 1;
    const int ysize = this->getYNodeSize() - 1;
    const int zsize = this->getZNodeSize() - 1;

    boost::multi_array<float, 3> zones(boost::extents[xsize][ysize][zsize], boost::fortran_storage_order());
    const auto size = static_cast<float>(Normalizer::unnormalizeDensity(Environment::getParticleType(pid)->getSize()));

    for(int pnum = 0; pnum < particles[pid].size(); ++pnum){
        const Particle& p = particles[pid][pnum];

        if(p.isValid) {
            Position pos(p);
            zones[pos.i - 1][pos.j - 1][pos.k - 1] += size;
        }
    }

    // RVO
    return zones;
}

//! Glueノードを含まないデータを生成
boost::multi_array<float, 3> Grid::getTrueNodes(const tdArray& x3D, const double unnorm) const {
    int xsize = this->getXNodeSize();
    int ysize = this->getYNodeSize();
    int zsize = this->getZNodeSize();

    boost::multi_array<float, 3> true_nodes(boost::extents[xsize][ysize][zsize], boost::fortran_storage_order());

    for(int k = 1; k < zsize + 1; ++k){
        for(int j = 1; j < ysize + 1; ++j){
            for(int i = 1; i < xsize + 1; ++i){
                true_nodes[i - 1][j - 1][k - 1] = static_cast<float>(x3D[i][j][k] * unnorm);
            }
        }
    }

    // RVO
    return true_nodes;
}

//! RhoArray用
boost::multi_array<float, 3> Grid::getTrueNodes(const RhoArray& rho, const int pid, const double unnorm) const {
    int xsize = this->getXNodeSize();
    int ysize = this->getYNodeSize();
    int zsize = this->getZNodeSize();

    boost::multi_array<float, 3> true_nodes(boost::extents[xsize][ysize][zsize], boost::fortran_storage_order());

    for(int k = 1; k < zsize + 1; ++k){
        for(int j = 1; j < ysize + 1; ++j){
            for(int i = 1; i < xsize + 1; ++i){
                true_nodes[i - 1][j - 1][k - 1] = static_cast<float>(rho[pid][i][j][k] * unnorm);
            }
        }
    }

    // RVO
    return true_nodes;
}

// -- DATA IO methods --
void Grid::putFieldData(HighFive::Group& group, const std::string& data_type_name, const std::string& i_timestamp) const {
    auto getGroup = [](auto& g, const std::string& group_name) {
        if (g.exist(group_name)) {
            return g.getGroup(group_name);
        } else {
            return g.createGroup(group_name);
        }
    };

    const std::string& level_str = "level" + std::to_string(level);
    HighFive::Group local_group = getGroup(group, level_str);

    if (data_type_name == "potential") {
        auto values = this->getTrueNodes(field->getPhi(), Normalizer::unnormalizePotential(1.0));
        auto dataset = local_group.createDataSet<float>(data_type_name, HighFive::DataSpace::From(values));
        dataset.write(values);
    } else if(data_type_name == "rho") {
        auto values = this->getTrueNodes(field->getRho(), 0, Normalizer::unnormalizeRho(1.0));
        auto dataset = local_group.createDataSet<float>(data_type_name, HighFive::DataSpace::From(values));
        dataset.write(values);
    } else if(data_type_name == "efield") {
        HighFive::Group data_type_group = getGroup(local_group, data_type_name);

        const std::array<std::string, 3> axis{{"ex", "ey", "ez"}};
        for(const auto& group_name : axis) {
            if (group_name == "ex") {
                auto values = this->getTrueNodes(field->getExRef(), Normalizer::unnormalizeEfield(1.0));
                auto dataset = local_group.createDataSet<float>(group_name, HighFive::DataSpace::From(values));
                dataset.write(values);
            } else if (group_name == "ey") {
                auto values = this->getTrueNodes(field->getEyRef(), Normalizer::unnormalizeEfield(1.0));
                auto dataset = local_group.createDataSet<float>(group_name, HighFive::DataSpace::From(values));
                dataset.write(values);
            } else if (group_name == "ez") {
                auto values = this->getTrueNodes(field->getEzRef(), Normalizer::unnormalizeEfield(1.0));
                auto dataset = local_group.createDataSet<float>(group_name, HighFive::DataSpace::From(values));
                dataset.write(values);
            }
        }
    } else if(data_type_name == "bfield") {
        HighFive::Group data_type_group = getGroup(local_group, data_type_name);
        const std::array<std::string, 3> axis{{"bx", "by", "bz"}};

        for(const auto& group_name : axis) {
            if (group_name == "bx") {
                auto values = this->getTrueNodes(field->getBxRef(), Normalizer::unnormalizeBfield(1.0));
                auto dataset = local_group.createDataSet<float>(group_name, HighFive::DataSpace::From(values));
                dataset.write(values);
            } else if (group_name == "by") {
                auto values = this->getTrueNodes(field->getByRef(), Normalizer::unnormalizeBfield(1.0));
                auto dataset = local_group.createDataSet<float>(group_name, HighFive::DataSpace::From(values));
                dataset.write(values);
            } else if (group_name == "bz") {
                auto values = this->getTrueNodes(field->getBzRef(), Normalizer::unnormalizeBfield(1.0));
                auto dataset = local_group.createDataSet<float>(group_name, HighFive::DataSpace::From(values));
                dataset.write(values);
            }
        }
    } else if(data_type_name == "density") {
        HighFive::Group data_type_group = getGroup(local_group, data_type_name);
        for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
            const std::string& pname = Environment::getParticleType(pid)->getName();
            auto values = this->getDensity(pid);
            auto dataset = data_type_group.createDataSet<float>(pname, HighFive::DataSpace::From(values));
            dataset.write(values);
        }
    }

    if (data_type_name == "potential" || data_type_name == "rho") {
        for (const auto& child : children) {
            child->putFieldData(group, data_type_name, i_timestamp);
        }
    }
}
