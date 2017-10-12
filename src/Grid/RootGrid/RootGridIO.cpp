#include "grid.hpp"
#include "utils.hpp"
#include "dataio.hpp"

#define H5_USE_BOOST
#include <highfive/H5File.hpp>

#define USE_BOOST
#include "simple_vtk.hpp"

void RootGrid::insertAMRBlockInfo(SimpleVTK& vtk_gen, const std::string& data_type_name, const std::string& i_timestamp) const {
    vtk_gen.beginBlock(0);
        vtk_gen.beginDataSet(id);
        vtk_gen.setAMRBoxNode(from_ix, from_ix + this->getXNodeSize() - 1, from_iy, from_iy + this->getYNodeSize() - 1, from_iz, from_iz + this->getZNodeSize() - 1);
        vtk_gen.setFile("raw_data/" + data_type_name + "_id_" + std::to_string(id) + "_" + i_timestamp + ".vti");
        vtk_gen.endDataSet();
    vtk_gen.endBlock();

    for(const auto& child : children) {
        child->insertAMRBlockInfo(vtk_gen, data_type_name, i_timestamp);
    }
}

void RootGrid::loadResumeData() {
    using H5F = HighFive::File;

    const std::string file_name = "resume/grid_id_" + std::to_string(id) + ".h5";

    if (Utils::isExistingFile(file_name)) {
        H5F file(file_name, H5F::ReadWrite);

        this->loadResumeParticleData(file);
        this->loadResumeFieldData(file);
        this->loadResumeObjectData(file);
    }
}

void RootGrid::saveResumeData() {
    using H5F = HighFive::File;

    const std::string file_name = "resume/grid_id_" + std::to_string(id) + ".h5";
    H5F file(file_name, H5F::ReadWrite | H5F::Create | H5F::Truncate);

    //! 不要な粒子は削除しておく
    this->removeInvalidParticles();
    this->saveResumeParticleData(file);
    this->saveResumeFieldData(file);
    this->saveResumeObjectData(file);
}

void RootGrid::loadResumeParticleData(HighFive::File& file) {
    for (size_t pid = 0; pid < particles.size(); ++pid) {
        const std::string pname = Environment::getParticleType(pid)->getName();
        auto group = file.getGroup(pname);

        auto& parray = particles[pid];

        //! 初期化時に入力された粒子は全てerase
        parray.erase(parray.begin(), parray.end());

        size_t size;
        {
            auto dataset = group.getDataSet("size");
            dataset.read(size);
        }

        if (size > 0) {
            //! 擬似乱数生成カウント
            {
                std::vector<size_t> counts;
                auto data_set = group.getDataSet("generated_counts");
                data_set.read(counts);
                Environment::getParticleType(pid)->proceedGeneratedCounts(counts);
            }

            std::vector<double> particle_x(size);
            std::vector<double> particle_y(size);
            std::vector<double> particle_z(size);
            {
                auto dataset = group.getDataSet("x");
                dataset.read(particle_x);
            }

            {
                auto dataset = group.getDataSet("y");
                dataset.read(particle_y);
            }

            {
                auto dataset = group.getDataSet("z");
                dataset.read(particle_z);
            }

            std::vector<double> particle_vx(size);
            std::vector<double> particle_vy(size);
            std::vector<double> particle_vz(size);

            {
                auto dataset = group.getDataSet("vx");
                dataset.read(particle_vx);
            }

            {
                auto dataset = group.getDataSet("vy");
                dataset.read(particle_vy);
            }

            {
                auto dataset = group.getDataSet("vz");
                dataset.read(particle_vz);
            }

            parray.resize(size);

            #pragma omp parallel for
            for(size_t pnum = 0; pnum < size; ++pnum) {
                parray[pnum].x = particle_x[pnum];
                parray[pnum].y = particle_y[pnum];
                parray[pnum].z = particle_z[pnum];
                parray[pnum].vx = particle_vx[pnum];
                parray[pnum].vy = particle_vy[pnum];
                parray[pnum].vz = particle_vz[pnum];
                parray[pnum].makeValid();
            }
        }
    }
}

void RootGrid::saveResumeParticleData(HighFive::File& file) const {
    //! HighFiveはユーザー定義型の格納ができないようなので、
    //! 各粒子の位置速度の配列を抽出して格納する
    for (size_t pid = 0; pid < particles.size(); ++pid) {
        const std::string pname = Environment::getParticleType(pid)->getName();
        auto group = file.createGroup(pname);

        const auto& parray = particles[pid];
        const auto size = parray.size();
        std::vector<double> particle_x(size);
        std::vector<double> particle_y(size);
        std::vector<double> particle_z(size);
        std::vector<double> particle_vx(size);
        std::vector<double> particle_vy(size);
        std::vector<double> particle_vz(size);

        #pragma omp parallel for
        for(size_t pnum = 0; pnum < size; ++pnum) {
            particle_x[pnum] = parray[pnum].x;
            particle_y[pnum] = parray[pnum].y;
            particle_z[pnum] = parray[pnum].z;
            particle_vx[pnum] = parray[pnum].vx;
            particle_vy[pnum] = parray[pnum].vy;
            particle_vz[pnum] = parray[pnum].vz;
        }

        //! 粒子数
        auto data_set = group.createDataSet<double>("size", HighFive::DataSpace::From(size));
        data_set.write(size);

        //! 擬似乱数生成カウント
        auto counts = Environment::getParticleType(pid)->getGeneratedCounts();
        data_set = group.createDataSet<size_t>("generated_counts", HighFive::DataSpace::From(counts));
        data_set.write(counts);

        data_set = group.createDataSet<double>("x", HighFive::DataSpace::From(particle_x));
        data_set.write(particle_x);
        data_set = group.createDataSet<double>("y", HighFive::DataSpace::From(particle_y));
        data_set.write(particle_y);
        data_set = group.createDataSet<double>("z", HighFive::DataSpace::From(particle_z));
        data_set.write(particle_z);

        data_set = group.createDataSet<double>("vx", HighFive::DataSpace::From(particle_vx));
        data_set.write(particle_vx);
        data_set = group.createDataSet<double>("vy", HighFive::DataSpace::From(particle_vy));
        data_set.write(particle_vy);
        data_set = group.createDataSet<double>("vz", HighFive::DataSpace::From(particle_vz));
        data_set.write(particle_vz);
    }
}

void RootGrid::loadResumeFieldData(HighFive::File& file) {
    {
        auto phi_group = file.getGroup("Potential");
        {
            auto& phi = field->getPhi();
            auto data_set = phi_group.getDataSet("phi");
            data_set.read(phi);
        }
    }

    {
        auto efield_group = file.getGroup("Efield");
        {
            auto& ex = field->getEx();
            auto data_set = efield_group.getDataSet("ex");
            data_set.read(ex);
        }

        {
            auto& ey = field->getEy();
            auto data_set = efield_group.getDataSet("ey");
            data_set.read(ey);
        }

        {
            auto& ez = field->getEz();
            auto data_set = efield_group.getDataSet("ez");
            data_set.read(ez);
        }
    }

    {
        auto bfield_group = file.getGroup("Bfield");

        {
            auto& bx = field->getBx();
            auto data_set = bfield_group.getDataSet("bx");
            data_set.read(bx);
        }

        {
            auto& by = field->getBy();
            auto data_set = bfield_group.getDataSet("by");
            data_set.read(by);
        }

        {
            auto& bz = field->getBz();
            auto data_set = bfield_group.getDataSet("bz");
            data_set.read(bz);
        }
    }
}

void RootGrid::saveResumeFieldData(HighFive::File& file) const {
    {
        auto phi_group = file.createGroup("Potential");

        auto& phi = field->getPhi();
        auto data_set = phi_group.createDataSet<double>("phi", HighFive::DataSpace::From(phi));
        data_set.write(phi);
    }

    {
        auto efield_group = file.createGroup("Efield");

        auto& ex = field->getEx();
        auto& ey = field->getEy();
        auto& ez = field->getEz();

        auto data_set = efield_group.createDataSet<double>("ex", HighFive::DataSpace::From(ex));
        data_set.write(ex);
        data_set = efield_group.createDataSet<double>("ey", HighFive::DataSpace::From(ey));
        data_set.write(ey);
        data_set = efield_group.createDataSet<double>("ez", HighFive::DataSpace::From(ez));
        data_set.write(ez);
    }

    {
        auto bfield_group = file.createGroup("Bfield");

        auto& bx = field->getBx();
        auto& by = field->getBy();
        auto& bz = field->getBz();
        auto data_set = bfield_group.createDataSet<double>("bx", HighFive::DataSpace::From(bx));
        data_set.write(bx);
        data_set = bfield_group.createDataSet<double>("by", HighFive::DataSpace::From(by));
        data_set.write(by);
        data_set = bfield_group.createDataSet<double>("bz", HighFive::DataSpace::From(bz));
        data_set.write(bz);
    }
}

void RootGrid::loadResumeObjectData(HighFive::File& file) {
    for(auto& object : objects) {
        auto group = file.getGroup(object.getName());

        auto data_set = group.getDataSet("charge_map");
        data_set.read(object.getChargeMap());
    }
}

void RootGrid::saveResumeObjectData(HighFive::File& file) const {
    for(const auto& object : objects) {
        auto group = file.createGroup(object.getName());

        auto data_set = group.createDataSet<double>("charge_map", HighFive::DataSpace::From(object.getChargeMap()));
        data_set.write(object.getChargeMap());
    }
}