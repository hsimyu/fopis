#include "grid.hpp"
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

}

void RootGrid::saveResumeData() {
    using H5F = HighFive::File;

    const std::string file_name = "resume/grid_id_" + std::to_string(id) + ".h5";
    H5F file(file_name, H5F::ReadWrite | H5F::Create | H5F::Truncate);

    //! 不要な粒子は削除しておく
    this->removeInvalidParticles();
    this->saveResumeParticleData(file);
    this->saveResumeFieldData(file);
}

void RootGrid::saveResumeParticleData(HighFive::File& file) const {
    //! HighFiveはユーザー定義型の格納ができないようなので、
    //! 各粒子の位置速度の配列を抽出して格納する
    for (size_t pid = 0; pid < particles.size(); ++pid) {
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

        const std::string data_set_name = Environment::getParticleType(pid)->getName();

        auto data_set = file.createDataSet<double>(data_set_name + "_x", HighFive::DataSpace::From(particle_x));
        data_set.write(particle_x);
        data_set = file.createDataSet<double>(data_set_name + "_y", HighFive::DataSpace::From(particle_y));
        data_set.write(particle_y);
        data_set = file.createDataSet<double>(data_set_name + "_z", HighFive::DataSpace::From(particle_z));
        data_set.write(particle_z);

        data_set = file.createDataSet<double>(data_set_name + "_vx", HighFive::DataSpace::From(particle_vx));
        data_set.write(particle_vx);
        data_set = file.createDataSet<double>(data_set_name + "_vy", HighFive::DataSpace::From(particle_vy));
        data_set.write(particle_vy);
        data_set = file.createDataSet<double>(data_set_name + "_vz", HighFive::DataSpace::From(particle_vz));
        data_set.write(particle_vz);
    }
}

void RootGrid::saveResumeFieldData(HighFive::File& file) const {
    {
        auto& ex = field->getEx();
        auto& ey = field->getEy();
        auto& ez = field->getEz();

        auto data_set = file.createDataSet<double>("ex", HighFive::DataSpace::From(ex));
        data_set.write(ex);
        data_set = file.createDataSet<double>("ey", HighFive::DataSpace::From(ey));
        data_set.write(ey);
        data_set = file.createDataSet<double>("ez", HighFive::DataSpace::From(ez));
        data_set.write(ez);
    }

    {
        auto& bx = field->getBx();
        auto& by = field->getBy();
        auto& bz = field->getBz();
        auto data_set = file.createDataSet<double>("bx", HighFive::DataSpace::From(bx));
        data_set.write(bx);
        data_set = file.createDataSet<double>("by", HighFive::DataSpace::From(by));
        data_set.write(by);
        data_set = file.createDataSet<double>("bz", HighFive::DataSpace::From(bz));
        data_set.write(bz);
    }
}