#include "grid.hpp"
#include "normalizer.hpp"
#include "dataio.hpp"

#define USE_BOOST
#include "simple_vtk.hpp"

//! Glueセルを含まないデータを生成
boost::multi_array<float, 3> Grid::getTrueCells(const tdArray& x3D, const double unnorm) const {
    // ZONECENTなので-1する
    const int xsize = this->getXNodeSize() - 1;
    const int ysize = this->getYNodeSize() - 1;
    const int zsize = this->getZNodeSize() - 1;

    boost::multi_array<float, 3> true_cells(boost::extents[xsize][ysize][zsize]);

    for(int k = 1; k < zsize + 1; ++k){
        for(int j = 1; j < ysize + 1; ++j){
            for(int i = 1; i < xsize + 1; ++i){
                true_cells[i - 1][j - 1][k - 1] = static_cast<float>(x3D[i][j][k] * unnorm);
            }
        }
    }

    // RVO
    return true_cells;
}

//! Glueノードを含まないデータを生成
boost::multi_array<float, 3> Grid::getTrueNodes(const tdArray& x3D, const double unnorm) const {
    int xsize = this->getXNodeSize();
    int ysize = this->getYNodeSize();
    int zsize = this->getZNodeSize();

    boost::multi_array<float, 3> true_nodes(boost::extents[xsize][ysize][zsize]);

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

//! Glueノードを含まないデータを生成
boost::multi_array<float, 4> Grid::getTrueNodeVectors(const tdArray& xvector, const tdArray& yvector, const tdArray& zvector, const double unnorm) const {
    int xsize = this->getXNodeSize();
    int ysize = this->getYNodeSize();
    int zsize = this->getZNodeSize();

    boost::multi_array<float, 4> true_nodes(boost::extents[3][xsize][ysize][zsize]);

    for(int k = 1; k < zsize + 1; ++k){
        for(int j = 1; j < ysize + 1; ++j){
            for(int i = 1; i < xsize + 1; ++i){
                true_nodes[0][i - 1][j - 1][k - 1] = static_cast<float>(xvector[i][j][k] * unnorm);
                true_nodes[1][i - 1][j - 1][k - 1] = static_cast<float>(yvector[i][j][k] * unnorm);
                true_nodes[2][i - 1][j - 1][k - 1] = static_cast<float>(zvector[i][j][k] * unnorm);
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

    boost::multi_array<float, 3> true_nodes(boost::extents[xsize][ysize][zsize]);

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

void Grid::plotFieldData(const std::string& data_type_name, const std::string& i_timestamp) const {
    SimpleVTK gen;
    gen.enableExtentManagement();
    gen.changeBaseExtent(0, this->getXNodeSize() - 1, 0, this->getYNodeSize() - 1, 0, this->getZNodeSize() - 1);
    gen.changeBaseOrigin(Normalizer::unnormalizeLength(base_x), Normalizer::unnormalizeLength(base_y), Normalizer::unnormalizeLength(base_z));
    const auto base_spacing = Normalizer::unnormalizeLength(dx);
    gen.changeBaseSpacing(base_spacing, base_spacing, base_spacing);
    gen.setInnerElementPerLine(100);

    gen.beginVTK("ImageData");
    gen.setVersion("0.1");
    gen.setByteOrder("LittleEndian");
        gen.beginContentWithPiece();
            if (data_type_name == "potential" || data_type_name == "rho") {
                gen.beginPointData();
                gen.setScalars(data_type_name);
                    gen.beginDataArray(data_type_name, "Float32", "ascii");
                        if (data_type_name == "potential") {
                            auto values = this->getTrueNodes(field->getPhi(), Normalizer::unnormalizePotential(1.0));
                            gen.addMultiArray(values);
                        } else {
                            auto values = this->getTrueNodes(field->getRho(), 0, Normalizer::unnormalizeRho(1.0));
                            gen.addMultiArray(values);
                        }
                    gen.endDataArray();
                gen.endPointData();
            } else if (data_type_name == "efield" || data_type_name == "bfield" || data_type_name == "current") {
                gen.beginPointData();
                gen.setVectors(data_type_name);
                    gen.beginDataArray(data_type_name, "Float32", "ascii");
                    gen.setNumberOfComponents("3");
                        if (data_type_name == "efield") {
                            auto values = this->getTrueNodeVectors(field->getExRef(), field->getEyRef(), field->getEzRef(), Normalizer::unnormalizeEfield(1.0));
                            gen.addMultiArray(values);
                        } else if (data_type_name == "bfield") {
                            auto values = this->getTrueNodeVectors(field->getBxRef(), field->getByRef(), field->getBzRef(), Normalizer::unnormalizeBfield(1.0));
                            gen.addMultiArray(values);
                        } else {
                            auto values = this->getTrueNodeVectors(field->getJxRef(), field->getJyRef(), field->getJzRef(), Normalizer::unnormalizeCurrent(1.0));
                            gen.addMultiArray(values);
                        }
                    gen.endDataArray();
                gen.endPointData();
            } else if (data_type_name == "density") {
                gen.beginCellData();
                std::string pnames = "";
                for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
                    pnames += Environment::getParticleType(pid)->getName();

                    if (pid != Environment::num_of_particle_types - 1) pnames += " ";
                }
                gen.setScalars(pnames);
                const auto& density = field->getDensity();
                const auto unnorm = Normalizer::unnormalizeDensity(1.0);

                for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
                    gen.beginDataArray(Environment::getParticleType(pid)->getName(), "Float32", "ascii");
                        auto values = this->getTrueCells(density[pid], unnorm);
                        gen.addMultiArray(values);
                    gen.endDataArray();
                }
                gen.endCellData();
            }
        gen.endContentWithPiece();
    gen.endVTK();

    std::string file_name = "data/raw_data/" + data_type_name + "_id_" + std::to_string(id) + "_" + i_timestamp;
    gen.generate(file_name);

    for(const auto& child : children) {
        child->plotFieldData(data_type_name, i_timestamp);
    }
}