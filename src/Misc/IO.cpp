#include "global.hpp"
#include "environment.hpp"
#include "dataio.hpp"
#include "grid.hpp"
#include "mpiw.hpp"
#include "normalizer.hpp"
#include "utils.hpp"

#define H5_USE_BOOST
#define USE_BOOST
#include <highfive/H5File.hpp>
#include <simple_xdmf.hpp>
#include <simple_vtk.hpp>

namespace IO {
    void writeDataInParallel(std::shared_ptr<const Grid> g, const int timestep, const std::string& data_type_name) {
        if (Environment::isRootNode) generateXdmf(timestep, data_type_name);

        const std::string i_timestamp = (format("%08d") % timestep).str();
        const std::string rank_str = "proc_" + (format("%08d") % MPIw::Environment::rank).str();
        const std::string file_name = "data/" + rank_str + ".h5";

        using H5F = HighFive::File;
        // initialize file pointer only once
        static H5F file(file_name, H5F::ReadWrite | H5F::Create | H5F::Truncate);

        HighFive::Group data_group;
        if (file.exist(i_timestamp)) {
            data_group = file.getGroup(i_timestamp);
        } else {
            data_group = file.createGroup(i_timestamp);
        }

        g->putFieldData(data_group, data_type_name, i_timestamp);
    }

    void writeCmatrixData(const Spacecraft& obj) {
        using H5F = HighFive::File;

        const std::string file_name = "data/objects/" + obj.getName() + ".h5";
        H5F file(file_name, H5F::ReadWrite | H5F::Create | H5F::Truncate);

        const auto num_of_cmatrix = obj.getCmatSize();
        boost::multi_array<double, 2> capacity_matrix_data(boost::extents[num_of_cmatrix][num_of_cmatrix]);

        for(unsigned int i = 0; i < num_of_cmatrix; ++i) {
            for(unsigned int j = 0; j < num_of_cmatrix; ++j) {
                capacity_matrix_data[i][j] = obj.getCmatValue(i, j);
            }
        }

        auto data_set = file.createDataSet<double>("capacity_matrix", HighFive::DataSpace::From(capacity_matrix_data));
        data_set.write(capacity_matrix_data);
    }

    bool loadCmatrixData(Spacecraft& obj) {
		using H5F = HighFive::File;
        const std::string file_name = "data/objects/" + obj.getName() + ".h5";

        if (MPIw::Environment::isRootNode(obj.getName())) {
            cout << "-- Cmatrix Loading from existing file: " << file_name << " --" << endl;
        }

        bool is_valid = false;
        if (Utils::isExistingFile(file_name)) {
            H5F file(file_name, H5F::ReadWrite);
            const std::string data_set_name = "capacity_matrix";

            if (file.exist(data_set_name)) {
                const auto num_of_cmatrix = obj.getCmatSize();
                auto data_set = file.getDataSet(data_set_name);
                boost::multi_array<double, 2> read_data;
                data_set.read(read_data);

                is_valid = (read_data.shape()[0] == num_of_cmatrix && read_data.shape()[1] == num_of_cmatrix);

                if (is_valid) {
                    if (MPIw::Environment::isRootNode(obj.getName())) {
                        cout << "-- [INFO] This is valid data. --" << endl;
                    }
                    for(unsigned int i = 0; i < num_of_cmatrix; ++i) {
                        for(unsigned int j = 0; j < num_of_cmatrix; ++j) {
                            obj.setCmatValue(i, j, read_data[i][j]);
                        }
                    }
                    obj.updateTotalCmatValue();
                }
            }
        }

        return is_valid;
    }

    void generateVTKHierarchicalBoxDataSet(std::shared_ptr<const RootGrid> root_grid, const int timestep, const std::string& data_type_name) {
        //! 1プロセスあたりのノード数
        const auto cx = Environment::cell_x;
        const auto cy = Environment::cell_y;
        const auto cz = Environment::cell_z;

        //! 最大グリッド幅
        const auto dx = static_cast<float>(Environment::dx);
        const auto dy = static_cast<float>(Environment::dx);
        const auto dz = static_cast<float>(Environment::dx);

        const std::string i_timestamp = (format("%08d") % timestep).str();
        const std::string f_timestamp = (format("%012.5e") % (timestep * Environment::dt)).str();

        //! VTK writer
        SimpleVTK gen;
        gen.enableExtentManagement();
        gen.changeBaseExtent(0, cx, 0, cy, 0, cz);
        gen.changeBaseOrigin(0.0, 0.0, 0.0);
        gen.changeBaseSpacing(dx, dy, dz);
        root_grid->insertAMRBlockInfo(gen, data_type_name, i_timestamp);
        std::string content = gen.getRawString();
        std::string result = MPIw::Environment::Comms["world"].gatherStringsTo(0, content);

        if (Environment::isRootNode) {
            gen.init();
            gen.beginVTK("vtkHierarchicalBoxDataSet");
            gen.setGridDescription("XYZ");
                gen.beginContent();
                    gen.addItem(result);
                gen.endContent();
            gen.endVTK();
            gen.generate(data_type_name + "_" + i_timestamp);
        }

        root_grid->plotFieldData(data_type_name, i_timestamp);
    }

    void generateXdmf(const int timestep, const std::string& data_type_name) {
        //! 1プロセスあたりのノード数
        const auto cx = Environment::cell_x;
        const auto cy = Environment::cell_y;
        const auto cz = Environment::cell_z;

        //! 最大グリッド幅
        const auto dx = static_cast<float>(Environment::dx);
        const auto dy = static_cast<float>(Environment::dx);
        const auto dz = static_cast<float>(Environment::dx);

        const std::string i_timestamp = (format("%08d") % timestep).str();
        const std::string f_timestamp = (format("%012.5e") % (timestep * Environment::dt)).str();
        //! XDMF writing
        SimpleXdmf gen;

        //! 実際のデータ用のDomain
        gen.beginDomain(data_type_name);
            gen.beginGrid("", "Collection");
                gen.beginTime("");
                gen.setValue(f_timestamp);
                gen.endTime();

                //! 全空間
                gen.begin3DStructuredGrid("Entire Space", "3DCoRectMesh", Environment::nx, Environment::ny, Environment::nz);
                    gen.add3DGeometryOrigin("", 0.0f, 0.0f, 0.0f, dx, dy, dz);
                gen.end3DStructuredGrid();

                //! 各グリッド
                for(size_t xrank = 0; xrank < Environment::proc_x; ++xrank) {
                    //! 一番上の時以外は上側のGlueノードを含める
                    const auto xsize = (xrank == Environment::proc_x - 1) ? cx : cx + 1;
                    for(size_t yrank = 0; yrank < Environment::proc_y; ++yrank) {
                        const auto ysize = (yrank == Environment::proc_y - 1) ? cy : cy + 1;
                        for(size_t zrank = 0; zrank < Environment::proc_z; ++zrank) {
                            const auto zsize = (zrank == Environment::proc_z - 1) ? cz : cz + 1;
                            const auto rank = Environment::getSpecifiedRankFromXYZRanks(xrank, yrank, zrank);
                            const std::string rankString = "proc_" + (format("%08d") % rank).str();

                            gen.begin3DStructuredGrid(rankString, "3DCoRectMesh", xsize, ysize, zsize);
                                gen.add3DGeometryOrigin("", dx * cx * xrank, dy * cy * yrank, dz * cz * zrank, dx, dy, dz);

                                if (data_type_name == "potential" || data_type_name == "rho") {
                                    gen.beginAttribute(data_type_name);
                                    gen.setCenter("Node");
                                        gen.beginDataItem();
                                        gen.setFormat("HDF");
                                        gen.setDimensions(xsize, ysize, zsize);
                                            gen.addItem(rankString + ".h5:/" + i_timestamp + "/level0/" + data_type_name);
                                        gen.endDataItem();
                                    gen.endAttribute();
                                } else if (data_type_name == "efield" || data_type_name == "bfield") {
                                    gen.beginAttribute(data_type_name, "Vector");
                                    gen.setCenter("Node");
                                        gen.beginDataItem();
                                        gen.setFormat("HDF");
                                        gen.setDimensions(xsize, ysize, zsize);
                                            gen.addItem(rankString + ".h5:/" + i_timestamp + "/level0/" + data_type_name);
                                        gen.endDataItem();
                                    gen.endAttribute();
                                } else if (data_type_name == "density") {
                                    gen.beginAttribute(data_type_name);
                                    gen.setCenter("Cell");
                                        gen.beginDataItem();
                                        gen.setFormat("HDF");
                                        gen.setDimensions(xsize - 1, ysize - 1, zsize - 1);
                                            gen.addItem(rankString + ".h5:/" + i_timestamp + "/level0/" + data_type_name);
                                        gen.endDataItem();
                                    gen.endAttribute();
                                }

                            gen.end3DStructuredGrid();
                        }
                    }
                }
            gen.endGrid();
        gen.endDomain();
        gen.generate("data/" + data_type_name + "_" + i_timestamp + ".xmf");

    }

    void putHeader(const std::string& filename, const std::string& header) {
        auto openmode = (Environment::timestep == 1) ? std::ios::out : std::ios::app;
        std::ofstream ofs(filename, openmode);

        ofs << "# " << format("%8s %16s %s") % "Timestep" % "Datatime" % header << endl; 
    }

    void putLog(const std::string& filename, const std::string& log_entry) {
        const auto datatime = Environment::getDataTime();
        std::ofstream ofs(filename, std::ios::app);
        ofs << format("%10d %16.7e %s") % Environment::timestep % datatime % log_entry << endl;
    }

    void plotEnergy(std::shared_ptr<const Grid> g, int timestep) {
        const auto receivedParticleEnergy = MPIw::Environment::Comms["world"].sum(g->getParticleEnergy(), 0);
        const auto receivedEFieldEnergy = MPIw::Environment::Comms["world"].sum(g->getEFieldEnergy(), 0);
        const auto receivedBFieldEnergy = MPIw::Environment::Comms["world"].sum(g->getBFieldEnergy(), 0);

        if(Environment::isRootNode) {
            std::string filename = "data/energy.txt";

            if(timestep == 1) {
                const auto header = (format("%16s %16s %16s %16s") % "Energy [J]" % "Particle [J]" % "EField [J]" % "BField [J]").str();
                putHeader(filename, header);
            }

            std::string entry = (format("%16.7e %16.7e %16.7e %16.7e") %
                Normalizer::unnormalizeEnergy(receivedParticleEnergy + receivedEFieldEnergy + receivedBFieldEnergy) %
                Normalizer::unnormalizeEnergy(receivedParticleEnergy) %
                Normalizer::unnormalizeEnergy(receivedEFieldEnergy) %
                Normalizer::unnormalizeEnergy(receivedBFieldEnergy)).str();

            putLog(filename, entry);
        }
    }

    void plotValidParticleNumber(std::shared_ptr<const Grid> g) {
        std::vector<size_t> total_pnums;

        size_t total_pnum = 0;
        for(int id = 0; id < Environment::num_of_particle_types; ++id) {
            total_pnums.push_back(
                MPIw::Environment::Comms["world"].sum(g->getValidParticleNumber(id))
            );
            total_pnum += total_pnums[id];
        }

        if (Environment::isRootNode) {
            std::string filename = "data/total_particle_number.txt";

            if(Environment::timestep == 1) {
                std::string format_string = "%16s";
                for(int id = 0; id < Environment::num_of_particle_types; ++id) {
                    format_string += " %16s";
                }
                auto format_base = format(format_string);
                format_base = format_base % "Total";
                for(int id = 0; id < Environment::num_of_particle_types; ++id) {
                    format_base = format_base % (Environment::getParticleType(id)->getName());
                }
                const auto header = format_base.str();
                putHeader(filename, header);
            }

            std::string format_string = "%16d";
            for(int id = 0; id < Environment::num_of_particle_types; ++id) {
                format_string += " %16d";
            }

            auto format_base = format(format_string);
            format_base = format_base % total_pnum;
            for(int id = 0; id < Environment::num_of_particle_types; ++id) {
                format_base = format_base % total_pnums[id];
            }

            std::string entry = format_base.str();
            putLog(filename, entry);
        }
    }

    void plotObjectsData(std::shared_ptr<const Grid> g) {
        for(const auto& object : g->getObjects()) {
            std::string filename = "data/object_" + object.getName() + ".txt";

            if ( MPIw::Environment::isRootNode(object.getName()) ) {
                if(Environment::timestep == 1) {
                    const auto header = object.getLogHeader();
                    putHeader(filename, header);
                }

                const std::string& entry = object.getLogEntry();
                putLog(filename, entry);
            }
        }
    }

    // x, y, zは全体に対する座標を指定する
    /*
    void plotEfieldAt(Grid const& g, int x, int y, int z, const std::string filename_header){
        const double datatime = timestep * Environment::dt;

        if(Environment::isRootNode) {
            std::string filename = "data/energy.txt";
            auto openmode = (timestep == 1) ? std::ios::out : std::ios::app;
            std::ofstream ofs(filename, openmode);

            if(timestep == 1) {
                ofs << "# " << format("%8s %15s %15s %15s %15s") % "timestep" % "time" % "Energy [J]" % "Particle [J]" % "Field [J]" << endl;
            }

            ofs << format("%10d %15.7e %15.7e %15.7e %15.7e") % timestep % datatime %
                Utils::Normalizer::unnormalizeEnergy(receivedParticleEnergy + receivedFieldEnergy) %
                Utils::Normalizer::unnormalizeEnergy(receivedParticleEnergy) %
                Utils::Normalizer::unnormalizeEnergy(receivedFieldEnergy) << endl;
        }
    }
    */


    void plotParticleVelocityDistribution(ParticleArray const& particles, const std::string filename_header) {
        plotParticleDistribution(particles, "velocity", filename_header);
    }

    void plotParticleEnergyDistribution(ParticleArray const& particles, const std::string filename_header) {
        plotParticleDistribution(particles, "energy", filename_header);
    }

    void plotParticleDistribution(ParticleArray const& particles, const std::string type, const std::string filename_header) {
        constexpr int dist_size = 200;
        const int ptypes = particles.size();
        std::vector<double> max_value(ptypes);

        //! 最大値を取得
        for(int pid = 0; pid < ptypes; ++pid) {
            for(int i = 0; i < particles[pid].size(); ++i){
                if(particles[pid][i].isValid) {
                    double val;
                    if(type == "velocity") {
                        val = particles[pid][i].getMagnitudeOfVelocity();
                    } else {
                        val = particles[pid][i].getEnergy();
                    }
                    max_value[pid] = std::max(max_value[pid], val);
                }
            }
        }

        std::vector<double> unit_value(ptypes);

        if(type == "velocity") {
            for(int pid = 0; pid < ptypes; ++pid) {
                unit_value[pid] = max_value[pid]/dist_size; //! m/s
            }
        } else {
            for(int pid = 0; pid < ptypes; ++pid) {
                unit_value[pid] = max_value[pid]/dist_size; //! eV
            }
        }

        std::vector< std::vector<int> > pdist(ptypes);

        for(int pid = 0; pid < ptypes; ++pid) {
            pdist[pid].resize(dist_size + 1);
            for(int i = 0; i < particles[pid].size(); ++i){
                if(particles[pid][i].isValid) {
                    double val;
                    if(type == "velocity") {
                        val = particles[pid][i].getMagnitudeOfVelocity();
                    } else {
                        val = particles[pid][i].getEnergy();
                    }
                    int index = floor(val/unit_value[pid]);
                    pdist[pid][index] += 1;
                }
            }
        }

        std::string filename;
        std::string header;
        filename = (format("data/%s_distribution_%04d_%04d.csv") % (filename_header + type) % MPIw::Environment::rank % Environment::timestep).str();

        if(type == "velocity") {
            header = "Velocity [km/s]";
        } else {
            header = "Energy [eV]";
        }

        std::ofstream ofs(filename, std::ios::out);

        for(int pid = 0; pid < ptypes; ++pid) {
            ofs << format("# %s") % Environment::getParticleType(pid)->getName() << endl;
            ofs << format("# %11s %11s") % header % "Ratio" << endl;

            int maxElement = *std::max_element(pdist[pid].begin(), pdist[pid].end());

            double raw_unit_value;
            if(type == "velocity") {
                raw_unit_value = Normalizer::unnormalizeVelocity(unit_value[pid]) * 1e-3; // km/s単位
            } else {
                raw_unit_value = Normalizer::unnormalizeEnergy(unit_value[pid]) / e; // eV単位
            }

            for(int i = 0; i < pdist[pid].size(); ++i){
                double ratio = static_cast<double>(pdist[pid][i]) / maxElement;
                ofs << format("  %11.4f %11.4f") % (raw_unit_value * (i+1)) % ratio << endl;
            }

            ofs << endl << endl;
        }
    }

    void print3DArray(const tdArray& data){
        const int nx = data.shape()[0];
        const int ny = data.shape()[1];
        const int nz = data.shape()[2];

        for (int k = 0; k < nz; ++k ) {
            cout << "[z:" << k << "] " << endl;

            for ( int i = 0 ; i < nx; ++i ) {
                    if(i == 0) {
                        cout << "     [x/y]";
                        for ( int j = 0 ; j < ny; ++j ) {
                            cout << "[" << j << "]";
                        }
                        cout << endl;
                    }
                    cout << "     [" << i << "]  ";
                for ( int j = 0 ; j < ny; ++j ) {
                    cout << " " << data[i][j][k] << " ";
                }
                cout << endl;
            }
            cout << endl;
        }
    }

    void outputParticlePositions(const ParticleArray& parray, std::string filename){
        std::ofstream ofs(filename);

        for(int id = 0; id < Environment::num_of_particle_types; ++id){

            ofs << "## " << Environment::getParticleType(id)->getName() << endl;

            for(int i = 0; i < parray[id].size(); ++i){
                ofs << format("%9.4f %9.4f %9.4f") % parray[id][i].x % parray[id][i].y % parray[id][i].z;
                ofs << format("%13.4e %13.4e %13.4e") % parray[id][i].vx % parray[id][i].vy % parray[id][i].vz;
                ofs << endl;
            }

            ofs << endl << endl;
        }
    }
}
