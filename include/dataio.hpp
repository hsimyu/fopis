#ifndef __TDPIC_DATAIO_H_INCLUDED__
#define __TDPIC_DATAIO_H_INCLUDED__
#include "global.hpp"
#include "particle.hpp"
class Grid;
class Spacecraft;

namespace IO {
    // 時系列データ出力用のベース関数
    void putHeader(const std::string& filename, const std::string& header);
    void putLog(const std::string& filename, const std::string& log_entry);

    void plotEnergy(std::shared_ptr<const Grid>, int);
    void plotEfieldAt(std::shared_ptr<const Grid>, int, int, int, const std::string filename_header = "");
    void plotBfieldAt(std::shared_ptr<const Grid>, int, int, int, const std::string filename_header = "");

    void plotValidParticleNumber(std::shared_ptr<const Grid>);
    void plotParticleVelocityDistribution(ParticleArray const&, const std::string filename_header = "");
    void plotParticleEnergyDistribution(ParticleArray const&, const std::string filename_header = "");
    void plotParticleDistribution(ParticleArray const&, const std::string, const std::string);

    //! オブジェクトデータ入出力
    void plotObjectsData(std::shared_ptr<const Grid>);
    void writeCmatrixData(const Spacecraft& obj);
    bool loadCmatrixData(Spacecraft& obj);

    // 場のデータ等の出力用関数
    void writeDataInParallel(std::shared_ptr<const Grid>, const int, const std::string&);
    void generateXdmf(const int timestep, const std::string& data_type_name);

    void print3DArray(const tdArray&);
    void outputParticlePositions(const ParticleArray&, std::string filename = "data/particlePositions.csv");

}

#endif
