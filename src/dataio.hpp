#ifndef __TDPIC_DATAIO_H_INCLUDED__
#define __TDPIC_DATAIO_H_INCLUDED__
#include "global.hpp"
#include "particle.hpp"
class Grid;

namespace IO {
    void plotEnergy(Grid const&, int);
    void plotParticleVelocityDistribution(ParticleArray const&, const std::string filename_header = "");
    void plotParticleEnergyDistribution(ParticleArray const&, const std::string filename_header = "");
    void plotParticleDistribution(ParticleArray const&, const std::string, const std::string);

    void writeDataInParallel(Grid&, int, std::string);
    void print3DArray(const tdArray&);
    void outputParticlePositions(const ParticleArray&, std::string filename = "data/particlePositions.csv");
}

#endif
