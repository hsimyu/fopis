#ifndef __TDPIC_DATAIO_H_INCLUDED__
#define __TDPIC_DATAIO_H_INCLUDED__
#include "global.hpp"
#include "particle.hpp"
class Grid;

namespace IO {
    void plotEnergy(Grid const&, int);
    void plotParticleVelocityDistribution(Grid const&);
    void plotParticleEnergyDistribution(Grid const&);
    void plotParticleDistribution(Grid const&, const std::string);

    void writeDataInParallel(Grid*, int, std::string);
    void print3DArray(const tdArray&);
    void outputParticlePositions(const ParticleArray&, std::string filename = "data/particlePositions.csv");
}

#endif
