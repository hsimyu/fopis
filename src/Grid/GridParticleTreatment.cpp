#include "position.hpp"
#include "grid.hpp"
#include "particle.hpp"
#include "environment.hpp"
#include "mpiw.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

void Grid::addParticle(const Particle& p) {
    particles[p.typeId].push_back( p );
}

void Grid::addParticle(Particle&& p) {
    particles[p.typeId].push_back( std::move(p) );
}

void Grid::updateParticleVelocity(void) {
    for(auto& child : children) {
        child->updateParticleVelocity();
    }

    if (Environment::solver_type == "ES") {
        this->updateParticleVelocityES();
    } else {
        this->updateParticleVelocityEM();
    }
}

void Grid::updateParticlePosition(void) {
    for(auto& child : children) {
        child->updateParticlePosition();
    }

    if (Environment::solver_type == "ES") {
        this->updateParticlePositionES();
    } else {
        this->updateParticlePositionEM();
    }
}

void Grid::moveParticlesIntoChildren() {
    size_t moved_num = 0;
    for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
        for(auto& p : particles[pid]) {
            if (p.isValid) {
                int index = this->getChildIndexIfCovered(p.getPosition());

                if (index >= 0) {
                    this->moveParticleToChild(index, p);
                    ++moved_num;
                }
            }
        }
    }
    cout << format("[%d-%d] %d particles move to children.") % level % id % moved_num << endl;
}

void Grid::moveParticlesIntoSpecifiedChild(const int index) {
    size_t moved_num = 0;
    for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
        for(auto& p : particles[pid]) {
            if (p.isValid) {
                if (this->checkSpecifiedChildDoesCoverThisPosition(index, p.getPosition())) {
                    this->moveParticleToChild(index, p);
                    ++moved_num;
                }
            }
        }
    }
    cout << format("[%d-%d] %d particles move to children.") % level % id % moved_num << endl;
}

void Grid::moveParticleToChild(int child_index, Particle& p) {
    Particle new_particle = p; // コピー演算

    auto& child = children[child_index];
    new_particle.x = 2.0 * (new_particle.x - static_cast<double>(child->getFromIX()) + 1.0);
    new_particle.y = 2.0 * (new_particle.y - static_cast<double>(child->getFromIY()) + 1.0);
    new_particle.z = 2.0 * (new_particle.z - static_cast<double>(child->getFromIZ()) + 1.0);
    new_particle.vx /= 2.0;
    new_particle.vy /= 2.0;
    new_particle.vz /= 2.0;

    //! 親グリッド上の粒子をinvalidに
    p.makeInvalid();

    //! 子グリッド上の粒子をpush
    child->addParticle( std::move(new_particle) );
}

void Grid::removeInvalidParticles() {
    auto remove_func = [](const Particle& p) {
        return (p.isValid == 0);
    };

    for(size_t pid = 0; pid < Environment::num_of_particle_types; ++pid) {
        auto erase_iterator = std::remove_if(particles[pid].begin(), particles[pid].end(), remove_func);
        particles[pid].erase(erase_iterator, particles[pid].end());
    }
}

double Grid::getParticleEnergy(void) const {
    double res = 0.0;

    for(int pid = 0; pid < Environment::num_of_particle_types; ++pid) {
        double eachEnergy = 0.0;
        for(int i = 0; i < particles[pid].size(); ++i){
            if(particles[pid][i].isValid) {
                eachEnergy += particles[pid][i].getSquaredMagnitudeOfVelocity();
            }
        }
        res += 0.5 * Environment::getParticleType(pid)->getMass() * eachEnergy * Environment::getParticleType(pid)->getSize();
    }

    return res;
}

size_t Grid::getValidParticleNumber(const int pid) const {
    size_t count = 0;

    for(const auto& p : particles[pid]) {
        if (p.isValid) ++count;
    }

    return count;
}
