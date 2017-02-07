#include "global.hpp"
#include "mpiw.hpp"
#include "environment.hpp"
#include "particle.hpp"
#include <cstddef>

namespace MPIw {
    enum TAG {
        SEND_PARTICLE = 1,
        SEND_PARTICLE_LENGTH,
        RECV_PARTICLE,
        RECV_PARTICLE_LENGTH,
        SENDRECV_PARTICLE,
        SENDRECV_PARTICLE_LENGTH,
    };

    // Environmentのstatic変数の実体
    int Environment::rank = -1;
    int Environment::numprocs = -1;
    int Environment::xrank = -1;
    int Environment::yrank = -1;
    int Environment::zrank = -1;

    // 各方向への隣接プロセスランク
    int Environment::adj[6];

    std::map<std::string, Communicator*> Environment::Comms;
    MPI_Datatype Environment::MPI_PARTICLE;

    Environment::Environment(int argc, char* argv[]) {
#ifndef BUILD_TEST
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        Comms["world"] = new Communicator();
#else
        rank = 0;
        numprocs = 1;
        xrank = 0;
        yrank = 0;
        zrank = 0;
#endif
        MPI_PARTICLE = registerParticleType();
    }

    void Environment::finalize(void) {
#ifndef BUILD_TEST
        static bool finalized = false;
        if(!finalized) {
            MPI_Finalize();
            finalized = true;
        }
#endif
    }

    void Environment::sendRecvParticles(std::vector< std::vector<Particle> > const& pbuff, std::vector< std::vector<Particle> >& pbuffRecv, const int prev, const int next, std::string commName){
        Comms[commName]->sendRecvVector(pbuff[prev], pbuffRecv[next], adj[prev], adj[next]);
        Comms[commName]->sendRecvVector(pbuff[next], pbuffRecv[prev], adj[next], adj[prev]);
    }

    void Environment::sendRecvParticlesX(std::vector< std::vector<Particle> > const& pbuff, std::vector< std::vector<Particle> >& pbuffRecv){
        //! X-axis
        int prev = 0; int next = 1;
        sendRecvParticles(pbuff, pbuffRecv, prev, next, "world");
    }

    void Environment::sendRecvParticlesY(std::vector< std::vector<Particle> > const& pbuff, std::vector< std::vector<Particle> >& pbuffRecv){
        //! Y-axis
        int prev = 2; int next = 3;
        sendRecvParticles(pbuff, pbuffRecv, prev, next, "world");
    }

    void Environment::sendRecvParticlesZ(std::vector< std::vector<Particle> > const& pbuff, std::vector< std::vector<Particle> >& pbuffRecv){
        //! Z-axis
        int prev = 4; int next = 5;
        sendRecvParticles(pbuff, pbuffRecv, prev, next, "world");
    }

    MPI_Datatype registerParticleType() {
        const size_t num_members = 8;
        int lengths[num_members] = {1, 1, 1, 1, 1, 1, 1, 1};

        MPI_Aint offsets[num_members] = {
            offsetof(Particle, typeId),
            offsetof(Particle, isValid),
            offsetof(Particle, x),
            offsetof(Particle, y),
            offsetof(Particle, z),
            offsetof(Particle, vx),
            offsetof(Particle, vy),
            offsetof(Particle, vz)
        };
        MPI_Datatype types[num_members] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};

        MPI_Datatype newType;
        MPI_Type_create_struct(num_members, lengths, offsets, types, &newType);
        MPI_Type_commit(&newType);

        return newType;
    }

    void deregisterMpiType(MPI_Datatype type){
        MPI_Type_free(&type);
    }

    void Communicator::send(Particle const& p, const int dest) {
        MPI_Send(&p, 1, Environment::MPI_PARTICLE, dest, TAG::SEND_PARTICLE, comm);
    }

    void Communicator::recv(Particle& p, const int src) {
        MPI_Status s;
        MPI_Recv(&p, 1, Environment::MPI_PARTICLE, src, TAG::RECV_PARTICLE, comm, &s);
    }

    void Communicator::sendVector(std::vector<Particle> const& parray, const int dest) {
        unsigned int len = parray.size();
        MPI_Send(&len, 1, MPI_UNSIGNED, dest, TAG::SEND_PARTICLE_LENGTH, comm);

        if(len != 0) {
            MPI_Send(parray.data(), len, Environment::MPI_PARTICLE, dest, TAG::SEND_PARTICLE, comm);
        }
    }

    void Communicator::recvVector(std::vector<Particle>& parray, const int src) {
        MPI_Status s;
        unsigned int len;
        MPI_Recv(&len, 1, MPI_UNSIGNED, src, TAG::RECV_PARTICLE_LENGTH, comm, &s);

        if(len != 0) {
            parray.resize(len);
            MPI_Recv(parray.data(), len, Environment::MPI_PARTICLE, src, TAG::RECV_PARTICLE, comm, &s);
        } else {
            parray.clear();
            parray.shrink_to_fit();
        }
    }

    void Communicator::sendRecvVector(std::vector<Particle> const& sendArray, std::vector<Particle>& recvArray, const int dest, const int src) {
        unsigned int sendlen = sendArray.size();
        unsigned int recvlen;
        MPI_Status s;
        MPI_Sendrecv(&sendlen, 1, MPI_UNSIGNED, dest, TAG::SENDRECV_PARTICLE_LENGTH, &recvlen, 1, MPI_UNSIGNED, src, TAG::SENDRECV_PARTICLE_LENGTH, comm, &s);

        if(sendlen != 0 && recvlen != 0) {
            recvArray.resize(recvlen);
            MPI_Sendrecv(sendArray.data(), sendlen, Environment::MPI_PARTICLE, dest, TAG::SENDRECV_PARTICLE, recvArray.data(), recvlen, Environment::MPI_PARTICLE, src, TAG::SENDRECV_PARTICLE, comm, &s);
        } else if(sendlen != 0) {
            MPI_Send(sendArray.data(), sendlen, Environment::MPI_PARTICLE, dest, TAG::SENDRECV_PARTICLE, comm);
            recvArray.clear();
            recvArray.shrink_to_fit();
        } else if(recvlen != 0) {
            recvArray.resize(recvlen);
            MPI_Recv(recvArray.data(), recvlen, Environment::MPI_PARTICLE, src, TAG::SENDRECV_PARTICLE, comm, &s);
        } else {
            recvArray.clear();
            recvArray.shrink_to_fit();
        }
    }

    double Communicator::sum(double value, const int target_rank) {
        double res = 0.0;
        MPI_Reduce(&value, &res, 1, MPI_DOUBLE, MPI_SUM, target_rank, comm);
        return res;
    }

    int Communicator::sum(int value, const int target_rank) {
        int res = 0;
        MPI_Reduce(&value, &res, 1, MPI_INT, MPI_SUM, target_rank, comm);
        return res;
    }
}
