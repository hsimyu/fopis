#include <tdpic.h>
#include <cstddef>

namespace MPIw {
    // Environmentのstatic変数の実体
    int Environment::rank = -1;
    int Environment::numprocs = -1;
    int Environment::xrank = -1;
    int Environment::yrank = -1;
    int Environment::zrank = -1;
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

    // 破棄時にFinalize()する
    Environment::~Environment() {
        finalize();
    }

    void Environment::exitWithFinalize(int code) {
        finalize();
        exit(code);
    }

    std::string Environment::rankStr(void) {
        return (format("[RANK P%04d] ") % rank).str();
    }

    MPI_Datatype registerParticleType() {
        const size_t num_members = 7;
        int lengths[num_members] = {1, 1, 1, 1, 1, 1, 1};

        MPI_Aint offsets[num_members] = {
            offsetof(Particle, typeId),
            offsetof(Particle, x),
            offsetof(Particle, y),
            offsetof(Particle, z),
            offsetof(Particle, vx),
            offsetof(Particle, vy),
            offsetof(Particle, vz)
        };
        MPI_Datatype types[num_members] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};

        MPI_Datatype newType;
        MPI_Type_create_struct(num_members, lengths, offsets, types, &newType);
        MPI_Type_commit(&newType);

        return newType;
    }

    void deregisterMpiType(MPI_Datatype type){
        MPI_Type_free(&type);
    }

    // Communicator
    Communicator::Communicator(void) {
        comm = MPI_COMM_WORLD;
    }

    Communicator::Communicator(MPI_Comm _comm) {
        comm = _comm;
    }

    void Communicator::send(Particle const& p, const int dest, const int tag) {
        MPI_Send(&p, 1, Environment::MPI_PARTICLE, dest, tag, comm);
    }

    void Communicator::recv(Particle& p, const int src, const int tag) {
        MPI_Status s;
        MPI_Recv(&p, 1, Environment::MPI_PARTICLE, src, tag, comm, &s);
    }

    void Communicator::sendVector(std::vector<Particle> const& parray, const int dest, const int tag) {
        unsigned int len = parray.size();
        MPI_Send(&len, 1, MPI_UNSIGNED, dest, tag, comm);

        if(len != 0) {
            MPI_Send(parray.data(), len, Environment::MPI_PARTICLE, dest, tag, comm);
        }
    }

    void Communicator::recvVector(std::vector<Particle>& parray, const int src, const int tag) {
        MPI_Status s;
        unsigned int len;
        MPI_Recv(&len, 1, MPI_UNSIGNED, src, tag, comm, &s);

        if(len != 0) {
            parray.resize(len);
            MPI_Recv(parray.data(), len, Environment::MPI_PARTICLE, src, tag, comm, &s);
        } else {
            parray.clear();
            parray.shrink_to_fit();
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

    void Communicator::barrier(void) {
        MPI_Barrier(comm);
    }
}
