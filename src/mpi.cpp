#include <mpi_wrapper.hpp>

namespace MPIw {
    // Environmentのstatic変数の実体
    int Environment::rank = -1;
    int Environment::numprocs = -1;
    int Environment::xrank = -1;
    int Environment::yrank = -1;
    int Environment::zrank = -1;
    std::map<std::string, Communicator*> Environment::Comms;

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

    // Communicator
    Communicator::Communicator(void) {
        comm = MPI_COMM_WORLD;
    }

    Communicator::Communicator(MPI_Comm _comm) {
        comm = _comm;
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
