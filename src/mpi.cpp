#include <mpi_wrapper.hpp>

namespace MPI {
    // Environmentのstatic変数の実体
    int Environment::rank = -1;
    int Environment::numprocs = -1;

    Environment::Environment(int argc, char* argv[]) {
        MPI_Init(&argc, &argv);

        MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    }

    void Environment::finalize(void) {
        static bool finalized = false;
        if(!finalized) {
            MPI_Finalize();
            finalized = true;
        }
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
        return (format("[P%04d] ") % rank).str();
    }

    // Communicator
    Communicator::Communicator(void) {
        comm = MPI_COMM_WORLD;
    }

    Communicator::Communicator(MPI_Comm _comm) {
        comm = _comm;
    }
}
