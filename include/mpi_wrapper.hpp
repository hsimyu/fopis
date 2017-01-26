#include <iostream>
#include <string>
#include <map>
#include <boost/format.hpp>
#include <mpi.h>

using std::cout;
using std::endl;
using boost::format;

namespace MPIw {
    class Communicator {
        private:
            MPI_Comm comm;

        public:
            Communicator();
            Communicator(MPI_Comm);
            //デストラクタは特に必要ない

            MPI_Comm getComm(void);
            void setComm(MPI_Comm);

            // communicate methods
            void barrier();
            double sum(double, const int);
            int sum(int, const int);
    };

    //! MPI ランク保持用のクラス
    class Environment {
        private:
            static void finalize(void);

        public:
            Environment(int, char**);
            ~Environment();

            //! MPIのランクとプロセス数はstaticに持つ
            //! int rank = MPI::Environment::rank; でアクセスする
            static int rank;
            static int numprocs;
            static int xrank, yrank, zrank;

            //! コミュニケータのリスト
            static std::map<std::string, Communicator*> Comms;

            static std::string rankStr(void);
            static void exitWithFinalize(int);
    };

};
