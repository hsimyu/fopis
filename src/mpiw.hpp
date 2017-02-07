#ifndef __TDPIC_MPIW_H_INCLUDED__
#define __TDPIC_MPIW_H_INCLUDED__
#include <mpi.h>
#include <vector>
#include <string>
#include <map>

class Particle;

namespace MPIw {
    //! MPI通信をラップするためのクラス
    class Communicator {
        //! 内部的に持つcommunicatorを使って通信する
        private:
            MPI_Comm comm;

        public:
            // Communicator
            Communicator(void) {
                comm = MPI_COMM_WORLD;
            }

            Communicator(MPI_Comm _comm) {
                comm = _comm;
            }

            MPI_Comm getComm(void);
            void setComm(MPI_Comm);

            // -- communicate methods --
            void barrier(void) {
                MPI_Barrier(comm);
            }

            // summation
            double sum(double, const int);
            int sum(int, const int);

            // send
            void send(Particle const&, const int);
            void recv(Particle&, const int);
            void sendVector(std::vector<Particle> const&, const int);
            void recvVector(std::vector<Particle>&, const int);
            void sendRecvVector(std::vector<Particle> const&, std::vector<Particle>&, const int, const int);
    };

    MPI_Datatype registerParticleType(void);
    void deregisterMpiType(MPI_Datatype);

    //! MPIランク保持用のクラス
    class Environment {
        private:
            static void finalize(void);

        public:
            Environment(int, char**);

            //! 破棄時にFinalize()する
            ~Environment() {
                finalize();
            }

            //! MPIのランクとプロセス数はstaticに持つ
            //! int rank = MPI::Environment::rank; でアクセスする
            static int rank;
            static int numprocs;
            static int xrank, yrank, zrank;
            static int adj[6];

            //! 通信用の型
            static MPI_Datatype MPI_PARTICLE;

            //! コミュニケータのリスト
            static std::map<std::string, Communicator*> Comms;

            //! Utility member functions
            static std::string rankStr(void);
            static void exitWithFinalize(int);
            static void sendRecvParticlesX(std::vector< std::vector<Particle> > const&, std::vector< std::vector<Particle> >&);
            static void sendRecvParticlesY(std::vector< std::vector<Particle> > const&, std::vector< std::vector<Particle> >&);
            static void sendRecvParticlesZ(std::vector< std::vector<Particle> > const&, std::vector< std::vector<Particle> >&);
            static void sendRecvParticles(std::vector< std::vector<Particle> > const&, std::vector< std::vector<Particle> >&, const int, const int, std::string);

    };

};

#endif
