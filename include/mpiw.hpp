#ifndef __TDPIC_MPIW_H_INCLUDED__
#define __TDPIC_MPIW_H_INCLUDED__
#include <mpi.h>
#include <vector>
#include <string>
#include <map>

class Particle;

namespace MPIw {
    enum TAG {
        SEND_PARTICLE = 1,
        SEND_PARTICLE_LENGTH,
        RECV_PARTICLE,
        RECV_PARTICLE_LENGTH,
        SENDRECV_PARTICLE,
        SENDRECV_PARTICLE_LENGTH,
        SENDRECV_FIELD,
        SENDRECV_SCALAR_UPPER,
        SENDRECV_SCALAR_LOWER,
        PARTICIPATE_NEW_COMM,
        GATHER_STRINGS,
    };

    //! MPI通信をラップするためのクラス
    class Communicator {
        //! 内部的に持つcommunicatorを使って通信する
        private:
            MPI_Comm comm;
            int root;

        public:
            // CommunicatorのデフォルトCommはMPI_COMM_WORLDとする
            Communicator(void) {
                comm = MPI_COMM_WORLD;
                root = 0;
            }

            Communicator(MPI_Comm _comm);

            MPI_Comm getComm(void) const { return comm; }
            int getRoot(void) const { return root; }

            // -- communicate methods --
            void barrier(void) {
                MPI_Barrier(comm);
            }

            // summation only on target
            double sum(double, const int);
            int sum(int, const int);
            size_t sum(size_t, const int);

            // all reducer
            double sum(double);
            int    sum(int);
            size_t sum(size_t);
            double max(double);
            double min(double);
            int    min(int);

            std::vector<float> sum(std::vector<float>&);
            boost::multi_array<double, 2> sum(boost::multi_array<double, 2>&);

            // treats strings
            std::string gatherStringsTo(int target_rank, std::string& content) const;

            // send
            void send(Particle const&, const int);
            void recv(Particle&, const int);
            void sendVector(std::vector<Particle> const&, const int);
            void recvVector(std::vector<Particle>&, const int);
            void sendRecvVector(std::vector<Particle> const&, std::vector<Particle>&, const int, const int);
            void sendRecvFieldX(tdArray&, const int, const int);
            void sendRecvFieldY(tdArray&, const int, const int);
            void sendRecvFieldZ(tdArray&, const int, const int);
            void sendRecvScalarX(tdArray&, const int, const int, const int);
            void sendRecvScalarY(tdArray&, const int, const int, const int);
            void sendRecvScalarZ(tdArray&, const int, const int, const int);
    };

    MPI_Datatype registerParticleType(void);
    void deregisterMpiType(MPI_Datatype);

    //! MPIランク保持用のクラス
    class Environment {
        public:
            Environment(int, char**);

            ~Environment() {
                finalize();
            }

            //! 正常終了時
            static void finalize(void);

            //! 例外時などの終了処理のために呼び出す
            static void abort(const int errorcode);

            //! MPIランクを接頭辞として出力する時のための関数
            static std::string rankStr(void) {
                return (format("[RANK P%04d] ") % rank).str();
            }

            static int getRootNode(const std::string& comm_name) {
                if (Comms.count(comm_name) == 0) return -1;

                return Comms[comm_name].getRoot();
            }

            static bool isRootNode(const std::string& comm_name) {
                if (Comms.count(comm_name) == 0) return false;

                return rank == Comms[comm_name].getRoot();
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
            static std::map<std::string, Communicator> Comms;
            static void addNewComm(const std::string& new_comm_name, const MPI_Comm new_comm);
            static void makeNewComm(const std::string& new_comm_name, const bool is_not_empty_comm);

            //! MPI通信をラップするためのメンバ関数群
            static void sendRecvParticlesX(std::vector< std::vector<Particle> > const&, std::vector< std::vector<Particle> >&);
            static void sendRecvParticlesY(std::vector< std::vector<Particle> > const&, std::vector< std::vector<Particle> >&);
            static void sendRecvParticlesZ(std::vector< std::vector<Particle> > const&, std::vector< std::vector<Particle> >&);
            static void sendRecvParticles(std::vector< std::vector<Particle> > const&, std::vector< std::vector<Particle> >&, const int, const int, std::string);

            static void sendRecvField(tdArray&);
            static void sendRecvNodeScalar(tdArray&);
            static void sendRecvCellScalar(tdArray&);
    };

};

#endif
