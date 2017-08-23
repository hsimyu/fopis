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
        SENDRECV_FIELD,
        PARTICIPATE_NEW_COMM,
    };

    // Environmentのstatic変数の実体
    int Environment::rank = -1;
    int Environment::numprocs = -1;
    int Environment::xrank = -1;
    int Environment::yrank = -1;
    int Environment::zrank = -1;

    // 各方向への隣接プロセスランク
    int Environment::adj[6];

    std::map<std::string, Communicator> Environment::Comms;
    MPI_Datatype Environment::MPI_PARTICLE;

    Environment::Environment(int argc, char* argv[]) {
#ifndef BUILD_TEST
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        addNewComm("world", MPI_COMM_WORLD);
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

    void Environment::abort(const int code) {
        MPI_Abort(MPI_COMM_WORLD, code);
    }

    void Environment::addNewComm(const std::string& new_comm_name, const MPI_Comm new_comm) {
        Comms.emplace(std::piecewise_construct, std::make_tuple(new_comm_name), std::make_tuple(new_comm));
    }

    void Environment::makeNewComm(const std::string& new_comm_name, const bool is_not_empty_comm) {
        auto source_comm = Comms["world"].getComm();
        MPI_Comm new_comm;
        int color = (is_not_empty_comm ? 0 : MPI_UNDEFINED);
        MPI_Comm_split(source_comm, color, rank, &new_comm);
        if (is_not_empty_comm) addNewComm(new_comm_name, new_comm);
    }

    void Environment::sendRecvField(tdArray& x3D){
        int prev, next;

        // 対応する方向の proc 数が 1 かつ周期境界でない場合には通信しなくてよい
        if( (::Environment::proc_x > 1) || ::Environment::isNotBoundary(AXIS::x, AXIS_SIDE::low) ) {
            prev = 0; next = 1;
            Comms["world"].sendRecvFieldX(x3D, adj[prev], adj[next]);
        }

        if( (::Environment::proc_y > 1) || ::Environment::isNotBoundary(AXIS::y, AXIS_SIDE::low) ) {
            prev = 2; next = 3;
            Comms["world"].sendRecvFieldY(x3D, adj[prev], adj[next]);
        }

        if( (::Environment::proc_z > 1) || ::Environment::isNotBoundary(AXIS::z, AXIS_SIDE::low) ) {
            prev = 4; next = 5;
            Comms["world"].sendRecvFieldZ(x3D, adj[prev], adj[next]);
        }
    }

    //! -- Particle Communication --
    void Environment::sendRecvParticles(std::vector< std::vector<Particle> > const& pbuff, std::vector< std::vector<Particle> >& pbuffRecv, const int prev, const int next, std::string commName){
        Comms[commName].sendRecvVector(pbuff[prev], pbuffRecv[next], adj[prev], adj[next]);
        Comms[commName].sendRecvVector(pbuff[next], pbuffRecv[prev], adj[next], adj[prev]);
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
        const unsigned int sendlen = sendArray.size();
        unsigned int recvlen;

        if ( (src == Environment::rank) && (dest == Environment::rank) ) {
            //! 自分との通信
            if (sendlen != 0) {
                recvArray = sendArray; // 単にコピーする
            } else {
                recvArray.clear();
                recvArray.shrink_to_fit();
            }
        } else {
            MPI_Status s;
            MPI_Sendrecv(&sendlen, 1, MPI_UNSIGNED, dest, TAG::SENDRECV_PARTICLE_LENGTH, &recvlen, 1, MPI_UNSIGNED, src, TAG::SENDRECV_PARTICLE_LENGTH, comm, &s);

            if(sendlen != 0 && recvlen != 0) {
                //! ここでParticleのコンストラクタが呼ばれてしまうことに注意
                recvArray.resize(recvlen);
                MPI_Sendrecv(sendArray.data(), sendlen, Environment::MPI_PARTICLE, dest, TAG::SENDRECV_PARTICLE, recvArray.data(), recvlen, Environment::MPI_PARTICLE, src, TAG::SENDRECV_PARTICLE, comm, &s);
            } else if(sendlen != 0) {
                MPI_Send(sendArray.data(), sendlen, Environment::MPI_PARTICLE, dest, TAG::SENDRECV_PARTICLE, comm);
                recvArray.clear();
                recvArray.shrink_to_fit();
            } else if(recvlen != 0) {
                //! ここでParticleのコンストラクタが呼ばれてしまうことに注意
                recvArray.resize(recvlen);
                MPI_Recv(recvArray.data(), recvlen, Environment::MPI_PARTICLE, src, TAG::SENDRECV_PARTICLE, comm, &s);
            } else {
                recvArray.clear();
                recvArray.shrink_to_fit();
            }
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

    double Communicator::sum(double value) {
        double res = 0.0;
        MPI_Allreduce(&value, &res, 1, MPI_DOUBLE, MPI_SUM, comm);
        return res;
    }

    int Communicator::sum(int value) {
        int res = 0;
        MPI_Allreduce(&value, &res, 1, MPI_INT, MPI_SUM, comm);
        return res;
    }

    double Communicator::max(double value) {
        double res = 0.0;
        MPI_Allreduce(&value, &res, 1, MPI_DOUBLE, MPI_MAX, comm);
        return res;
    }

    double Communicator::min(double value) {
        double res = 0.0;
        MPI_Allreduce(&value, &res, 1, MPI_DOUBLE, MPI_MIN, comm);
        return res;
    }

    //! @note: prev == nextの場合はわざわざ通信しなくてよいので効率を上げられる
    void Communicator::sendRecvFieldX(tdArray& tdValue, const int prev, const int next) {
        MPI_Status s;
        boost::multi_array<double, 2> buff(boost::extents[ tdValue.shape()[1] ][ tdValue.shape()[2] ]);

        //! 下側の値をバッファへ詰める
        for(int j = 0; j < tdValue.shape()[1]; ++j) {
            for(int k = 0; k < tdValue.shape()[2]; ++k) {
                buff[j][k] = tdValue[1][j][k];
            }
        }

        if ( (prev == Environment::rank) && (next == Environment::rank) ) {
            //! 自分自身への通信
            //! std::swapのが速いかも
            for(int j = 0; j < tdValue.shape()[1]; ++j) {
                for(int k = 0; k < tdValue.shape()[2]; ++k) {
                    tdValue[0][j][k] = tdValue[tdValue.shape()[0] - 2][j][k];
                    tdValue[tdValue.shape()[0] - 1][j][k] = buff[j][k];
                }
            }
        } else {
            int length = tdValue.shape()[1] * tdValue.shape()[2];

            MPI_Sendrecv_replace(buff.data(), length, MPI_DOUBLE, prev, TAG::SENDRECV_FIELD, next, TAG::SENDRECV_FIELD, comm, &s);

            for(int j = 0; j < tdValue.shape()[1]; ++j) {
                for(int k = 0; k < tdValue.shape()[2]; ++k) {
                    tdValue[tdValue.shape()[0] - 1][j][k] = buff[j][k];

                    // 反対側の値を buff に詰める
                    buff[j][k] = tdValue[tdValue.shape()[0] - 2][j][k];
                }
            }

            MPI_Sendrecv_replace(buff.data(), length, MPI_DOUBLE, next, TAG::SENDRECV_FIELD, prev, TAG::SENDRECV_FIELD, comm, &s);

            for(int j = 0; j < tdValue.shape()[1]; ++j) {
                for(int k = 0; k < tdValue.shape()[2]; ++k) {
                    tdValue[0][j][k] = buff[j][k];
                }
            }
        }
    }

    void Communicator::sendRecvFieldY(tdArray& tdValue, const int prev, const int next) {
        MPI_Status s;
        boost::multi_array<double, 2> buff(boost::extents[ tdValue.shape()[0] ][ tdValue.shape()[2] ]);

        for(int i = 0; i < tdValue.shape()[0]; ++i) {
            for(int k = 0; k < tdValue.shape()[2]; ++k) {
                buff[i][k] = tdValue[i][1][k];
            }
        }

        if ( (prev == Environment::rank) && (next == Environment::rank) ) {
            for(int i = 0; i < tdValue.shape()[0]; ++i) {
                for(int k = 0; k < tdValue.shape()[2]; ++k) {
                    tdValue[i][0][k] = tdValue[i][tdValue.shape()[1] - 2][k];
                    tdValue[i][tdValue.shape()[1] - 1][k] = buff[i][k];
                }
            }
        } else {
            int length = tdValue.shape()[0] * tdValue.shape()[2];

            MPI_Sendrecv_replace(buff.data(), length, MPI_DOUBLE, prev, TAG::SENDRECV_FIELD, next, TAG::SENDRECV_FIELD, comm, &s);

            for(int i = 0; i < tdValue.shape()[0]; ++i) {
                for(int k = 0; k < tdValue.shape()[2]; ++k) {
                    tdValue[i][tdValue.shape()[1] - 1][k] = buff[i][k];
                    buff[i][k] = tdValue[i][tdValue.shape()[1] - 2][k];
                }
            }

            MPI_Sendrecv_replace(buff.data(), length, MPI_DOUBLE, next, TAG::SENDRECV_FIELD, prev, TAG::SENDRECV_FIELD, comm, &s);

            for(int i = 0; i < tdValue.shape()[0]; ++i) {
                for(int k = 0; k < tdValue.shape()[2]; ++k) {
                    tdValue[i][0][k] = buff[i][k];
                }
            }
        }
    }

    void Communicator::sendRecvFieldZ(tdArray& tdValue, const int prev, const int next) {
        MPI_Status s;
        boost::multi_array<double, 2> buff(boost::extents[ tdValue.shape()[0] ][ tdValue.shape()[1] ]);

        for(int i = 0; i < tdValue.shape()[0]; ++i) {
            for(int j = 0; j < tdValue.shape()[1]; ++j) {
                buff[i][j] = tdValue[i][j][1];
            }
        }

        if ( (prev == Environment::rank) && (next == Environment::rank) ) {
            for(int i = 0; i < tdValue.shape()[0]; ++i) {
                for(int j = 0; j < tdValue.shape()[1]; ++j) {
                    tdValue[i][j][0] = tdValue[i][j][tdValue.shape()[2] - 2];
                    tdValue[i][j][tdValue.shape()[2] - 1] = buff[i][j];
                }
            }
        } else {

            int length = tdValue.shape()[0] * tdValue.shape()[1];

            MPI_Sendrecv_replace(buff.data(), length, MPI_DOUBLE, prev, TAG::SENDRECV_FIELD, next, TAG::SENDRECV_FIELD, comm, &s);

            for(int i = 0; i < tdValue.shape()[0]; ++i) {
                for(int j = 0; j < tdValue.shape()[1]; ++j) {
                    tdValue[i][j][tdValue.shape()[2] - 1] = buff[i][j];
                    buff[i][j] = tdValue[i][j][tdValue.shape()[2] - 2];
                }
            }

            MPI_Sendrecv_replace(buff.data(), length, MPI_DOUBLE, next, TAG::SENDRECV_FIELD, prev, TAG::SENDRECV_FIELD, comm, &s);

            for(int i = 0; i < tdValue.shape()[0]; ++i) {
                for(int j = 0; j < tdValue.shape()[1]; ++j) {
                    tdValue[i][j][0] = buff[i][j];
                }
            }
        }
    }
}
