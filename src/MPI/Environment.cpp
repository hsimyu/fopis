#include "global.hpp"
#include "mpiw.hpp"
#include "environment.hpp"
#include "particle.hpp"

namespace MPIw {
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
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        addNewComm("world", MPI_COMM_WORLD);
        MPI_PARTICLE = registerParticleType();
    }

    void Environment::finalize(void) {
        static bool finalized = false;

        if(!finalized) {
            MPI_Finalize();
            finalized = true;
        }
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

    void Environment::sendRecvNodeScalar(tdArray& x3D){
        int prev, next;
        constexpr int node_glue_size = 2;

        // 対応する方向の proc 数が 1 かつ周期境界でない場合には通信しなくてよい
        if( (::Environment::proc_x > 1) || ::Environment::isNotBoundary(AXIS::x, AXIS_SIDE::low) ) {
            prev = 0; next = 1;
            Comms["world"].sendRecvScalarX(x3D, adj[prev], adj[next], node_glue_size);
        }

        if( (::Environment::proc_y > 1) || ::Environment::isNotBoundary(AXIS::y, AXIS_SIDE::low) ) {
            prev = 2; next = 3;
            Comms["world"].sendRecvScalarY(x3D, adj[prev], adj[next], node_glue_size);
        }

        if( (::Environment::proc_z > 1) || ::Environment::isNotBoundary(AXIS::z, AXIS_SIDE::low) ) {
            prev = 4; next = 5;
            Comms["world"].sendRecvScalarZ(x3D, adj[prev], adj[next], node_glue_size);
        }
    }

    void Environment::sendRecvEdgeScalar(tdArray& x3D, AXIS axis){
        int prev, next;
        constexpr int node_glue_size = 2;
        constexpr int cell_glue_size = 1;

        //! Edgeの場合、どの向きのベクトルかによってGlueサイズが変わる
        if (axis == AXIS::x) {
            //! x方向のみGlueサイズ = 1
            if( (::Environment::proc_x > 1) || ::Environment::isNotBoundary(AXIS::x, AXIS_SIDE::low) ) {
                prev = 0; next = 1;
                Comms["world"].sendRecvScalarX(x3D, adj[prev], adj[next], node_glue_size);
            }

            //! その他方向はGlueサイズ = 2
            if( (::Environment::proc_y > 1) || ::Environment::isNotBoundary(AXIS::y, AXIS_SIDE::low) ) {
                prev = 2; next = 3;
                Comms["world"].sendRecvScalarY(x3D, adj[prev], adj[next], cell_glue_size);
            }

            if( (::Environment::proc_z > 1) || ::Environment::isNotBoundary(AXIS::z, AXIS_SIDE::low) ) {
                prev = 4; next = 5;
                Comms["world"].sendRecvScalarZ(x3D, adj[prev], adj[next], cell_glue_size);
            }
        } else if (axis == AXIS::y) {
            if( (::Environment::proc_x > 1) || ::Environment::isNotBoundary(AXIS::x, AXIS_SIDE::low) ) {
                prev = 0; next = 1;
                Comms["world"].sendRecvScalarX(x3D, adj[prev], adj[next], cell_glue_size);
            }

            //! y方向のみGlueサイズ = 1
            if( (::Environment::proc_y > 1) || ::Environment::isNotBoundary(AXIS::y, AXIS_SIDE::low) ) {
                prev = 2; next = 3;
                Comms["world"].sendRecvScalarY(x3D, adj[prev], adj[next], node_glue_size);
            }

            if( (::Environment::proc_z > 1) || ::Environment::isNotBoundary(AXIS::z, AXIS_SIDE::low) ) {
                prev = 4; next = 5;
                Comms["world"].sendRecvScalarZ(x3D, adj[prev], adj[next], cell_glue_size);
            }
        } else {
            if( (::Environment::proc_x > 1) || ::Environment::isNotBoundary(AXIS::x, AXIS_SIDE::low) ) {
                prev = 0; next = 1;
                Comms["world"].sendRecvScalarX(x3D, adj[prev], adj[next], cell_glue_size);
            }

            //! y方向のみGlueサイズ = 1
            if( (::Environment::proc_y > 1) || ::Environment::isNotBoundary(AXIS::y, AXIS_SIDE::low) ) {
                prev = 2; next = 3;
                Comms["world"].sendRecvScalarY(x3D, adj[prev], adj[next], cell_glue_size);
            }

            if( (::Environment::proc_z > 1) || ::Environment::isNotBoundary(AXIS::z, AXIS_SIDE::low) ) {
                prev = 4; next = 5;
                Comms["world"].sendRecvScalarZ(x3D, adj[prev], adj[next], node_glue_size);
            }
        }
    }

    void Environment::sendRecvCellScalar(tdArray& x3D){
        int prev, next;
        constexpr int cell_glue_size = 1;

        // 対応する方向の proc 数が 1 かつ周期境界でない場合には通信しなくてよい
        if( (::Environment::proc_x > 1) || ::Environment::isNotBoundary(AXIS::x, AXIS_SIDE::low) ) {
            prev = 0; next = 1;
            Comms["world"].sendRecvScalarX(x3D, adj[prev], adj[next], cell_glue_size);
        }

        if( (::Environment::proc_y > 1) || ::Environment::isNotBoundary(AXIS::y, AXIS_SIDE::low) ) {
            prev = 2; next = 3;
            Comms["world"].sendRecvScalarY(x3D, adj[prev], adj[next], cell_glue_size);
        }

        if( (::Environment::proc_z > 1) || ::Environment::isNotBoundary(AXIS::z, AXIS_SIDE::low) ) {
            prev = 4; next = 5;
            Comms["world"].sendRecvScalarZ(x3D, adj[prev], adj[next], cell_glue_size);
        }
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

    // for PSOR
    void Environment::sendRecvPartialPhi(tdArray& phi, const size_t i_begin, const size_t i_end, const size_t j_begin, const size_t j_end, const size_t k_begin, const size_t k_end){
        if (i_begin == i_end) {
            if (i_begin == 1) {
                //! 下側へ送る
                if( (::Environment::proc_x > 1) || ::Environment::isNotBoundary(AXIS::x, AXIS_SIDE::low) ) {
                    Comms["world"].sendRecvPartialFieldX(phi, i_begin, j_begin, j_end, k_begin, k_end, adj[0], adj[1]);
                }
            } else {
                //! 上側へ送る
                if( (::Environment::proc_x > 1) || ::Environment::isNotBoundary(AXIS::x, AXIS_SIDE::up) ) {
                    Comms["world"].sendRecvPartialFieldX(phi, i_begin, j_begin, j_end, k_begin, k_end, adj[1], adj[0]);
                }
            }
        } else if (j_begin == j_end) {
            if (j_begin == 1) {
                if( (::Environment::proc_y > 1) || ::Environment::isNotBoundary(AXIS::y, AXIS_SIDE::low) ) {
                    Comms["world"].sendRecvPartialFieldY(phi, i_begin, i_end, j_begin, k_begin, k_end, adj[2], adj[3]);
                }
            } else {
                if( (::Environment::proc_y > 1) || ::Environment::isNotBoundary(AXIS::y, AXIS_SIDE::up) ) {
                    Comms["world"].sendRecvPartialFieldY(phi, i_begin, i_end, j_begin, k_begin, k_end, adj[3], adj[2]);
                }
            }
        } else if (k_begin == k_end) {
            if (k_begin == 1) {
                if( (::Environment::proc_z > 1) || ::Environment::isNotBoundary(AXIS::z, AXIS_SIDE::low) ) {
                    Comms["world"].sendRecvPartialFieldZ(phi, i_begin, i_end, j_begin, j_end, k_begin, adj[4], adj[5]);
                }
            } else {
                if( (::Environment::proc_z > 1) || ::Environment::isNotBoundary(AXIS::z, AXIS_SIDE::up) ) {
                    Comms["world"].sendRecvPartialFieldZ(phi, i_begin, i_end, j_begin, j_end, k_begin, adj[5], adj[4]);
                }
            }
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

}
