#include "global.hpp"
#include "mpiw.hpp"
#include "particle.hpp"
#include <cstddef>

namespace MPIw {
    // Communicator クラス
    Communicator::Communicator(MPI_Comm _comm) {
        comm = _comm;
        root = this->min(Environment::rank);
    }

    void Communicator::send(Particle const& p, const int dest) {
        MPI_Send(&p, 1, Environment::MPI_PARTICLE, dest, TAG::SEND_PARTICLE, comm);
    }

    void Communicator::recv(Particle& p, const int src) {
        MPI_Status s;
        MPI_Recv(&p, 1, Environment::MPI_PARTICLE, src, TAG::RECV_PARTICLE, comm, &s);
    }

    void Communicator::sendVector(std::vector<Particle> const& parray, const int dest) {
        int len = static_cast<int>(parray.size());
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
        const int sendlen = static_cast<int>(sendArray.size());
        int recvlen;

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

//! size_t通信用の型判定
#if SIZE_MAX == UCHAR_MAX
    #define MY_MPI_SIZE_T MPI_UNSIGNED_CHAR
 #elif SIZE_MAX == USHRT_MAX
    #define MY_MPI_SIZE_T MPI_UNSIGNED_SHORT
 #elif SIZE_MAX == UINT_MAX
    #define MY_MPI_SIZE_T MPI_UNSIGNED
 #elif SIZE_MAX == ULONG_MAX
    #define MY_MPI_SIZE_T MPI_UNSIGNED_LONG
 #elif SIZE_MAX == ULLONG_MAX
    #define MY_MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
 #endif

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

    size_t Communicator::sum(size_t value, const int target_rank) {
        size_t res = 0;
        MPI_Reduce(&value, &res, 1, MY_MPI_SIZE_T, MPI_SUM, target_rank, comm);
        return res;
    }

    double Communicator::sum(double value) {
        double res = 0.0;
        MPI_Allreduce(&value, &res, 1, MPI_DOUBLE, MPI_SUM, comm);
        return res;
    }

    std::vector<float> Communicator::sum(std::vector<float>& values) {
        const int size = static_cast<int>(values.size());

        std::vector<float> res(size);
        MPI_Allreduce(values.data(), res.data(), size, MPI_FLOAT, MPI_SUM, comm);

        return res;
    }

    boost::multi_array<double, 2> Communicator::sum(boost::multi_array<double, 2>& values) {
        boost::multi_array<double, 2> results(boost::extents[ values.shape()[0] ][ values.shape()[1] ]);
        const int size = static_cast<int>(values.num_elements());

        MPI_Allreduce(values.data(), results.data(), size, MPI_DOUBLE, MPI_SUM, comm);
        return results;
    }

    int Communicator::sum(int value) {
        int res = 0;
        MPI_Allreduce(&value, &res, 1, MPI_INT, MPI_SUM, comm);
        return res;
    }

    size_t Communicator::sum(size_t value) {
        size_t res = 0;
        MPI_Allreduce(&value, &res, 1, MY_MPI_SIZE_T, MPI_SUM, comm);
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

    int Communicator::min(int value) {
        int res = 0;
        MPI_Allreduce(&value, &res, 1, MPI_INT, MPI_MIN, comm);
        return res;
    }

    std::string Communicator::gatherStringsTo(int target_rank, std::string& content) const {
        std::string result = "";
        if (Environment::rank == target_rank) {
            for (int source_rank = 0; source_rank < Environment::numprocs; ++source_rank) {
                if (source_rank == target_rank) {
                    result += content;
                    continue;
                }

                MPI_Status status;
                MPI_Probe(source_rank, TAG::GATHER_STRINGS, comm, &status);

                int num_of_received_elements;
                MPI_Get_count(&status, MPI_CHAR, &num_of_received_elements);

                char* buff = new char[ num_of_received_elements + 1 ];
                MPI_Recv(buff, num_of_received_elements, MPI_CHAR, source_rank, TAG::GATHER_STRINGS, comm, &status);
                buff[num_of_received_elements] = '\0';
                std::string temporary_result(buff);
                delete [] buff;

                result += temporary_result;
            }
        } else {
            MPI_Send(content.c_str(), static_cast<int>(content.size()), MPI_CHAR, target_rank, TAG::GATHER_STRINGS, comm);
        }
        return result;
    }

    void Communicator::sendRecvPartialFieldX(tdArray& tdValue, const size_t i_index, const size_t j_begin, const size_t j_end, const size_t k_begin, const size_t k_end, const int send_target, const int recv_target) {
        MPI_Status s;
        const size_t j_size = j_end - j_begin + 1;
        const size_t k_size = k_end - k_begin + 1;
        const size_t length = j_size * k_size;

        boost::multi_array<double, 2> buff(boost::extents[j_size][k_size]);

        for(int j = 0; j < j_size; ++j) {
            for(int k = 0; k < k_size; ++k) {
                buff[j][k] = tdValue[i_index][j + j_begin][k + k_begin];
            }
        }

        MPI_Sendrecv_replace(buff.data(), length, MPI_DOUBLE, send_target, TAG::SENDRECV_FIELD, recv_target, TAG::SENDRECV_FIELD, comm, &s);

        const size_t i_index_inv = (i_index == 1) ? tdValue.shape()[0] - 1 : 0;
        for(int j = 0; j < j_size; ++j) {
            for(int k = 0; k < k_size; ++k) {
                tdValue[i_index_inv][j + j_begin][k + k_begin] = buff[j][k];
            }
        }
    }

    void Communicator::sendRecvPartialFieldY(tdArray& tdValue, const size_t i_begin, const size_t i_end, const size_t j_index, const size_t k_begin, const size_t k_end, const int send_target, const int recv_target) {
        MPI_Status s;
        const size_t i_size = i_end - i_begin + 1;
        const size_t k_size = k_end - k_begin + 1;
        const size_t length = i_size * k_size;

        boost::multi_array<double, 2> buff(boost::extents[i_size][k_size]);

        for(int i = 0; i < i_size; ++i) {
            for(int k = 0; k < k_size; ++k) {
                buff[i][k] = tdValue[i + i_begin][j_index][k + k_begin];
            }
        }

        MPI_Sendrecv_replace(buff.data(), length, MPI_DOUBLE, send_target, TAG::SENDRECV_FIELD, recv_target, TAG::SENDRECV_FIELD, comm, &s);

        const size_t j_index_inv = (j_index == 1) ? tdValue.shape()[1] - 1 : 0;
        for(int i = 0; i < i_size; ++i) {
            for(int k = 0; k < k_size; ++k) {
                tdValue[i + i_begin][j_index_inv][k + k_begin] = buff[i][k];
            }
        }
    }

    void Communicator::sendRecvPartialFieldZ(tdArray& tdValue, const size_t i_begin, const size_t i_end, const size_t j_begin, const size_t j_end, const size_t k_index, const int send_target, const int recv_target) {
        MPI_Status s;
        const size_t i_size = i_end - i_begin + 1;
        const size_t j_size = j_end - j_begin + 1;
        const size_t length = i_size * j_size;

        boost::multi_array<double, 2> buff(boost::extents[i_size][j_size]);

        for(int i = 0; i < i_size; ++i) {
            for(int j = 0; j < j_size; ++j) {
                buff[i][j] = tdValue[i + i_begin][j + j_begin][k_index];
            }
        }

        MPI_Sendrecv_replace(buff.data(), length, MPI_DOUBLE, send_target, TAG::SENDRECV_FIELD, recv_target, TAG::SENDRECV_FIELD, comm, &s);

        const size_t k_index_inv = (k_index == 1) ? tdValue.shape()[2] - 1 : 0;
        for(int i = 0; i < i_size; ++i) {
            for(int j = 0; j < j_size; ++j) {
                tdValue[i + i_begin][j + j_begin][k_index_inv] = buff[i][j];
            }
        }
    }

    void Communicator::sendRecvFieldX(tdArray& tdValue, const int prev, const int next) {
        MPI_Status s;
        boost::multi_array<double, 2> buff(boost::extents[ tdValue.shape()[1] ][ tdValue.shape()[2] ]);

        #pragma omp parallel shared(s, tdValue, buff)
        {
            //! 下側の値をバッファへ詰める
            #pragma omp for
            for(int j = 0; j < tdValue.shape()[1]; ++j) {
                for(int k = 0; k < tdValue.shape()[2]; ++k) {
                    buff[j][k] = tdValue[1][j][k];
                }
            }

            if ( (prev == Environment::rank) && (next == Environment::rank) ) {
                //! 自分自身への通信
                //! std::swapのが速いかも
                #pragma omp barrier
                #pragma omp for
                for(int j = 0; j < tdValue.shape()[1]; ++j) {
                    for(int k = 0; k < tdValue.shape()[2]; ++k) {
                        tdValue[0][j][k] = tdValue[tdValue.shape()[0] - 2][j][k];
                        tdValue[tdValue.shape()[0] - 1][j][k] = buff[j][k];
                    }
                }
            } else {
                int length = static_cast<int>(tdValue.shape()[1] * tdValue.shape()[2]);

                #pragma omp barrier
                #pragma omp single
                {
                    MPI_Sendrecv_replace(buff.data(), length, MPI_DOUBLE, prev, TAG::SENDRECV_FIELD, next, TAG::SENDRECV_FIELD, comm, &s);
                }

                #pragma omp for
                for(int j = 0; j < tdValue.shape()[1]; ++j) {
                    for(int k = 0; k < tdValue.shape()[2]; ++k) {
                        tdValue[tdValue.shape()[0] - 1][j][k] = buff[j][k];

                        // 反対側の値を buff に詰める
                        buff[j][k] = tdValue[tdValue.shape()[0] - 2][j][k];
                    }
                }

                #pragma omp barrier
                #pragma omp single
                {
                    MPI_Sendrecv_replace(buff.data(), length, MPI_DOUBLE, next, TAG::SENDRECV_FIELD, prev, TAG::SENDRECV_FIELD, comm, &s);
                }

                #pragma omp for
                for(int j = 0; j < tdValue.shape()[1]; ++j) {
                    for(int k = 0; k < tdValue.shape()[2]; ++k) {
                        tdValue[0][j][k] = buff[j][k];
                    }
                }
            }
        }
    }

    void Communicator::sendRecvFieldY(tdArray& tdValue, const int prev, const int next) {
        MPI_Status s;
        boost::multi_array<double, 2> buff(boost::extents[ tdValue.shape()[0] ][ tdValue.shape()[2] ]);

        #pragma omp parallel shared(s, tdValue, buff)
        {
            #pragma omp for
            for(int i = 0; i < tdValue.shape()[0]; ++i) {
                for(int k = 0; k < tdValue.shape()[2]; ++k) {
                    buff[i][k] = tdValue[i][1][k];
                }
            }

            if ( (prev == Environment::rank) && (next == Environment::rank) ) {
                #pragma omp barrier
                #pragma omp for
                for(int i = 0; i < tdValue.shape()[0]; ++i) {
                    for(int k = 0; k < tdValue.shape()[2]; ++k) {
                        tdValue[i][0][k] = tdValue[i][tdValue.shape()[1] - 2][k];
                        tdValue[i][tdValue.shape()[1] - 1][k] = buff[i][k];
                    }
                }
            } else {
                int length = static_cast<int>(tdValue.shape()[0] * tdValue.shape()[2]);

                #pragma omp barrier
                #pragma omp single
                {
                    MPI_Sendrecv_replace(buff.data(), length, MPI_DOUBLE, prev, TAG::SENDRECV_FIELD, next, TAG::SENDRECV_FIELD, comm, &s);
                }

                #pragma omp for
                for(int i = 0; i < tdValue.shape()[0]; ++i) {
                    for(int k = 0; k < tdValue.shape()[2]; ++k) {
                        tdValue[i][tdValue.shape()[1] - 1][k] = buff[i][k];
                        buff[i][k] = tdValue[i][tdValue.shape()[1] - 2][k];
                    }
                }

                #pragma omp barrier
                #pragma omp single
                {
                    MPI_Sendrecv_replace(buff.data(), length, MPI_DOUBLE, next, TAG::SENDRECV_FIELD, prev, TAG::SENDRECV_FIELD, comm, &s);
                }

                #pragma omp for
                for(int i = 0; i < tdValue.shape()[0]; ++i) {
                    for(int k = 0; k < tdValue.shape()[2]; ++k) {
                        tdValue[i][0][k] = buff[i][k];
                    }
                }
            }
        }
    }

    void Communicator::sendRecvFieldZ(tdArray& tdValue, const int prev, const int next) {
        MPI_Status s;
        boost::multi_array<double, 2> buff(boost::extents[ tdValue.shape()[0] ][ tdValue.shape()[1] ]);

        #pragma omp parallel shared(s, tdValue, buff)
        {
            #pragma omp for
            for(int i = 0; i < tdValue.shape()[0]; ++i) {
                for(int j = 0; j < tdValue.shape()[1]; ++j) {
                    buff[i][j] = tdValue[i][j][1];
                }
            }

            if ( (prev == Environment::rank) && (next == Environment::rank) ) {
                #pragma omp barrier
                #pragma omp for
                for(int i = 0; i < tdValue.shape()[0]; ++i) {
                    for(int j = 0; j < tdValue.shape()[1]; ++j) {
                        tdValue[i][j][0] = tdValue[i][j][tdValue.shape()[2] - 2];
                        tdValue[i][j][tdValue.shape()[2] - 1] = buff[i][j];
                    }
                }
            } else {
                int length = static_cast<int>(tdValue.shape()[0] * tdValue.shape()[1]);

                #pragma omp barrier
                #pragma omp single
                {
                    MPI_Sendrecv_replace(buff.data(), length, MPI_DOUBLE, prev, TAG::SENDRECV_FIELD, next, TAG::SENDRECV_FIELD, comm, &s);
                }

                #pragma omp for
                for(int i = 0; i < tdValue.shape()[0]; ++i) {
                    for(int j = 0; j < tdValue.shape()[1]; ++j) {
                        tdValue[i][j][tdValue.shape()[2] - 1] = buff[i][j];
                        buff[i][j] = tdValue[i][j][tdValue.shape()[2] - 2];
                    }
                }

                #pragma omp barrier
                #pragma omp single
                {
                    MPI_Sendrecv_replace(buff.data(), length, MPI_DOUBLE, next, TAG::SENDRECV_FIELD, prev, TAG::SENDRECV_FIELD, comm, &s);
                }

                #pragma omp for
                for(int i = 0; i < tdValue.shape()[0]; ++i) {
                    for(int j = 0; j < tdValue.shape()[1]; ++j) {
                        tdValue[i][j][0] = buff[i][j];
                    }
                }
            }
        }
    }

    void Communicator::sendRecvScalarX(tdArray& tdValue, const int prev, const int next, const int glue_size) {
        MPI_Status s;

        const auto size_x = tdValue.shape()[0];
        const auto size_y = tdValue.shape()[1];
        const auto size_z = tdValue.shape()[2];

        boost::multi_array<double, 3> buff_upper(boost::extents[glue_size][size_y][size_z]);
        boost::multi_array<double, 3> buff_lower(boost::extents[glue_size][size_y][size_z]);

        #pragma omp parallel
        {
            #pragma omp for
            for(int j = 0; j < size_y; ++j) {
                for(int k = 0; k < size_z; ++k) {
                    for (int glue_index = 0; glue_index < glue_size; ++glue_index) {
                        //! 下側の値をバッファへ詰める
                        buff_lower[glue_index][j][k] = tdValue[(glue_size - 1) - glue_index][j][k];

                        //! 上側の値をバッファへ詰める
                        buff_upper[glue_index][j][k] = tdValue[size_x - 1 - glue_index][j][k];
                    }
                }
            }

            #pragma omp barrier

            if ( (prev == Environment::rank) && (next == Environment::rank) ) {
                //! 自分自身への通信
                //! std::swapのが速いかも
                #pragma omp for
                for(int j = 0; j < size_y; ++j) {
                    for(int k = 0; k < size_z; ++k) {
                        for (int glue_index = 0; glue_index < glue_size; ++glue_index) {
                            tdValue[(glue_size - 1) - glue_index][j][k] += buff_upper[glue_index][j][k];
                            tdValue[size_x - 1 - glue_index][j][k] += buff_lower[glue_index][j][k];
                        }
                    }
                }
            } else {
                int length = static_cast<int>(glue_size * size_y * size_z);

                #pragma omp single
                {
                    MPI_Sendrecv_replace(buff_upper.data(), length, MPI_DOUBLE, next, TAG::SENDRECV_SCALAR_UPPER, prev, TAG::SENDRECV_SCALAR_UPPER, comm, &s);
                }
                // buff_upper に上側からの値が入る (0->1の順で外側から)

                #pragma omp single
                {
                    MPI_Sendrecv_replace(buff_lower.data(), length, MPI_DOUBLE, prev, TAG::SENDRECV_SCALAR_LOWER, next, TAG::SENDRECV_SCALAR_LOWER, comm, &s);
                }
                // buff_lower に上側からの値が入る (0->1の順で外側から)

                if (::Environment::isNotBoundary(AXIS::x, AXIS_SIDE::up)) {
                    #pragma omp for
                    for(int j = 0; j < size_y; ++j) {
                        for(int k = 0; k < size_z; ++k) {
                            for (int glue_index = 0; glue_index < glue_size; ++glue_index) {
                                tdValue[size_x - 1 - glue_index][j][k] += buff_lower[glue_index][j][k];
                            }
                        }
                    }
                }

                if (::Environment::isNotBoundary(AXIS::x, AXIS_SIDE::low)) {
                    #pragma omp for
                    for(int j = 0; j < size_y; ++j) {
                        for(int k = 0; k < size_z; ++k) {
                            for (int glue_index = 0; glue_index < glue_size; ++glue_index) {
                                tdValue[(glue_size - 1) - glue_index][j][k] += buff_upper[glue_index][j][k];
                            }
                        }
                    }
                }
            }
        }
    }

    void Communicator::sendRecvScalarY(tdArray& tdValue, const int prev, const int next, const int glue_size) {
        MPI_Status s;

        const auto size_x = tdValue.shape()[0];
        const auto size_y = tdValue.shape()[1];
        const auto size_z = tdValue.shape()[2];

        boost::multi_array<double, 3> buff_upper(boost::extents[glue_size][size_x][size_z]);
        boost::multi_array<double, 3> buff_lower(boost::extents[glue_size][size_x][size_z]);

        #pragma omp parallel
        {
            #pragma omp for
            for(int i = 0; i < size_x; ++i) {
                for(int k = 0; k < size_z; ++k) {
                    for (int glue_index = 0; glue_index < glue_size; ++glue_index) {
                        buff_lower[glue_index][i][k] = tdValue[i][(glue_size - 1) - glue_index][k];
                        buff_upper[glue_index][i][k] = tdValue[i][size_y - 1 - glue_index][k];
                    }
                }
            }

            if ( (prev == Environment::rank) && (next == Environment::rank) ) {
                //! 自分自身への通信
                //! std::swapのが速いかも
                #pragma omp barrier
                #pragma omp for
                for(int i = 0; i < size_x; ++i) {
                    for(int k = 0; k < size_z; ++k) {
                        for (int glue_index = 0; glue_index < glue_size; ++glue_index) {
                            tdValue[i][(glue_size - 1) - glue_index][k] += buff_upper[glue_index][i][k];
                            tdValue[i][size_y - 1 - glue_index][k] += buff_lower[glue_index][i][k];
                        }
                    }
                }
            } else {
                int length = static_cast<int>(glue_size * size_x * size_z);

                #pragma omp barrier
                #pragma omp single
                {
                    MPI_Sendrecv_replace(buff_upper.data(), length, MPI_DOUBLE, next, TAG::SENDRECV_SCALAR_UPPER, prev, TAG::SENDRECV_SCALAR_UPPER, comm, &s);
                }
                // buff_upper に下側からの値が入る

                #pragma omp single
                {
                    MPI_Sendrecv_replace(buff_lower.data(), length, MPI_DOUBLE, prev, TAG::SENDRECV_SCALAR_LOWER, next, TAG::SENDRECV_SCALAR_LOWER, comm, &s);
                }
                // buff_lower に上側からの値が入る

                if (::Environment::isNotBoundary(AXIS::y, AXIS_SIDE::up)) {
                    #pragma omp for
                    for(int i = 0; i < size_x; ++i) {
                        for(int k = 0; k < size_z; ++k) {
                            for (int glue_index = 0; glue_index < glue_size; ++glue_index) {
                                tdValue[i][size_y - 1 - glue_index][k] += buff_lower[glue_index][i][k];
                            }
                        }
                    }
                }

                if (::Environment::isNotBoundary(AXIS::y, AXIS_SIDE::low)) {
                    #pragma omp for
                    for(int i = 0; i < size_x; ++i) {
                        for(int k = 0; k < size_z; ++k) {
                            for (int glue_index = 0; glue_index < glue_size; ++glue_index) {
                                tdValue[i][(glue_size - 1) - glue_index][k] += buff_upper[glue_index][i][k];
                            }
                        }
                    }
                }
            }
        }
    }

    void Communicator::sendRecvScalarZ(tdArray& tdValue, const int prev, const int next, const int glue_size) {
        MPI_Status s;

        const auto size_x = tdValue.shape()[0];
        const auto size_y = tdValue.shape()[1];
        const auto size_z = tdValue.shape()[2];

        boost::multi_array<double, 3> buff_upper(boost::extents[glue_size][size_x][size_y]);
        boost::multi_array<double, 3> buff_lower(boost::extents[glue_size][size_x][size_y]);

        #pragma omp parallel
        {
            #pragma omp for
            for(int i = 0; i < size_x; ++i) {
                for(int j = 0; j < size_y; ++j) {
                    for (int glue_index = 0; glue_index < glue_size; ++glue_index) {
                        buff_lower[glue_index][i][j] = tdValue[i][j][(glue_size - 1) - glue_index];
                        buff_upper[glue_index][i][j] = tdValue[i][j][size_z - 1 - glue_index];
                    }
                }
            }

            #pragma omp barrier

            if ( (prev == Environment::rank) && (next == Environment::rank) ) {
                #pragma omp for
                for(int i = 0; i < size_x; ++i) {
                    for(int j = 0; j < size_y; ++j) {
                        for (int glue_index = 0; glue_index < glue_size; ++glue_index) {
                            tdValue[i][j][(glue_size - 1) - glue_index] += buff_upper[glue_index][i][j];
                            tdValue[i][j][size_z - 1 - glue_index] += buff_lower[glue_index][i][j];
                        }
                    }
                }
            } else {
                int length = static_cast<int>(glue_size * size_x * size_y);

                #pragma omp single
                {
                    MPI_Sendrecv_replace(buff_upper.data(), length, MPI_DOUBLE, next, TAG::SENDRECV_SCALAR_UPPER, prev, TAG::SENDRECV_SCALAR_UPPER, comm, &s);
                }
                // buff_upper に下側からの値が入る

                #pragma omp single
                {
                    MPI_Sendrecv_replace(buff_lower.data(), length, MPI_DOUBLE, prev, TAG::SENDRECV_SCALAR_LOWER, next, TAG::SENDRECV_SCALAR_LOWER, comm, &s);
                }
                // buff_lower に上側からの値が入る

                if (::Environment::isNotBoundary(AXIS::z, AXIS_SIDE::up)) {
                    #pragma omp for
                    for(int i = 0; i < size_x; ++i) {
                        for(int j = 0; j < size_y; ++j) {
                            for (int glue_index = 0; glue_index < glue_size; ++glue_index) {
                                tdValue[i][j][size_z - 1 - glue_index] += buff_lower[glue_index][i][j];
                            }
                        }
                    }
                }

                if (::Environment::isNotBoundary(AXIS::z, AXIS_SIDE::low)) {
                    #pragma omp for
                    for(int i = 0; i < size_x; ++i) {
                        for(int j = 0; j < size_y; ++j) {
                            for (int glue_index = 0; glue_index < glue_size; ++glue_index) {
                                tdValue[i][j][(glue_size - 1) - glue_index] += buff_upper[glue_index][i][j];
                            }
                        }
                    }
                }
            }
        }
    }
}