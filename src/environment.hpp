#ifndef __TDPIC_ENVIRONMENT_H_INCLUDED__
#define __TDPIC_ENVIRONMENT_H_INCLUDED__
class ParticleType;

struct Environment {
    public:
        static int max_particle_num;
        static int num_of_particle_types;
        static double dx;
        static double dt;
        static int nx, ny, nz;
        static int proc_x, proc_y, proc_z;
        static int cell_x, cell_y, cell_z;
        static int max_iteration;
        static bool isRootNode;

        static bool onLowXedge, onHighXedge;
        static bool onLowYedge, onHighYedge;
        static bool onLowZedge, onHighZedge;

        static std::string jobtype;
        static std::string solver_type;
        static std::string boundary;
        static std::string dimension;

        static ParticleType* ptype;

        //! 中身はMPI::Environment::getRankStr()と同様
        static std::string rankStr(void);
        static void printInfo(void);
};
#endif
