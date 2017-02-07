#ifndef __TDPIC_FIELD_H_INCLUDED__
#define __TDPIC_FIELD_H_INCLUDED__
#include "global.hpp"
#include <mkl.h>

struct Poisson {
    public:
        ~Poisson();

        MKL_INT nx;
        MKL_INT ny;
        MKL_INT nz;

        MKL_INT* ipar;
        double* dpar;
        MKL_INT stat;
        DFTI_DESCRIPTOR_HANDLE* xhandle;
        DFTI_DESCRIPTOR_HANDLE* yhandle;
        double* b_lx;
        double* b_ly;
        double* b_lz;

        double* rho1D;
};

class Field {
    private:
        tdArray rho;
        tdArray phi;

        tdArray ex;
        tdArray ey;
        tdArray ez;
        tdArray exref;
        tdArray eyref;
        tdArray ezref;

        tdArray bx;
        tdArray by;
        tdArray bz;

        Poisson* psn = nullptr;
    public:
        Field();
        ~Field();

        // scalars
        tdArray& getPhi();
        void setPhi(tdArray&);
        tdArray& getRho();
        void setRho(tdArray&);
        tdArray& getScalar(std::string);

        // vectors
        tdArray& getEx();
        void setEx(tdArray&);
        tdArray& getEy();
        void setEy(tdArray&);
        tdArray& getEz();
        void setEz(tdArray&);

        // reference efield on nodes
        tdArray& getExRef();
        tdArray& getEyRef();
        tdArray& getEzRef();
        void setExRef(tdArray&);
        void setEyRef(tdArray&);
        void setEzRef(tdArray&);

        tdArray& getBx();
        void setBx(tdArray&);
        tdArray& getBy();
        void setBy(tdArray&);
        tdArray& getBz();
        void setBz(tdArray&);

        void initializePoisson(const int, const int, const int);
        void updateEfield(const int, const int, const int);
        void updateBfield(const int, const int, const int);

        //! solver の実体
        void solvePoissonMKL(const int, const int, const int);
        void solvePoissonSOR(const int, const double);
        void solvePoisson(const int, const double);
};


#endif
