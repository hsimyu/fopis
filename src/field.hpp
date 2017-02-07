#ifndef __TDPIC_FIELD_H_INCLUDED__
#define __TDPIC_FIELD_H_INCLUDED__
#include "global.hpp"
#include <mkl.h>

struct Poisson {
    public:
        ~Poisson(void){
            delete [] ipar;
            delete [] dpar;
            delete [] b_lx;
            delete [] b_ly;
            delete [] b_lz;
            delete [] xhandle;
            delete [] yhandle;
            delete [] rho1D;
        }

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
        Field(void) :
            rho(boost::extents[0][0][0], boost::fortran_storage_order()),
            phi(boost::extents[0][0][0], boost::fortran_storage_order()),
            ex(boost::extents[0][0][0], boost::fortran_storage_order()),
            ey(boost::extents[0][0][0], boost::fortran_storage_order()),
            ez(boost::extents[0][0][0], boost::fortran_storage_order()),
            exref(boost::extents[0][0][0], boost::fortran_storage_order()),
            eyref(boost::extents[0][0][0], boost::fortran_storage_order()),
            ezref(boost::extents[0][0][0], boost::fortran_storage_order()),
            bx(boost::extents[0][0][0], boost::fortran_storage_order()),
            by(boost::extents[0][0][0], boost::fortran_storage_order()),
            bz(boost::extents[0][0][0], boost::fortran_storage_order()) {}

        // destructor
        ~Field(){
            delete psn;
        }

        // potential
        void setPhi(tdArray& _phi){
            phi = _phi;
        }
        tdArray& getPhi(){
            return phi;
        }

        // charge density
        void setRho(tdArray& _rho){
            rho = _rho;
        }
        tdArray& getRho(){
            return rho;
        }

        tdArray& getScalar(std::string varname){
            if(varname == "potential" || varname == "phi") {
                return phi;
            } else if(varname == "rho" || varname == "charge") {
                return rho;
            } else {
                throw std::invalid_argument("[ERROR] Invalid varname is passed to getScalarField():" + varname);
            }
        }

        // electric fields
        void setEx(tdArray& _ex){
            ex = _ex;
        }
        tdArray& getEx(){
            return ex;
        }

        void setEy(tdArray& _ey){
            ey = _ey;
        }
        tdArray& getEy(){
            return ey;
        }

        void setEz(tdArray& _ez){
            ez = _ez;
        }
        tdArray& getEz(){
            return ez;
        }

        //! reference efield on nodes
        tdArray& getExRef(){
            return exref;
        }

        tdArray& getEyRef(){
            return eyref;
        }

        tdArray& getEzRef(){
            return ezref;
        }

        // magnetic fields
        void setBx(tdArray& _bx){
            bx = _bx;
        }
        tdArray& getBx(){
            return bx;
        }

        void setBy(tdArray& _by){
            by = _by;
        }
        tdArray& getBy(){
            return by;
        }

        void setBz(tdArray& _bz){
            bz = _bz;
        }
        tdArray& getBz(){
            return bz;
        }

        void initializePoisson(const int, const int, const int);
        void updateEfield(const double, const int, const int, const int);
        void updateBfield(const double, const int, const int, const int);

        //! solver の実体
        void solvePoissonMKL(const int, const int, const int);
        void solvePoissonSOR(const int, const double);
        void solvePoisson(const int, const double);

        double getEnergy(const int, const int, const int);
};


#endif
