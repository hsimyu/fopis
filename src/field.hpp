#ifndef __TDPIC_FIELD_H_INCLUDED__
#define __TDPIC_FIELD_H_INCLUDED__
#include "global.hpp"

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

        void setBoundaryConditionPhi(void);
        void setDirichletPhi(const std::string, const std::string);

        //! solver の実体
        void solvePoissonPSOR(const int, const double);
        void solvePoissonJacobi(const int, const double);
        double checkPhiResidual(void);
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
        ~Field(){}

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

        void updateEfield(const double);
        void updateBfield(const double, const int, const int, const int);
        void solvePoisson(const int, const double);

        double getEnergy(const int, const int, const int);
};


#endif
