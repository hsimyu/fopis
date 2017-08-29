#ifndef __TDPIC_FIELD_H_INCLUDED__
#define __TDPIC_FIELD_H_INCLUDED__
#include "global.hpp"

using RhoArray = std::vector<tdArray>;

class Field {
    private:
        //! rho のみ粒子種で必要な行列数が変わる
        RhoArray rho;

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
        tdArray bxref;
        tdArray byref;
        tdArray bzref;
        tdArray jx;
        tdArray jy;
        tdArray jz;

        void setBoundaryConditionPhi(void);
        void setDirichletPhi(const AXIS, const AXIS_SIDE);
        void setNeumannPhi(const AXIS, const AXIS_SIDE);

        //! solver の実体
        void solvePoissonPSOR(const int, const double);
        void solvePoissonJacobi(const int, const double);
        double checkPhiResidual(void);

    public:
        Field(void) :
            rho{},
            phi(boost::extents[0][0][0], boost::fortran_storage_order()),
            ex(boost::extents[0][0][0], boost::fortran_storage_order()),
            ey(boost::extents[0][0][0], boost::fortran_storage_order()),
            ez(boost::extents[0][0][0], boost::fortran_storage_order()),
            exref(boost::extents[0][0][0], boost::fortran_storage_order()),
            eyref(boost::extents[0][0][0], boost::fortran_storage_order()),
            ezref(boost::extents[0][0][0], boost::fortran_storage_order()),
            bx(boost::extents[0][0][0], boost::fortran_storage_order()),
            by(boost::extents[0][0][0], boost::fortran_storage_order()),
            bz(boost::extents[0][0][0], boost::fortran_storage_order()),
            bxref(boost::extents[0][0][0], boost::fortran_storage_order()),
            byref(boost::extents[0][0][0], boost::fortran_storage_order()),
            bzref(boost::extents[0][0][0], boost::fortran_storage_order()),
            jx(boost::extents[0][0][0], boost::fortran_storage_order()),
            jy(boost::extents[0][0][0], boost::fortran_storage_order()),
            jz(boost::extents[0][0][0], boost::fortran_storage_order()) {}

        // destructor
        ~Field(){}

        // potential
        void setPhi(tdArray& _phi){ phi = _phi; }
        tdArray& getPhi(){ return phi; }

        // charge density
        void setRho(RhoArray& _rho){ rho = _rho; }
        RhoArray& getRho(){ return rho; }

        auto& getScalar(const std::string varname) {
            if(varname == "potential" || varname == "phi") {
                return phi;
            } else if(varname == "rho" || varname == "charge") {
                return phi;
            } else {
                throw std::invalid_argument("[ERROR] Invalid varname is passed to getScalar():" + varname);
            }
        }

        // electric fields
        void setEx(tdArray& _ex){ ex = _ex; }
        void setEy(tdArray& _ey){ ey = _ey; }
        void setEz(tdArray& _ez){ ez = _ez; }
        tdArray& getEx(){ return ex; }
        tdArray& getEy(){ return ey; }
        tdArray& getEz(){ return ez; }

        //! reference efield on nodes
        tdArray& getExRef(){ return exref; }
        tdArray& getEyRef(){ return eyref; }
        tdArray& getEzRef(){ return ezref; }

        // magnetic fields
        void setBx(tdArray& _bx){ bx = _bx; }
        void setBy(tdArray& _by){ by = _by; }
        void setBz(tdArray& _bz){ bz = _bz; }
        tdArray& getBx(){ return bx; }
        tdArray& getBy(){ return by; }
        tdArray& getBz(){ return bz; }

        //! reference bfield on nodes
        tdArray& getBxRef(){ return bxref; }
        tdArray& getByRef(){ return byref; }
        tdArray& getBzRef(){ return bzref; }

        //! current density
        tdArray& getJx(){ return jx; }
        tdArray& getJy(){ return jy; }
        tdArray& getJz(){ return jz; }

        void updateEfield(const double dx);
        void updateEfieldFDTD(const double dx, const double dt);
        void updateReferenceEfield(void);
        void setDampingBoundaryOnEfield(void);

        void updateBfield(const double dx, const int, const int, const int, const double dt);
        void solvePoisson(const int, const double);

        void initializeCurrent(const double dt);

        double getEfieldEnergy(void) const;
        double getBfieldEnergy(void) const;

        //! check functions
        void checkChargeConservation(const RhoArray&, const double, const double) const;
};

#endif