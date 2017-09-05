#ifndef __TDPIC_FIELD_H_INCLUDED__
#define __TDPIC_FIELD_H_INCLUDED__
#include "global.hpp"

using RhoArray = std::vector<tdArray>;

class Field {
    private:
        //! rho のみ粒子種で必要な行列数が変わる
        RhoArray rho;

        //! マルチグリッド Poisson 計算用
        tdArray poisson_residual;
        tdArray poisson_error;

        //! 場
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

        void setDirichletPhi(const AXIS, const AXIS_SIDE);
        void setNeumannPhi(const AXIS, const AXIS_SIDE);

    public:
        Field(void) :
            rho{},
            poisson_residual(boost::extents[0][0][0], boost::fortran_storage_order()),
            poisson_error(boost::extents[0][0][0], boost::fortran_storage_order()),
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
        tdArray& getPhi(){ return phi; }
        tdArray& getPoissonResidual(){ return poisson_residual; }
        tdArray& getPoissonError(){ return poisson_error; }
        void setBoundaryConditionPhi(void);

        double poissonOperator(const tdArray& phi, const int i, const int j, const int k) const {
            return (phi[i-1][j][k] + phi[i+1][j][k] + phi[i][j-1][k] + phi[i][j+1][k] + phi[i][j][k-1] + phi[i][j][k+1] - 6.0*phi[i][j][k]);
        }

        // charge density
        void setRho(RhoArray& _rho){ rho = _rho; }
        RhoArray& getRho(){ return rho; }

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
        void initializeCurrent(const double dt);

        double getEfieldEnergy(void) const;
        double getBfieldEnergy(void) const;

        //! check functions
        void checkChargeConservation(const RhoArray&, const double, const double) const;
};

#endif
