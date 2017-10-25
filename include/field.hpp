#ifndef __TDPIC_FIELD_H_INCLUDED__
#define __TDPIC_FIELD_H_INCLUDED__
#include "global.hpp"

using RhoArray = std::vector<tdArray>;

class Field {
    private:
        //! rhoとdensityのみ粒子種で必要な行列数が変わる
        RhoArray rho;
        RhoArray density;

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
        tdArray jxref;
        tdArray jyref;
        tdArray jzref;

        void setDirichletPhi(const AXIS, const AXIS_SIDE);
        void setNeumannPhi(const AXIS, const AXIS_SIDE);

    public:
        Field(void) :
            rho{}, density{},
            poisson_residual{}, poisson_error{}, phi{},
            ex{}, ey{}, ez{},
            exref{}, eyref{}, ezref{},
            bx{}, by{}, bz{},
            bxref{}, byref{}, bzref{},
            jx{}, jy{}, jz{},
            jxref{}, jyref{}, jzref{}
            {}

        // destructor
        ~Field(){}

        // potential
        tdArray& getPhi(){ return phi; }
        tdArray& getPoissonResidual(){ return poisson_residual; }
        tdArray& getPoissonError(){ return poisson_error; }

        double poissonOperator(const tdArray& phi, const size_t i, const size_t j, const size_t k) const {
            return (phi[i-1][j][k] + phi[i+1][j][k] + phi[i][j-1][k] + phi[i][j+1][k] + phi[i][j][k-1] + phi[i][j][k+1] - 6.0*phi[i][j][k]);
        }

        // charge density
        RhoArray& getRho(){ return rho; }

        // density
        RhoArray& getDensity(){ return density; }

        // electric fields
        tdArray& getEx(){ return ex; }
        tdArray& getEy(){ return ey; }
        tdArray& getEz(){ return ez; }

        //! reference efield on nodes
        tdArray& getExRef(){ return exref; }
        tdArray& getEyRef(){ return eyref; }
        tdArray& getEzRef(){ return ezref; }

        // magnetic fields
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

        //! reference current density (for data output)
        tdArray& getJxRef(){ return jxref; }
        tdArray& getJyRef(){ return jyref; }
        tdArray& getJzRef(){ return jzref; }

        //! boundary conditions
        void setBoundaryConditionPhi(void);
        void setAbsorbingBoundaryOnEfield(void);
        void setDampingBoundaryOnEfield(void);

        void initializeCurrent(const double dt);

        double getEfieldEnergy(void) const;
        double getBfieldEnergy(void) const;

        //! check functions
        void checkChargeConservation(const RhoArray&, const double, const double) const;
};

#endif
