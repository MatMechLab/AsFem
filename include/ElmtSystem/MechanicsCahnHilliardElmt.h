//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.11.21
//+++ Purpose: implement the residual and jacobian for general
//+++          stress-chan hilliard equilibrium equation
//+++          governing equation:
//+++          1) dc/dt=div(M*grad(mu))
//+++          2) mu   = dF/dc-kappa*grad(c)
//+++          3) div(sigma)=0
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/BulkElmtBase.h"


/**
 * This class implement the calculation of mechanically coupled CahnHilliard equation
 */
class MechanicsCahnHilliardElmt:public BulkElmtBase{
public:
    /**
     * Function for different calculate action
     */
    virtual void ComputeAll(const FECalcType &calctype,const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
            const LocalElmtSolution &soln,const LocalShapeFun &shp,
            const Materials &Mate,const Materials &MateOld,
            ScalarMateType &gpProj,
            MatrixXd &localK,VectorXd &localR) override;

private:
    /**
     * This function calculate the residual of the stress equilibrium and CahnHilliard equation. <br>
     */
    virtual void ComputeResidual(const LocalElmtInfo &elmtinfo,
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const Materials &Mate,const Materials &MateOld,
                                 VectorXd &localR) override;

    /**
     * This function calculate the jacobian matrix of the stress equilibrium and CahnHilliard quation. <br>
     */
    virtual void ComputeJacobian(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const Materials &Mate,const Materials &MateOld,
                                 MatrixXd &localK) override;
    
    /**
     */
    virtual void ComputeProjection(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                   const LocalElmtSolution &soln,
                                   const LocalShapeFun &shp,
                                   const Materials &Mate,const Materials &MateOld,
                                   ScalarMateType &gpProj) override;

private:
    RankTwoTensor Stress,dStressdC,dMudStrain;
    double M,dMdC;
    double dFdC,d2FdC2,kappa;
    double valx,valy,valz;
};
