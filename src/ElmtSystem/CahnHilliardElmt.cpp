//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.04.09
//+++ Purpose: implement the residual and jacobian for general
//+++          Cahn-Hilliard equation
//+++          dc/dt=div(M*grad(mu))
//+++             mu=df/dc-Kappa*Delta(c)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/CahnHilliardElmt.h"

void CahnHilliardElmt::ComputeAll(const FECalcType &calctype,const LocalElmtInfo &elmtinfo,
                                  const double (&ctan)[3],
                                  const LocalElmtSolution &soln,const LocalShapeFun &shp,
                                  const Materials &Mate,const Materials &MateOld,
                                  ScalarMateType &gpProj,
                                  MatrixXd &localK,VectorXd &localR) {
    if(calctype==FECalcType::ComputeResidual){
        ComputeResidual(elmtinfo,soln,shp,Mate,MateOld,localR);
    }
    else if(calctype==FECalcType::ComputeJacobian){
        ComputeJacobian(elmtinfo,ctan,soln,shp,Mate,MateOld,localK);
    }
    else if(calctype==FECalcType::Projection){
        ComputeProjection(elmtinfo,ctan,soln,shp,Mate,MateOld,gpProj);
    }
    else{
        MessagePrinter::PrintErrorTxt("unsupported calculation type in CahnHilliardElmt, please check your related code");
        MessagePrinter::AsFem_Exit();
    }
}
//*******************************************************************
void CahnHilliardElmt::ComputeResidual(const LocalElmtInfo &elmtinfo,
                                    const LocalElmtSolution &soln,
                                    const LocalShapeFun &shp,
                                    const Materials &Mate,const Materials &MateOld,
                                    VectorXd &localR) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.nDim||soln.gpU.size()||shp.test||Mate.GetScalarMate().size()||MateOld.GetScalarMate().size()){}
    // For R_c
    localR(1)=soln.gpV[1]*shp.test+Mate.ScalarMaterials("M")*(soln.gpGradU[2]*shp.grad_test);
    // For R_mu
    localR(2)=soln.gpU[2]*shp.test-Mate.ScalarMaterials("dFdc")*shp.test
            -Mate.ScalarMaterials("Kappa")*(soln.gpGradU[1]*shp.grad_test);
}
//***************************************************************************
void CahnHilliardElmt::ComputeJacobian(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                    const LocalElmtSolution &soln,
                                    const LocalShapeFun &shp,
                                    const Materials &Mate,const Materials &MateOld,
                                    MatrixXd &localK) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.dt||soln.gpU.size()||shp.test||Mate.GetVectorMate().size()||MateOld.GetVectorMate().size()){}
    // K_c,c
    localK(1,1)=shp.trial*shp.test*ctan[1]
                +Mate.ScalarMaterials("dMdc")*shp.trial*(soln.gpGradU[2]*shp.grad_test)*ctan[0];
    // K_c,mu
    localK(1,2)=Mate.ScalarMaterials("M")*shp.grad_trial*shp.grad_test;
    // K_mu,c
    localK(2,1)=-Mate.ScalarMaterials("d2Fdc2")*shp.trial*shp.test*ctan[0]
            -Mate.ScalarMaterials("Kappa")*shp.grad_trial*shp.grad_test*ctan[0];
    // K_mu,mu
    localK(2,2)=shp.trial*shp.test*ctan[0];
}
//***********************************************************
void CahnHilliardElmt::ComputeProjection(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                         const LocalElmtSolution &soln,
                                         const LocalShapeFun &shp,
                                         const Materials &Mate,const Materials &MateOld,
                                         ScalarMateType &gpProj) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.dt||ctan[0]||soln.gpU.size()||shp.test||Mate.GetScalarMate().size()||MateOld.GetScalarMate().size()||gpProj.size()){}
}
