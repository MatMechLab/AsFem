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
//+++ Date   : 2021.04.10
//+++ Purpose: implement the residual and jacobian for Poisson equation
//+++          Sigma*div(grad(phi))=F
//+++          Both Sigma and F could be nonlinear function !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/PoissonElmt.h"

void PoissonElmt::ComputeAll(const FECalcType &calctype,const LocalElmtInfo &elmtinfo,const double (&ctan)[2],
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
        MessagePrinter::PrintErrorTxt("unsupported calculation type in PoissonElmt, please check your related code");
        MessagePrinter::AsFem_Exit();
    }
}
//***************************************************************************
void PoissonElmt::ComputeResidual(const LocalElmtInfo &elmtinfo,
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const Materials &Mate,const Materials &MateOld,
                                 VectorXd &localR) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.dt||soln.gpU[0]||shp.test||Mate.GetScalarMate().size()||MateOld.GetScalarMate().size()) {}

    localR(1)=Mate.ScalarMaterials("sigma")*(soln.gpGradU[1]*shp.grad_test)+Mate.ScalarMaterials("f")*shp.test;

}
//*****************************************************************************
void PoissonElmt::ComputeJacobian(const LocalElmtInfo &elmtinfo,const double (&ctan)[2],
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const Materials &Mate,const Materials &MateOld,
                                 MatrixXd &localK) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.dt||ctan[0]||soln.gpU[0]||Mate.GetScalarMate().size()||MateOld.GetScalarMate().size()){}

    localK(1,1)=Mate.ScalarMaterials("dsigmadu")*shp.trial*(soln.gpGradU[1]*shp.grad_test)*ctan[0]
            +Mate.ScalarMaterials("sigma")*shp.grad_trial*shp.grad_test*ctan[0]
            +Mate.ScalarMaterials("dfdu")*shp.trial*shp.test*ctan[0];

}
//*******************************************************************************
void PoissonElmt::ComputeProjection(const LocalElmtInfo &elmtinfo,const double (&ctan)[2],
                                   const LocalElmtSolution &soln,
                                   const LocalShapeFun &shp,
                                   const Materials &Mate,const Materials &MateOld,
                                   ScalarMateType &gpProj) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.dt||ctan[0]||soln.gpU[0]||shp.test||Mate.GetScalarMate().size()||MateOld.GetScalarMate().size()||gpProj.size()){}
}
