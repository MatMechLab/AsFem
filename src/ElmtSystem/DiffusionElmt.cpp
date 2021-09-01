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
//+++ Purpose: implement the residual and jacobian for general
//+++          diffusion equation:
//+++          dc/dt=div(D*grad(c))
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/DiffusionElmt.h"

void DiffusionElmt::ComputeAll(const FECalcType &calctype,const LocalElmtInfo &elmtinfo,
                               const double (&ctan)[2],
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
        MessagePrinter::PrintErrorTxt("unsupported calculation type in the DiffusionElmt, please check your related code");
        MessagePrinter::AsFem_Exit();
    }
}
//****************************************************************
void DiffusionElmt::ComputeResidual(const LocalElmtInfo &elmtinfo,
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const Materials &Mate,const Materials &MateOld,
                                 VectorXd &localR) {
    //***********************************************************
    //*** get rid of unused warnings
    //***********************************************************
    if(elmtinfo.dt||soln.gpU.size()||shp.test||Mate.GetVectorMate().size()||MateOld.GetScalarMate().size()){}

    localR(1)=soln.gpV[1]*shp.test+Mate.ScalarMaterials("D")*(soln.gpGradU[1]*shp.grad_test);

}

//****************************************************************************
void DiffusionElmt::ComputeJacobian(const LocalElmtInfo &elmtinfo,const double (&ctan)[2],
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const Materials &Mate,const Materials &MateOld,
                                 MatrixXd &localK) {
    //***********************************************************
    //*** get rid of unused warnings
    //***********************************************************
    if(elmtinfo.dt||soln.gpU.size()||Mate.GetScalarMate().size()||MateOld.GetScalarMate().size()){}

    localK(1,1)=shp.trial*shp.test*ctan[1]+Mate.ScalarMaterials("dDdc")*shp.trial*(soln.gpGradU[1]*shp.grad_test)*ctan[0]
            +Mate.ScalarMaterials("D")*shp.grad_trial*shp.grad_test*ctan[0];

}

//**************************************************************************
void DiffusionElmt::ComputeProjection(const LocalElmtInfo &elmtinfo,const double (&ctan)[2],
                                   const LocalElmtSolution &soln,
                                   const LocalShapeFun &shp,
                                   const Materials &Mate,const Materials &MateOld,
                                   ScalarMateType &gpProj) {
    //***********************************************************
    //*** get rid of unused warnings
    //***********************************************************
    if(elmtinfo.dt||ctan[0]||soln.gpU.size()||shp.test||Mate.GetVectorMate().size()||MateOld.GetVectorMate().size()||gpProj.size()){}

}
