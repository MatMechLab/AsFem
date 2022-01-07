//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2022.01.07
//+++ Purpose: implement the residual and jacobian for wave equation
//+++          original one is: d2u/dt2=c*c*lap(u)+f
//+++          here we use:
//+++            1) dv/dt=c*c*lap(u)+f
//+++            2) du/dt=v
//+++          1st dof: v
//+++          2nd dof: u
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#include "ElmtSystem/WaveElmt.h"

void WaveElmt::ComputeAll(const FECalcType &calctype,const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
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
        MessagePrinter::PrintErrorTxt("unsupported calculation type in WaveElmt, please check your related code");
        MessagePrinter::AsFem_Exit();
    }
}
//***************************************************************************
void WaveElmt::ComputeResidual(const LocalElmtInfo &elmtinfo,
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const Materials &Mate,const Materials &MateOld,
                                 VectorXd &localR) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.dt||soln.gpU[0]||shp.test||Mate.GetScalarMate().size()||MateOld.GetScalarMate().size()) {}

    double c=Mate.ScalarMaterials("C");
    double f=Mate.ScalarMaterials("f");

    // for R_v
    localR(1)=soln.gpV[1]*shp.test
        +c*c*(soln.gpGradU[2]*shp.grad_test)
        -f*shp.test;

    // for R_u
    localR(2)=soln.gpV[2]*shp.test
        -soln.gpU[1]*shp.test;

}
//*****************************************************************************
void WaveElmt::ComputeJacobian(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const Materials &Mate,const Materials &MateOld,
                                 MatrixXd &localK) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.dt||ctan[0]||soln.gpU[0]||Mate.GetScalarMate().size()||MateOld.GetScalarMate().size()){}

    double c=Mate.ScalarMaterials("C");
    double dfdu=Mate.ScalarMaterials("dfdu");
    double dfdv=Mate.ScalarMaterials("dfdv");

    // for R_v,v
    localK(1,1)= shp.trial*shp.test*ctan[1]
        -dfdv*shp.trial*shp.test*ctan[0];
    // for R_v,u
    localK(1,2)= c*c*(shp.grad_trial*shp.grad_test)*ctan[0]
        -dfdu*shp.trial*shp.test*ctan[0];

    // for R_u,u
    localK(2,2)= shp.trial*shp.test*ctan[1];
    //for R_u,v
    localK(2,1)=-shp.trial*shp.test*ctan[0];
}
//*******************************************************************************
void WaveElmt::ComputeProjection(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                   const LocalElmtSolution &soln,
                                   const LocalShapeFun &shp,
                                   const Materials &Mate,const Materials &MateOld,
                                   ScalarMateType &gpProj) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.dt||ctan[0]||soln.gpU[0]||shp.test||Mate.GetScalarMate().size()||MateOld.GetScalarMate().size()||gpProj.size()){}
}
