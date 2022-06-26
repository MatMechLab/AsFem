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
//+++ Date   : 2021.11.14
//+++ Purpose: implement the residual and jacobian for general
//+++          stress-diffusion equation
//+++          1) dc/dt=div(D*grad(c)+Dc*omega*grad(sigma_H))
//+++          2) div(sigma)=0
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/StressDiffusionElmt.h"

void StressDiffusionElmt::ComputeAll(const FECalcType &calctype,const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
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
        MessagePrinter::PrintErrorTxt("unsupported calculation type in StressDiffusionElmt, please check your related code");
        MessagePrinter::AsFem_Exit();
    }
}
//*****************************************************************************
void StressDiffusionElmt::ComputeResidual(const LocalElmtInfo &elmtinfo,
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const Materials &Mate,const Materials &MateOld,
                                 VectorXd &localR) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.dt||soln.gpU.size()||shp.test||Mate.GetScalarMate().size()||MateOld.GetScalarMate().size()){}
    // calculate the residual contribution of Mechanics problem
    // For R_c
    D=Mate.ScalarMaterials("D");
    Omega=Mate.ScalarMaterials("Omega");
    GradSigmaH=Mate.VectorMaterials("GradSigmaH");
    localR(1)=soln.gpV[1]*shp.test
        +D*soln.gpGradU[1]*shp.grad_test
        +D*soln.gpU[1]*Omega*GradSigmaH*shp.grad_test;
    //***************************************************
    // For mechanics part
    Stress=Mate.Rank2Materials("stress");
    // For R_ux
    localR(2)=Stress.IthRow(1)*shp.grad_test;
    if(elmtinfo.nDim>=2){
        // For R_uy
        localR(3)=Stress.IthRow(2)*shp.grad_test;
        if(elmtinfo.nDim==3){
            // For R_uz
            localR(4)=Stress.IthRow(3)*shp.grad_test;
        }
    }
}
//******************************************************************************
void StressDiffusionElmt::ComputeJacobian(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const Materials &Mate,const Materials &MateOld,
                                 MatrixXd &localK) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.dt||ctan[0]||soln.gpU.size()||shp.test||Mate.GetScalarMate().size()||MateOld.GetScalarMate().size()){}
    // for the stiffness matrix of mechanics problem
    // For diffusion equation
    D=Mate.ScalarMaterials("D");
    Omega=Mate.ScalarMaterials("Omega");
    GradSigmaH=Mate.VectorMaterials("GradSigmaH");
    dSigmaHdC=Mate.ScalarMaterials("dSigmaHdC");
    // K_c,c
    localK(1,1)=shp.trial*shp.test*ctan[1]
        +D*shp.grad_trial*shp.grad_test*ctan[0]
        +D*shp.trial*Omega*GradSigmaH*shp.grad_test*ctan[0]
        +D*soln.gpU[1]*Omega*dSigmaHdC*shp.grad_trial*shp.grad_test*ctan[0];
    // K_c,ux
    localK(1,2)=0.0;
    if(elmtinfo.nDim>=2) localK(1,3)=0.0;// K_c,uy
    if(elmtinfo.nDim==3) localK(1,4)=0.0;// K_c,uz

    //*********************************************
    //*** for stress equilibrium equation
    //*********************************************
    // K_ux,ux
    localK(2,2)=Mate.Rank4Materials("jacobian").GetIKjlComponent(1,1,shp.grad_test,shp.grad_trial)*ctan[0];
    // K_ux,c
    localK(2,1)=dStressdC.IthRow(1)*shp.grad_test*shp.trial*ctan[0];
    if(elmtinfo.nDim>=2){
        // K_ux,uy
        localK(2,3)=Mate.Rank4Materials("jacobian").GetIKjlComponent(1,2,shp.grad_test,shp.grad_trial)*ctan[0];
        //*************************
        // K_uy,c
        localK(3,1)=dStressdC.IthRow(2)*shp.grad_test*shp.trial*ctan[0];
        // K_uy,ux
        localK(3,2)=Mate.Rank4Materials("jacobian").GetIKjlComponent(2,1,shp.grad_test,shp.grad_trial)*ctan[0];
        // K_uy,uy
        localK(3,3)=Mate.Rank4Materials("jacobian").GetIKjlComponent(2,2,shp.grad_test,shp.grad_trial)*ctan[0];
        if(elmtinfo.nDim==3){
            // K_ux,uz
            localK(2,4)=Mate.Rank4Materials("jacobian").GetIKjlComponent(1,3,shp.grad_test,shp.grad_trial)*ctan[0];
            // K_uy,uz
            localK(3,4)=Mate.Rank4Materials("jacobian").GetIKjlComponent(2,3,shp.grad_test,shp.grad_trial)*ctan[0];
            // K_uz,c
            localK(4,1)=dStressdC.IthRow(3)*shp.grad_test*shp.trial*ctan[0];
            // K_uz,ux
            localK(4,2)=Mate.Rank4Materials("jacobian").GetIKjlComponent(3,1,shp.grad_test,shp.grad_trial)*ctan[0];
            // K_uz,uy
            localK(4,3)=Mate.Rank4Materials("jacobian").GetIKjlComponent(3,2,shp.grad_test,shp.grad_trial)*ctan[0];
            // K_uz,uz
            localK(4,4)=Mate.Rank4Materials("jacobian").GetIKjlComponent(3,3,shp.grad_test,shp.grad_trial)*ctan[0];
        }
    }
}
//*************************************************
void StressDiffusionElmt::ComputeProjection(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                   const LocalElmtSolution &soln,
                                   const LocalShapeFun &shp,
                                   const Materials &Mate,const Materials &MateOld,
                                   ScalarMateType &gpProj) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.dt||ctan[0]||soln.gpU.size()||shp.test||
       Mate.GetScalarMate().size()||MateOld.GetScalarMate().size()||
       gpProj.size()){}

}
