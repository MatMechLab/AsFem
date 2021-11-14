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
//+++ Date   : 2021.11.14
//+++ Purpose: implement the residual and jacobian for general
//+++          stress-diffusion-fracture equation
//+++          1) dc/dt=div(D*grad(c)+Dc*omega*grad(sigma_H))
//+++          2) dd/dt=-delta F/delta d
//+++          3) div(sigma)=0
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/DiffusionFractureElmt.h"

void DiffusionFractureElmt::ComputeAll(const FECalcType &calctype,const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
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
        MessagePrinter::PrintErrorTxt("unsupported calculation type in DiffusionFractureElmt, please check your related code");
        MessagePrinter::AsFem_Exit();
    }
}
//*****************************************************************************
void DiffusionFractureElmt::ComputeResidual(const LocalElmtInfo &elmtinfo,
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const Materials &Mate,const Materials &MateOld,
                                 VectorXd &localR) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.dt||soln.gpU.size()||shp.test||Mate.GetScalarMate().size()||MateOld.GetScalarMate().size()){}
    //*************************************************
    // for diffusion equation
    //*************************************************
    D=Mate.ScalarMaterials("D");
    Omega=Mate.ScalarMaterials("Omega");
    GradSigmaH=Mate.VectorMaterials("GradSigmaH");
    // R_c
    localR(1)=soln.gpV[1]*shp.test
        +D*soln.gpGradU[1]*shp.grad_test
        +D*soln.gpU[1]*Omega*GradSigmaH*shp.grad_test;
    //***************************************************
    // for damage field
    //***************************************************
    viscosity=Mate.ScalarMaterials("viscosity");
    Gc=Mate.ScalarMaterials("Gc");
    L=Mate.ScalarMaterials("L");
    Hist=Mate.ScalarMaterials("H");
    localR(2)=viscosity*soln.gpV[2]*shp.test
        +2*(soln.gpU[2]-1)*Hist*shp.test
        +(Gc/L)*soln.gpU[2]*shp.test
        +Gc*L*(soln.gpGradU[2]*shp.grad_test);
    //***************************************************
    // For mechanics part
    Stress=Mate.Rank2Materials("stress")-MateOld.Rank2Materials("stress");
    localR(3)=Stress.IthRow(1)*shp.grad_test;
    if(elmtinfo.nDim>=2){
        localR(4)=Stress.IthRow(2)*shp.grad_test;
        if(elmtinfo.nDim==3){
            localR(5)=Stress.IthRow(3)*shp.grad_test;
        }
    }
}
//******************************************************************************
void DiffusionFractureElmt::ComputeJacobian(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                 const LocalElmtSolution &soln,
                                 const LocalShapeFun &shp,
                                 const Materials &Mate,const Materials &MateOld,
                                 MatrixXd &localK) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.dt||ctan[0]||soln.gpU.size()||shp.test||Mate.GetScalarMate().size()||MateOld.GetScalarMate().size()){}
    //***************************************
    // for diffusion equation
    //***************************************
    D=Mate.ScalarMaterials("D");
    Omega=Mate.ScalarMaterials("Omega");
    GradSigmaH=Mate.VectorMaterials("GradSigmaH");
    dSigmaHdC=Mate.ScalarMaterials("dSigmaHdC");
    // K_c,c
    localK(1,1)=shp.trial*shp.test*ctan[1]
        +D*shp.grad_trial*shp.grad_test*ctan[0]
        +D*shp.trial*Omega*GradSigmaH*shp.grad_test*ctan[0]
        +D*soln.gpU[1]*Omega*dSigmaHdC*shp.grad_trial*shp.grad_test*ctan[0];
    // K_c,d
    localK(1,2)=0.0;
    // K_c,ux
    localK(1,3)=0.0;
    if(elmtinfo.nDim==2) localK(1,4)=0.0;// K_c,uy
    if(elmtinfo.nDim==3) localK(1,5)=0.0;// K_c,uz

    //*********************************************
    // for damage field
    //*********************************************
    viscosity=Mate.ScalarMaterials("viscosity");
    Gc=Mate.ScalarMaterials("Gc");
    L=Mate.ScalarMaterials("L");
    Hist=Mate.ScalarMaterials("H");
    dStressdD=Mate.Rank2Materials("dStressdD");
    dStressdC=Mate.Rank2Materials("dStressdC");
    dHdstrain=Mate.Rank2Materials("dHdstrain");
    // K_d,d 
    localK(2,2)=viscosity*shp.trial*shp.test*ctan[1]
               +2*shp.trial*Hist*shp.test*ctan[0]
               +(Gc/L)*shp.trial*shp.test*ctan[0]
               +Gc*L*(shp.grad_trial*shp.grad_test)*ctan[0];
    //*******************************
    valx=0.0;valy=0.0;valz=0.0;
    for(int i=1;i<=3;i++){
        valx+=0.5*(dHdstrain(1,i)+dHdstrain(i,1))*shp.grad_trial(i);
        valy+=0.5*(dHdstrain(2,i)+dHdstrain(i,2))*shp.grad_trial(i);
        valz+=0.5*(dHdstrain(3,i)+dHdstrain(i,3))*shp.grad_trial(i);
    }
    // K_d,c
    localK(2,1)=0.0;
    // K_d,ux
    localK(2,3)=2*(soln.gpU[2]-1)*valx*shp.test*ctan[0];
    // K_d,uy
    localK(2,4)=2*(soln.gpU[2]-1)*valy*shp.test*ctan[0];
    if(elmtinfo.nDim==3){
        // K_d,uz
        localK(2,5)=2*(soln.gpU[2]-1)*valz*shp.test*ctan[0];
    }
    //*********************************************
    //*** for stress equilibrium equation
    //*********************************************
    // K_ux,ux
    localK(3,3)=Mate.Rank4Materials("jacobian").GetIKjlComponent(1,1,shp.grad_test,shp.grad_trial)*ctan[0];
    // K_ux,c
    localK(3,1)=dStressdC.IthRow(1)*shp.grad_test*shp.trial*ctan[0];
    // K_ux,d
    localK(3,2)=dStressdD.IthRow(1)*shp.grad_test*shp.trial*ctan[0];
    if(elmtinfo.nDim>=2){
        // K_ux,uy
        localK(3,4)=Mate.Rank4Materials("jacobian").GetIKjlComponent(1,2,shp.grad_test,shp.grad_trial)*ctan[0];
        //*************************
        // K_uy,c
        localK(4,1)=dStressdC.IthRow(2)*shp.grad_test*shp.trial*ctan[0];
        // K_uy,d
        localK(4,2)=dStressdD.IthRow(2)*shp.grad_test*shp.trial*ctan[0];
        // K_uy,ux
        localK(4,3)=Mate.Rank4Materials("jacobian").GetIKjlComponent(2,1,shp.grad_test,shp.grad_trial)*ctan[0];
        // K_uy,uy
        localK(4,4)=Mate.Rank4Materials("jacobian").GetIKjlComponent(2,2,shp.grad_test,shp.grad_trial)*ctan[0];
        if(elmtinfo.nDim==3){
            // K_ux,uz
            localK(3,5)=Mate.Rank4Materials("jacobian").GetIKjlComponent(1,3,shp.grad_test,shp.grad_trial)*ctan[0];
            // K_uy,uz
            localK(4,5)=Mate.Rank4Materials("jacobian").GetIKjlComponent(2,3,shp.grad_test,shp.grad_trial)*ctan[0];
            // K_uz,c
            localK(5,1)=dStressdC.IthRow(3)*shp.grad_test*shp.trial*ctan[0];
            // K_uz,d
            localK(5,2)=dStressdD.IthRow(3)*shp.grad_test*shp.trial*ctan[0];
            // K_uz,ux
            localK(5,3)=Mate.Rank4Materials("jacobian").GetIKjlComponent(3,1,shp.grad_test,shp.grad_trial)*ctan[0];
            // K_uz,uy
            localK(5,4)=Mate.Rank4Materials("jacobian").GetIKjlComponent(3,2,shp.grad_test,shp.grad_trial)*ctan[0];
            // K_uz,uz
            localK(5,5)=Mate.Rank4Materials("jacobian").GetIKjlComponent(3,3,shp.grad_test,shp.grad_trial)*ctan[0];
        }
    }
}
//*************************************************
void DiffusionFractureElmt::ComputeProjection(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                   const LocalElmtSolution &soln,
                                   const LocalShapeFun &shp,
                                   const Materials &Mate,const Materials &MateOld,
                                   ScalarMateType &gpProj) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.dt||ctan[0]||soln.gpU.size()||shp.test||Mate.GetScalarMate().size()||MateOld.GetScalarMate().size()){}

    gpProj["reacforce_x"]=Mate.Rank2Materials("stress").IthRow(1)*shp.grad_test;
    gpProj["reacforce_y"]=Mate.Rank2Materials("stress").IthRow(2)*shp.grad_test;
    gpProj["reacforce_z"]=Mate.Rank2Materials("stress").IthRow(3)*shp.grad_test;
}
