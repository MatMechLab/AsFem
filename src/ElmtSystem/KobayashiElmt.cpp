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
//+++ Date   : 2021.11.12
//+++ Purpose: implement the residual and jacobian for Kobayashi's
//+++          dendrite model
//+++          the governing equations are:
//+++          1) deta/dt=L*div(k*dk*v)+L*div(k*k*\nabla\eta)-L*df/deta
//+++          2) dT/dt=Lap(T)+K*deta/dt             
//+++ Dofs   : 1. eta, 2, T
//+++ Ref    : https://doi.org/10.1016/0167-2789(93)90120-P
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#include "ElmtSystem/KobayashiElmt.h"

void KobayashiElmt::ComputeAll(const FECalcType &calctype,const LocalElmtInfo &elmtinfo,
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
        MessagePrinter::PrintErrorTxt("unsupported calculation type in KobayashiElmt, please check your related code");
        MessagePrinter::AsFem_Exit();
    }
}
//*******************************************************************
void KobayashiElmt::ComputeResidual(const LocalElmtInfo &elmtinfo,
                                    const LocalElmtSolution &soln,
                                    const LocalShapeFun &shp,
                                    const Materials &Mate,const Materials &MateOld,
                                    VectorXd &localR) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.nDim||soln.gpU.size()||shp.test||Mate.GetScalarMate().size()||MateOld.GetScalarMate().size()){}
    L=Mate.ScalarMaterials("L");
    K=Mate.ScalarMaterials("K");
    dK=Mate.ScalarMaterials("dK");
    dFdeta=Mate.ScalarMaterials("dFdeta");
    Latent=Mate.ScalarMaterials("Latent");
    // For R_eta
    V(1)=-soln.gpGradU[1](2);
    V(2)= soln.gpGradU[1](1);
    V(3)= 0.0;
    localR(1)=soln.gpV[1]*shp.test
        +L*K*dK*(V*shp.grad_test)
        +L*K*K*(soln.gpGradU[1]*shp.grad_test)
        +L*dFdeta*shp.test;
    // For R_T
    localR(2)=soln.gpV[2]*shp.test
        +soln.gpGradU[2]*shp.grad_test
        -Latent*soln.gpV[1]*shp.test;

}
//***************************************************************************
void KobayashiElmt::ComputeJacobian(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                    const LocalElmtSolution &soln,
                                    const LocalShapeFun &shp,
                                    const Materials &Mate,const Materials &MateOld,
                                    MatrixXd &localK) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.dt||soln.gpU.size()||shp.test||Mate.GetVectorMate().size()||MateOld.GetVectorMate().size()){}
    //************************************************
    // take the parameters from Material properties
    //************************************************
    L=Mate.ScalarMaterials("L");
    K=Mate.ScalarMaterials("K");
    dK=Mate.ScalarMaterials("dK");
    Latent=Mate.ScalarMaterials("Latent");
    d2Fdeta2=Mate.ScalarMaterials("d2Fdeta2");
    d2FdetadT=Mate.ScalarMaterials("d2FdetadT");
    dKdGradEta=Mate.VectorMaterials("dKdGradEta");
    ddKdGradEta=Mate.VectorMaterials("ddKdGradEta");
    //************************************************
    V(1)=-soln.gpGradU[1](2);
    V(2)= soln.gpGradU[1](1);
    V(3)= 0.0;
    //************************
    dV(1)=-shp.grad_trial(2);
    dV(2)= shp.grad_trial(1);
    dV(3)= 0.0;
    // K_eta,eta
    localK(1,1)=shp.trial*shp.test*ctan[1]
        +L*(dKdGradEta*shp.grad_trial)*dK*(V*shp.grad_test)*ctan[0]
        +L*K*(ddKdGradEta*shp.grad_trial)*(V*shp.grad_test)*ctan[0]
        +L*K*dK*(dV*shp.grad_test)*ctan[0]
        +L*2*K*(dKdGradEta*shp.grad_trial)*(soln.gpGradU[1]*shp.grad_test)*ctan[0]
        +L*K*K*(shp.grad_trial*shp.grad_test)*ctan[0]
        +L*d2Fdeta2*shp.trial*shp.test*ctan[0];
    // K_eta,T
    localK(1,2)=L*d2FdetadT*shp.trial*shp.test*ctan[0];
    //*************************************************************
    // K_T,eta
    localK(2,1)=-Latent*shp.trial*shp.test*ctan[1];
    // K_T,T
    localK(2,2)=shp.trial*shp.test*ctan[1]
        +shp.grad_trial*shp.grad_test*ctan[0];
}
//***********************************************************
void KobayashiElmt::ComputeProjection(const LocalElmtInfo &elmtinfo,const double (&ctan)[3],
                                         const LocalElmtSolution &soln,
                                         const LocalShapeFun &shp,
                                         const Materials &Mate,const Materials &MateOld,
                                         ScalarMateType &gpProj) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(elmtinfo.dt||ctan[0]||soln.gpU.size()||shp.test||Mate.GetScalarMate().size()||MateOld.GetScalarMate().size()||gpProj.size()){}
}
