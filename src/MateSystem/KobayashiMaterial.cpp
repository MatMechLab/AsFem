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
//+++ Date   : 2021.11.12
//+++ Purpose: Calculate the material properties required by Kobayashi
//+++          element. In this code, we can define:
//+++           1) K
//+++           2) dK
//+++           3) dKdGradEta
//+++           4) ddKdGradEta
//+++ Ref.   :  https://doi.org/10.1016/0167-2789(93)90120-P
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/KobayashiMaterial.h"

void KobayashiMaterial::InitMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, Materials &Mate){
    //***************************************************
    //*** get rid of unused warning
    //***************************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU[0]||Mate.GetScalarMate().size()){}

}

//********************************************************************
void KobayashiMaterial::ComputeMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, const Materials &MateOld, Materials &Mate){
    //**************************************************************
    //*** get rid of unused warning
    //**************************************************************
    if(InputParams.size()||elmtinfo.dt||elmtsoln.gpU[0]||MateOld.GetScalarMate().size()||Mate.GetScalarMate().size()){}
    // parameter size check
    if(InputParams.size()<4){
        MessagePrinter::PrintErrorTxt("for kobayashi material, four parameters are required, you should give: L, k0, delta, N");
        MessagePrinter::AsFem_Exit();
    }

    // K=K0*[1+delta*cos(N*theta)]
    Mate.ScalarMaterials("L")=InputParams[1-1];
    K0=InputParams[2-1];
    delta=InputParams[3-1];
    N=InputParams[4-1];
    Mate.ScalarMaterials("Latent")=InputParams[5-1];
    if(N<=1.0){
        MessagePrinter::PrintErrorTxt("the mode number must be a positive integer");
        MessagePrinter::AsFem_Exit();
    }
    // Dofs: 1-->eta (order parameter), 2-->T (temperature)
    eta=elmtsoln.gpU[1];
    T=elmtsoln.gpU[2];

    const double PI=3.141592653589793238462;
    // Now we start to calculate the double well free energy density(T-dependent!)
    m=0.9*atan(10.0*T)/PI;
    dmdT=0.09/(PI*T*T+0.01*PI);
    dfdeta=0.5*(eta-1.0)*eta*(2*eta+2*m-1.0);
    d2fdeta2=m*(2*eta-1.0)+3*eta*eta-3*eta+0.5;
    d2fdetadT=0.5*(eta-1.0)*eta*2*dmdT;
    // store them into scalar materials
    Mate.ScalarMaterials("dFdeta")=dfdeta;
    Mate.ScalarMaterials("d2Fdeta2")=d2fdeta2;
    Mate.ScalarMaterials("d2FdetadT")=d2fdetadT;

    // Now we start to calculate the anisotropic kappa
    const double tol=1.0e-9;
    double threshold;
    double norm,normsq,n;

    GradEta=elmtsoln.gpGradU[1];
    threshold=1.0-tol;
    normsq=GradEta(1)*GradEta(1)+GradEta(2)*GradEta(2);
    norm=sqrt(normsq);

    n=0.0;
    if(normsq>tol){
        n=GradEta(1)/norm;
    }

    if(n> threshold) n= threshold;
    if(n<-threshold) n=-threshold;

    double theta,dthetadn;
    Vector3d dndgradeta;

    theta   = acos(n)*Sign(GradEta(2));
    dthetadn=-Sign(GradEta(2))/sqrt(1.0-n*n);

    dndgradeta=0.0;
    if(normsq>tol){
        dndgradeta(1)= GradEta(2)*GradEta(2)/(norm*normsq);
        dndgradeta(2)=-GradEta(1)*GradEta(2)/(norm*normsq);
        dndgradeta(3)=0.0;
    }

    K  = K0*(1.0+delta*cos(N*theta));
    dK =-K0*delta*N*sin(N*theta);
    ddK=-K0*delta*N*N*cos(N*theta);

    dKdGradEta=dK*dthetadn*dndgradeta;
    ddKdGradEta=ddK*dthetadn*dndgradeta;

    // store the variables into different materials
    Mate.ScalarMaterials("K")=K;
    Mate.ScalarMaterials("dK")=dK;
    Mate.VectorMaterials("dKdGradEta")=dKdGradEta;
    Mate.VectorMaterials("ddKdGradEta")=ddKdGradEta;

    Mate.VectorMaterials("GradEta")=GradEta;

}
