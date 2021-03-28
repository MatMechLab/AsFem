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
//+++ Date   : 2021.01.18
//+++ Purpose: implement the residual and jacobian for User-Defined-element
//+++          1) dc/dt=div(D*grad(c))
//+++          2) div(\Sigma)=0
//+++          this is used for testing purpose
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/BulkElmtSystem.h"

void BulkElmtSystem::User1Elmt(const FECalcType &calctype, const int &nDim, const int &nNodes, const int &nDofs,
                               const double &t, const double &dt, const double (&ctan)[2], const Vector3d &gpCoords,
                               const vector<double> &gpU, const vector<double> &gpV, const vector<Vector3d> &gpGradU,
                               const vector<Vector3d> &gpGradV, const double &test, const double &trial,
                               const Vector3d &grad_test, const Vector3d &grad_trial,
                               const ScalarMateType &ScalarMaterials, const VectorMateType &VectorMaterials,
                               const Rank2MateType &Rank2Materials, const Rank4MateType &Rank4Materials,
                               vector<double> &gpHist, vector<double> &gpHistOld,map<string,double> &gpProj,
                               MatrixXd &localK, VectorXd &localR) {
    //*******************************************************
    //*** to get rid of the warning for unused variables  ***
    //*** for normal users, you dont need to do this      ***
    //*******************************************************
    if(nDim||nNodes||nDofs||t||dt||ctan[0]||gpCoords(1)||gpU.size()||gpV.size()||
       gpGradU.size()||gpGradV.size()||test||trial||grad_test(1)||grad_trial(1)||
       ScalarMaterials.size()||VectorMaterials.size()||Rank2Materials.size()||Rank4Materials.size()||
       gpHist.size()||gpHistOld.size()||gpProj.size()){}

    // Dofs: C -->1
    //       ux-->2
    //       uy-->3
    //       uz-->4
    double D=ScalarMaterials.at("D");
    switch (calctype){
        case FECalcType::ComputeResidual:
            // For R_d
            localR(1)=gpV[1]*test+D*(gpGradU[1]*grad_test);
            // For R_ux
            localR(2)=Rank2Materials.at("stress").IthRow(1)*grad_test;
            // For R_uy
            localR(3)=Rank2Materials.at("stress").IthRow(2)*grad_test;
            if(nDim==3){
                // For R_uz
                localR(4)=Rank2Materials.at("stress").IthRow(3)*grad_test;
            }
            break;
        case FECalcType::ComputeJacobian:
            // K_c,c
            localK(1,1)=trial*test*ctan[1]+D*(grad_trial*grad_test)*ctan[0];
            localK(1,2)=0.0;localK(1,3)=0.0;
            //*******************************
            // K_ux,c
            localK(2,1)=0.0;
            // K_ux,ux
            localK(2,2)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(1,1,grad_test,grad_trial)*ctan[0];
            // K_ux,uy
            localK(2,3)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(1,2,grad_test,grad_trial)*ctan[0];

            // K_uy,c
            localK(3,1)=0.0;
            // K_uy,uy
            localK(3,3)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(2,2,grad_test,grad_trial)*ctan[0];
            // K_uy,ux
            localK(3,2)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(2,1,grad_test,grad_trial)*ctan[0];
            if(nDim==3){
                // K_ux,uz
                localK(2,4)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(1,3,grad_test,grad_trial)*ctan[0];

                // K_uy,uz
                localK(3,4)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(2,3,grad_test,grad_trial)*ctan[0];

                // K_uz,c
                localK(4,1)=0.0;
                // K_uz,uz
                localK(4,4)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(3,3,grad_test,grad_trial)*ctan[0];
                // K_uz,ux
                localK(4,2)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(3,1,grad_test,grad_trial)*ctan[0];
                // K_uz,uy
                localK(4,3)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(3,2,grad_test,grad_trial)*ctan[0];
            }
            break;
        case FECalcType::InitHistoryVariable:
            fill(gpHist.begin(),gpHist.end(),0.0);
            break;
        case FECalcType::UpdateHistoryVariable:
            gpHistOld=gpHist;
            break;
        case FECalcType::Projection:
            break;
        default:
            MessagePrinter::PrintErrorTxt("unsupported FEM calculation type in User1 element");
            MessagePrinter::AsFem_Exit();
            break;
    }
}