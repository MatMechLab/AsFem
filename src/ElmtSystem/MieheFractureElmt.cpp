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
//+++ Date   : 2021.01.17
//+++ Purpose: implement the residual and jacobian for Miehe's
//+++          phase field fracture model
//+++          1) eta*dD/dt=2(1-D)H-(Gc/L)D+Gc*L*Lap(D) (see Eq. 47)
//+++          2) div(\Sigma)=0
//+++ Reference: A phase field model for rate-independent crack propagation:
//+++            Robust algorithmic implementation based on operator splits
//+++ DOI    : https://doi.org/10.1016/j.cma.2010.04.011
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/BulkElmtSystem.h"

void BulkElmtSystem::MieheFractureElmt(const FECalcType &calctype,
                                      const int &nDim, const int &nNodes,const int &nDofs,
                                      const double &t,const double &dt, const double (&ctan)[2],
                                      const Vector3d &gpCoords,
                                      const vector<double> &gpU, const vector<double> &gpV,
                                      const vector<Vector3d> &gpGradU, const vector<Vector3d> &gpGradV,
                                      const double &test, const double &trial, const Vector3d &grad_test,
                                      const Vector3d &grad_trial, const ScalarMateType &ScalarMaterials,
                                      const VectorMateType &VectorMaterials, const Rank2MateType &Rank2Materials,
                                      const Rank4MateType &Rank4Materials, vector<double> &gpHist,
                                      vector<double> &gpHistOld, vector<double> &gpProj, MatrixXd &localK,
                                      VectorXd &localR) {
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
    int i;
    double val;
    switch (calctype){
        case FECalcType::ComputeResidual:
            // For R_d
            localR(1)=ScalarMaterials.at("viscosity")*gpV[1]*test
                       +2*(gpU[1]-1)*ScalarMaterials.at("Hist")*test
                       +(ScalarMaterials.at("Gc")/ScalarMaterials.at("L"))*gpU[1]*test
                       +ScalarMaterials.at("Gc")*ScalarMaterials.at("L")*gpGradU[1]*grad_test;
            // For R_ux
            localR(2)=Rank2Materials.at("Stress").IthRow(1)*grad_test;
            // For R_uy
            localR(3)=Rank2Materials.at("Stress").IthRow(2)*grad_test;
            if(nDim==3){
                // For R_uz
                localR(4)=Rank2Materials.at("Stress").IthRow(3)*grad_test;
            }
            break;
        case FECalcType::ComputeJacobian:
            // K_d,d  (see Eq. 47)
            localK(1,1)=ScalarMaterials.at("viscosity")*trial*test*ctan[1]
                           +2*trial*ScalarMaterials.at("Hist")*test*ctan[0]
                           +(ScalarMaterials.at("Gc")/ScalarMaterials.at("L"))*trial*test*ctan[0]
                           +ScalarMaterials.at("Gc")*ScalarMaterials.at("L")*grad_trial*grad_test*ctan[0];
            // K_d,ux
            val=0.0;
            for(i=1;i<=nDim;i++){
                val+=0.5*(Rank2Materials.at("dHdstrain")(1,i)+Rank2Materials.at("dHdstrain")(i,1))*grad_trial(i);
            }
            localK(1,2)=2*(gpU[1]-1)*val*test;

            // K_d,uy
            val=0.0;
            for(i=1;i<=nDim;i++){
                val+=0.5*(Rank2Materials.at("dHdstrain")(2,i)+Rank2Materials.at("dHdstrain")(i,2))*grad_trial(i);
            }
            localK(1,3)=2*(gpU[1]-1)*val*test;

            if(nDim==3){
                // K_d,uz
                val=0.0;
                for(i=1;i<=nDim;i++){
                    val+=0.5*(Rank2Materials.at("dHdstrain")(3,i)+Rank2Materials.at("dHdstrain")(i,3))*grad_trial(i);
                }
                localK(1,4)=2*(gpU[1]-1)*val*test;
            }

            // K_ux,ux
            localK(2,2)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(1,1,grad_test,grad_trial)*ctan[0];
            // K_ux,uy
            localK(2,3)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(1,2,grad_test,grad_trial)*ctan[0];
            // K_ux,d
            localK(2,1)=Rank2Materials.at("dStressdD").IthRow(1)*grad_test*trial*ctan[0];

            // K_uy,uy
            localK(3,3)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(2,2,grad_test,grad_trial)*ctan[0];
            // K_uy,ux
            localK(3,2)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(2,1,grad_test,grad_trial)*ctan[0];
            // K_uy,d
            localK(3,1)=Rank2Materials.at("dStressdD").IthRow(2)*grad_test*trial*ctan[0];
            if(nDim==3){
                // K_ux,uz
                localK(2,4)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(1,3,grad_test,grad_trial)*ctan[0];

                // K_uy,uz
                localK(3,4)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(2,4,grad_test,grad_trial)*ctan[0];

                // K_uz,uz
                localK(4,4)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(3,3,grad_test,grad_trial)*ctan[0];
                // K_uz,ux
                localK(4,2)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(3,1,grad_test,grad_trial)*ctan[0];
                // K_uz,uy
                localK(4,3)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(3,2,grad_test,grad_trial)*ctan[0];
                // K_uz,d
                localK(4,1)=Rank2Materials.at("dStressdD").IthRow(3)*grad_test*trial*ctan[0];
            }
            break;
        case FECalcType::InitHistoryVariable:
            gpHist[0]=0.0;
            break;
        case FECalcType::UpdateHistoryVariable:
            gpHistOld=gpHist;
            break;
        case FECalcType::Projection:
            gpProj[0]=ScalarMaterials.at("vonMises");
            gpProj[1]=gpGradU[1](1);
            gpProj[2]=gpGradU[1](2);
            gpProj[3]=gpGradU[1](3);
            break;
        default:
            MessagePrinter::PrintErrorTxt("unsupported FEM calculation type in Miehe fracture element");
            MessagePrinter::AsFem_Exit();
            break;
    }
}