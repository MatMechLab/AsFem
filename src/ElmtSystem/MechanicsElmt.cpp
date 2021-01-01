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
//+++ Date   : 2020.12.31
//+++ Purpose: implement the residual and jacobian for general
//+++          stress equilibrium equation
//+++          div(sigma)=0
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/BulkElmtSystem.h"

void BulkElmtSystem::MechanicsElmt(const FECalcType &calctype, const int &nDim, const int &nNodes, const double &t,
                                   const double &dt, const double (&ctan)[2], const Vector3d &gpCoords,
                                   const vector<double> &gpU, const vector<double> &gpV,
                                   const vector<Vector3d> &gpGradU, const vector<Vector3d> &gpGradV, const double &test,
                                   const double &trial, const Vector3d &grad_test, const Vector3d &grad_trial,
                                   const ScalarMateType &ScalarMaterials, const VectorMateType &VectorMaterials,
                                   const Rank2MateType &Rank2Materials, const Rank4MateType &Rank4Materials,
                                   vector<double> &gpHist, vector<double> &gpHistOld, vector<double> &gpProj,
                                   MatrixXd &localK, VectorXd &localR) {
    //*******************************************************
    //*** to get rid of the warning for unused variables  ***
    //*** for normal users, you dont need to do this      ***
    //*******************************************************
    if(nDim||nNodes||t||dt||ctan[0]||gpCoords(1)||gpU.size()||gpV.size()||
       gpGradU.size()||gpGradV.size()||test||trial||grad_test(1)||grad_trial(1)||
       ScalarMaterials.size()||VectorMaterials.size()||Rank2Materials.size()||Rank4Materials.size()||
       gpHist.size()||gpHistOld.size()||gpProj.size()){}

    // Dofs: C -->1
    //       Mu-->2
    switch (calctype){
        case FECalcType::ComputeResidual:
            // For R_ux
            localR(1)=Rank2Materials.at("stress").IthRow(1)*grad_test;
            if(nDim>=2){
                // for R_uy
                localR(2)=Rank2Materials.at("stress").IthRow(2)*grad_test;
                if(nDim==3){
                    localR(3)=Rank2Materials.at("stress").IthRow(3)*grad_test;
                }
            }
            break;
        case FECalcType::ComputeJacobian:
            // K_ux,ux
            localK(1,1)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(1,1,grad_test,grad_trial)*ctan[0];
            if(nDim>=2){
                // K_ux,uy
                localK(1,2)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(1,2,grad_test,grad_trial)*ctan[0];

                // K_uy,ux
                localK(2,1)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(2,1,grad_test,grad_trial)*ctan[0];
                // K_uy,uy
                localK(2,2)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(2,2,grad_test,grad_trial)*ctan[0];
                if(nDim==3){
                    // K_ux,uz
                    localK(1,3)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(1,3,grad_test,grad_trial)*ctan[0];

                    // K_uy,uz
                    localK(2,3)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(2,3,grad_test,grad_trial)*ctan[0];

                    // K_uz,ux
                    localK(3,1)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(3,1,grad_test,grad_trial)*ctan[0];
                    // K_uz,uy
                    localK(3,2)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(3,2,grad_test,grad_trial)*ctan[0];
                    // K_uz,uz
                    localK(3,3)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(3,3,grad_test,grad_trial)*ctan[0];
                }
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
            gpProj[1]=Rank2Materials.at("stress")(1,1);//sigma_xx
            if(nDim==2){
                gpProj[2]=Rank2Materials.at("stress")(2,2);//sigma_yy
                gpProj[3]=Rank2Materials.at("stress")(1,2);//sigma_xy

                gpProj[4]=Rank2Materials.at("strain")(1,1);//epsilon_xx
                gpProj[5]=Rank2Materials.at("strain")(2,2);//epsilon_yy
                gpProj[6]=Rank2Materials.at("strain")(1,2);//epsilon_xy
            }
            else if(nDim==3){
                gpProj[2]=Rank2Materials.at("stress")(2,2);//sigma_yy
                gpProj[3]=Rank2Materials.at("stress")(3,3);//sigma_zz
                gpProj[4]=Rank2Materials.at("stress")(2,3);//sigma_yz
                gpProj[5]=Rank2Materials.at("stress")(1,3);//sigma_xz
                gpProj[6]=Rank2Materials.at("stress")(1,2);//sigma_xy

                gpProj[7] =Rank2Materials.at("strain")(1,1);//epsilon_xx
                gpProj[8] =Rank2Materials.at("strain")(2,2);//epsilon_yy
                gpProj[9] =Rank2Materials.at("strain")(3,3);//epsilon_zz
                gpProj[10]=Rank2Materials.at("strain")(2,3);//epsilon_yz
                gpProj[11]=Rank2Materials.at("strain")(1,3);//epsilon_xz
                gpProj[12]=Rank2Materials.at("strain")(1,2);//epsilon_xy
            }
            break;
        default:
            MessagePrinter::PrintErrorTxt("unsupported FEM calculation type in Mechanics element");
            MessagePrinter::AsFem_Exit();
            break;
    }
}