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
                                      vector<double> &gpHistOld,map<string,double> &gpProj, MatrixXd &localK,
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
    double valx,valy,valz;
    double viscosity=ScalarMaterials.at("viscosity");
    double Gc=ScalarMaterials.at("Gc");
    double L=ScalarMaterials.at("L");
    double Hist=ScalarMaterials.at("Hist");
    RankTwoTensor Stress=Rank2Materials.at("Stress");
    RankTwoTensor dStressdD=Rank2Materials.at("dStressdD");
    RankTwoTensor dHdstrain=Rank2Materials.at("dHdstrain");

    switch (calctype){
        case FECalcType::ComputeResidual:
            // For R_d
            localR(1)=viscosity*gpV[1]*test
                       +2*(gpU[1]-1)*Hist*test
                       +(Gc/L)*gpU[1]*test
                       +Gc*L*(gpGradU[1]*grad_test);
            // For R_ux
            localR(2)=Stress.IthRow(1)*grad_test;
            // For R_uy
            localR(3)=Stress.IthRow(2)*grad_test;
            if(nDim==3){
                // For R_uz
                localR(4)=Stress.IthRow(3)*grad_test;
            }
            break;
        case FECalcType::ComputeJacobian:
            // K_d,d  (see Eq. 47)
            localK(1,1)=viscosity*trial*test*ctan[1]
                           +2*trial*Hist*test*ctan[0]
                           +(Gc/L)*trial*test*ctan[0]
                           +Gc*L*(grad_trial*grad_test)*ctan[0];
            //*******************************
            valx=0.0;valy=0.0;valz=0.0;
            for(i=1;i<=3;i++){
                valx+=0.5*(dHdstrain(1,i)+dHdstrain(i,1))*grad_trial(i);
                valy+=0.5*(dHdstrain(2,i)+dHdstrain(i,2))*grad_trial(i);
                valz+=0.5*(dHdstrain(3,i)+dHdstrain(i,3))*grad_trial(i);
            }
            // K_d,ux
            localK(1,2)=2*(gpU[1]-1)*valx*test*ctan[0];
            // K_d,uy
            localK(1,3)=2*(gpU[1]-1)*valy*test*ctan[0];
            if(nDim==3){
                // K_d,uz
                localK(1,4)=2*(gpU[1]-1)*valz*test*ctan[0];
            }

            // K_ux,ux
            localK(2,2)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(1,1,grad_test,grad_trial)*ctan[0];
            // K_ux,uy
            localK(2,3)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(1,2,grad_test,grad_trial)*ctan[0];
            // K_ux,d
            localK(2,1)=dStressdD.IthRow(1)*grad_test*trial*ctan[0];

            // K_uy,uy
            localK(3,3)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(2,2,grad_test,grad_trial)*ctan[0];
            // K_uy,ux
            localK(3,2)=Rank4Materials.at("elasticity_tensor").GetIKjlComponent(2,1,grad_test,grad_trial)*ctan[0];
            // K_uy,d
            localK(3,1)=dStressdD.IthRow(2)*grad_test*trial*ctan[0];
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
                localK(4,1)=dStressdD.IthRow(3)*grad_test*trial*ctan[0];
            }
            break;
        case FECalcType::InitHistoryVariable:
            fill(gpHist.begin(),gpHist.end(),0.0);
            break;
        case FECalcType::UpdateHistoryVariable:
            gpHistOld=gpHist;
            break;
        case FECalcType::Projection:
            gpProj["fx"]=Stress.IthRow(1)*grad_test;// reaction force-x
            gpProj["fy"]=Stress.IthRow(2)*grad_test;// reaction force-y
            gpProj["fz"]=Stress.IthRow(3)*grad_test;// reaction force-z
            break;
        default:
            MessagePrinter::PrintErrorTxt("unsupported FEM calculation type in Miehe fracture element");
            MessagePrinter::AsFem_Exit();
            break;
    }
}