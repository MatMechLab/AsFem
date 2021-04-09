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
//+++ Date   : 2021.04.09
//+++ Purpose: implement the residual and jacobian for general
//+++          stress equilibrium equation
//+++          div(sigma)=0
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/MechanicsElmt.h"

void MechanicsElmt::ComputeAll(const FECalcType &calctype, const int &nDim, const int &nNodes, const int &nDofs,
                               const double &t, const double &dt, const double (&ctan)[2], const Vector3d &gpCoords,
                               const vector<double> &gpU, const vector<double> &gpUold, const vector<double> &gpV,
                               const vector<double> &gpVold, const vector<Vector3d> &gpGradU,
                               const vector<Vector3d> &gpGradUold, const vector<Vector3d> &gpGradV,
                               const vector<Vector3d> &gpGradVold, const double &test, const double &trial,
                               const Vector3d &grad_test, const Vector3d &grad_trial, const Materials &Mate,
                               const Materials &MateOld, map<string, double> &gpProj, MatrixXd &localK,
                               VectorXd &localR) {
    if(calctype==FECalcType::ComputeResidual){
        ComputeResidual(nDim,nNodes,nDofs,t,dt,gpCoords,gpU,gpUold,gpV,gpVold,gpGradU,gpGradUold,gpGradV,gpGradVold,test,grad_test,
                        Mate,MateOld,localR);
    }
    else if(calctype==FECalcType::ComputeJacobian){
        ComputeJacobian(nDim,nNodes,nDofs,t,dt,ctan,gpCoords,gpU,gpUold,gpV,gpVold,gpGradU,gpGradUold,gpGradV,gpGradVold,
                        test,trial,grad_test,grad_trial,Mate,MateOld,localK);
    }
    else if(calctype==FECalcType::Projection){
        ComputeProjection(nDim,nNodes,nDofs,t,dt,ctan,gpCoords,gpU,gpUold,gpV,gpVold,gpGradU,gpGradUold,gpGradV,gpGradVold,
                          test,grad_test,Mate,MateOld,gpProj);
    }
    else{
        MessagePrinter::PrintErrorTxt("unsupported calculation type in MechanicsElmt, please check your related code");
        MessagePrinter::AsFem_Exit();
    }
}
//*****************************************************************************
void MechanicsElmt::ComputeResidual(const int &nDim, const int &nNodes, const int &nDofs, const double &t,
                                    const double &dt, const Vector3d &gpCoords, const vector<double> &gpU,
                                    const vector<double> &gpUold, const vector<double> &gpV,
                                    const vector<double> &gpVold, const vector<Vector3d> &gpGradU,
                                    const vector<Vector3d> &gpGradUold, const vector<Vector3d> &gpGradV,
                                    const vector<Vector3d> &gpGradVold, const double &test, const Vector3d &grad_test,
                                    const Materials &Mate, const Materials &MateOld, VectorXd &localR) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(nDim||nNodes||nDofs||t||dt||gpCoords(1)||gpU[0]||gpUold[0]||gpV[0]||gpVold[0]||
       gpGradU[0](1)||gpGradUold[0](1)||gpGradV[0](1)||gpGradVold[0](1)||
       test||grad_test(1)||
       Mate.ScalarMaterials.size()||MateOld.ScalarMaterials.size()){}

    // calculate the residual contribution of Mechanics problem
    // For R_ux
    localR(1)=Mate.Rank2Materials.at("stress").IthRow(1)*grad_test;
    if(nDim>=2){
        localR(2)=Mate.Rank2Materials.at("stress").IthRow(2)*grad_test;
        if(nDim==3){
            localR(3)=Mate.Rank2Materials.at("stress").IthRow(3)*grad_test;
        }
    }
}
//******************************************************************************
void MechanicsElmt::ComputeJacobian(const int &nDim, const int &nNodes, const int &nDofs, const double &t,
                                    const double &dt, const double (&ctan)[2], const Vector3d &gpCoords,
                                    const vector<double> &gpU, const vector<double> &gpUold, const vector<double> &gpV,
                                    const vector<double> &gpVold, const vector<Vector3d> &gpGradU,
                                    const vector<Vector3d> &gpGradUold, const vector<Vector3d> &gpGradV,
                                    const vector<Vector3d> &gpGradVold, const double &test, const double &trial,
                                    const Vector3d &grad_test, const Vector3d &grad_trial, const Materials &Mate,
                                    const Materials &MateOld, MatrixXd &localK) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(nDim||nNodes||nDofs||t||dt||ctan[0]||gpCoords(1)||gpU[0]||gpUold[0]||gpV[0]||gpVold[0]||
       gpGradU[0](1)||gpGradUold[0](1)||gpGradV[0](1)||gpGradVold[0](1)||
       test||trial||grad_test(1)||grad_trial(1)||
       Mate.VectorMaterials.size()||MateOld.ScalarMaterials.size()){}

    // for the stiffness matrix of mechanics problem
    // K_ux,ux
    localK(1,1)=Mate.Rank4Materials.at("elasticity_tensor").GetIKjlComponent(1,1,grad_test,grad_trial)*ctan[0];
    if(nDim>=2){
        // K_ux,uy
        localK(1,2)=Mate.Rank4Materials.at("elasticity_tensor").GetIKjlComponent(1,2,grad_test,grad_trial)*ctan[0];
        // K_uy,ux
        localK(2,1)=Mate.Rank4Materials.at("elasticity_tensor").GetIKjlComponent(2,1,grad_test,grad_trial)*ctan[0];
        // K_uy,uy
        localK(2,2)=Mate.Rank4Materials.at("elasticity_tensor").GetIKjlComponent(2,2,grad_test,grad_trial)*ctan[0];
        if(nDim==3){
            // K_ux,uz
            localK(1,3)=Mate.Rank4Materials.at("elasticity_tensor").GetIKjlComponent(1,3,grad_test,grad_trial)*ctan[0];
            // K_uy,uz
            localK(2,3)=Mate.Rank4Materials.at("elasticity_tensor").GetIKjlComponent(2,3,grad_test,grad_trial)*ctan[0];
            // K_uz,ux
            localK(3,1)=Mate.Rank4Materials.at("elasticity_tensor").GetIKjlComponent(3,1,grad_test,grad_trial)*ctan[0];
            // K_uz,uy
            localK(3,2)=Mate.Rank4Materials.at("elasticity_tensor").GetIKjlComponent(3,2,grad_test,grad_trial)*ctan[0];
            // K_uz,uz
            localK(3,3)=Mate.Rank4Materials.at("elasticity_tensor").GetIKjlComponent(3,3,grad_test,grad_trial)*ctan[0];
        }
    }
}
//*************************************************
void MechanicsElmt::ComputeProjection(const int &nDim, const int &nNodes, const int &nDofs, const double &t,
                                      const double &dt, const double (&ctan)[2], const Vector3d &gpCoords,
                                      const vector<double> &gpU, const vector<double> &gpUold,
                                      const vector<double> &gpV, const vector<double> &gpVold,
                                      const vector<Vector3d> &gpGradU, const vector<Vector3d> &gpGradUold,
                                      const vector<Vector3d> &gpGradV, const vector<Vector3d> &gpGradVold,
                                      const double &test, const Vector3d &grad_test, const Materials &Mate,
                                      const Materials &MateOld, map<string, double> &gpProj) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(nDim||nNodes||nDofs||t||dt||ctan[0]||gpCoords(1)||gpU[0]||gpUold[0]||gpV[0]||gpVold[0]||
       gpGradU[0](1)||gpGradUold[0](1)||gpGradV[0](1)||gpGradVold[0](1)||
       test||
       Mate.VectorMaterials.size()||MateOld.ScalarMaterials.size()){}

    gpProj["reacforce_x"]=Mate.Rank2Materials.at("stress").IthRow(1)*grad_test;
    gpProj["reacforce_y"]=Mate.Rank2Materials.at("stress").IthRow(2)*grad_test;
    gpProj["reacforce_z"]=Mate.Rank2Materials.at("stress").IthRow(3)*grad_test;
}