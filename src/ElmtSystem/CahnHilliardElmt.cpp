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
//+++          Cahn-Hilliard equation
//+++          dc/dt=div(M*grad(mu))
//+++             mu=d2f/dc-Kappa*Delta(c)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/CahnHilliardElmt.h"

void CahnHilliardElmt::ComputeAll(const FECalcType &calctype, const int &nDim, const int &nNodes, const int &nDofs,
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
        MessagePrinter::PrintErrorTxt("unsupported calculation type in CahnHilliardElmt, please check your related code");
        MessagePrinter::AsFem_Exit();
    }
}
//*******************************************************************
void CahnHilliardElmt::ComputeResidual(const int &nDim, const int &nNodes, const int &nDofs, const double &t,
                                       const double &dt, const Vector3d &gpCoords, const vector<double> &gpU,
                                       const vector<double> &gpUold, const vector<double> &gpV,
                                       const vector<double> &gpVold, const vector<Vector3d> &gpGradU,
                                       const vector<Vector3d> &gpGradUold, const vector<Vector3d> &gpGradV,
                                       const vector<Vector3d> &gpGradVold, const double &test,
                                       const Vector3d &grad_test, const Materials &Mate, const Materials &MateOld,
                                       VectorXd &localR) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(nDim||nNodes||nDofs||t||dt||gpCoords(1)||gpU[0]||gpUold[0]||gpV[0]||gpVold[0]||
       gpGradU[0](1)||gpGradUold[0](1)||gpGradV[0](1)||gpGradVold[0](1)||
       test||grad_test(1)||
       Mate.ScalarMaterials.size()||MateOld.ScalarMaterials.size()){}

    // For R_c
    localR(1)=gpV[1]*test+Mate.ScalarMaterials.at("M")*(gpGradU[2]*grad_test);
    // For R_mu
    localR(2)=gpU[2]*test-Mate.ScalarMaterials.at("dFdc")*test
            -Mate.ScalarMaterials.at("Kappa")*(gpGradU[1]*grad_test);
}
//***************************************************************************
void CahnHilliardElmt::ComputeJacobian(const int &nDim, const int &nNodes, const int &nDofs, const double &t,
                                       const double &dt, const double (&ctan)[2], const Vector3d &gpCoords,
                                       const vector<double> &gpU, const vector<double> &gpUold,
                                       const vector<double> &gpV, const vector<double> &gpVold,
                                       const vector<Vector3d> &gpGradU, const vector<Vector3d> &gpGradUold,
                                       const vector<Vector3d> &gpGradV, const vector<Vector3d> &gpGradVold,
                                       const double &test, const double &trial, const Vector3d &grad_test,
                                       const Vector3d &grad_trial, const Materials &Mate, const Materials &MateOld,
                                       MatrixXd &localK) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(nDim||nNodes||nDofs||t||dt||ctan[0]||gpCoords(1)||gpU[0]||gpUold[0]||gpV[0]||gpVold[0]||
       gpGradU[0](1)||gpGradUold[0](1)||gpGradV[0](1)||gpGradVold[0](1)||
       test||trial||grad_test(1)||grad_trial(1)||
       Mate.VectorMaterials.size()||MateOld.ScalarMaterials.size()){}

    // K_c,c
    localK(1,1)=trial*test*ctan[1]
                +Mate.ScalarMaterials.at("dMdc")*trial*(gpGradU[2]*grad_test)*ctan[0];
    // K_c,mu
    localK(1,2)=Mate.ScalarMaterials.at("M")*grad_trial*grad_test;
    // K_mu,c
    localK(2,1)=-Mate.ScalarMaterials.at("d2Fdc2")*trial*test*ctan[0]
            -Mate.ScalarMaterials.at("Kappa")*grad_trial*grad_test*ctan[0];
    // K_mu,mu
    localK(2,2)=trial*test*ctan[0];
}
//***********************************************************
void CahnHilliardElmt::ComputeProjection(const int &nDim, const int &nNodes, const int &nDofs, const double &t,
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
       test||grad_test(1)||gpProj.size()||
       Mate.VectorMaterials.size()||MateOld.ScalarMaterials.size()){}
}