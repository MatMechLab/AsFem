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
//+++ Date   : 2021.04.10
//+++ Purpose: implement the residual and jacobian for Poisson equation
//+++          Sigma*div(grad(phi))=F
//+++          Both Sigma and F could be nonlinear function !
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/PoissonElmt.h"

void PoissonElmt::ComputeAll(const FECalcType &calctype, const int &nDim, const int &nNodes, const int &nDofs,
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
        MessagePrinter::PrintErrorTxt("unsupported calculation type in PoissonElmt, please check your related code");
        MessagePrinter::AsFem_Exit();
    }
}
//***************************************************************************
void PoissonElmt::ComputeResidual(const int &nDim, const int &nNodes, const int &nDofs, const double &t,
                                  const double &dt, const Vector3d &gpCoords, const vector<double> &gpU,
                                  const vector<double> &gpUold, const vector<double> &gpV, const vector<double> &gpVold,
                                  const vector<Vector3d> &gpGradU, const vector<Vector3d> &gpGradUold,
                                  const vector<Vector3d> &gpGradV, const vector<Vector3d> &gpGradVold,
                                  const double &test, const Vector3d &grad_test, const Materials &Mate,
                                  const Materials &MateOld, VectorXd &localR) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(nDim||nNodes||nDofs||t||dt||gpCoords(1)||gpU[0]||gpUold[0]||gpV[0]||gpVold[0]||
       gpGradU[0](1)||gpGradUold[0](1)||gpGradV[0](1)||gpGradVold[0](1)||
       test||grad_test(1)||
       Mate.ScalarMaterials.size()||MateOld.ScalarMaterials.size()){}

    localR(1)=Mate.ScalarMaterials.at("sigma")*(gpGradU[1]*grad_test)+Mate.ScalarMaterials.at("f")*test;
}
//*****************************************************************************
void PoissonElmt::ComputeJacobian(const int &nDim, const int &nNodes, const int &nDofs, const double &t,
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

    localK(1,1)=Mate.ScalarMaterials.at("dsigmadu")*trial*(gpGradU[1]*grad_test)*ctan[0]
            +Mate.ScalarMaterials.at("sigma")*grad_trial*grad_test*ctan[0]
            +Mate.ScalarMaterials.at("dfdu")*trial*test*ctan[0];
}
//*******************************************************************************
void PoissonElmt::ComputeProjection(const int &nDim, const int &nNodes, const int &nDofs, const double &t,
                                    const double &dt, const double (&ctan)[2], const Vector3d &gpCoords,
                                    const vector<double> &gpU, const vector<double> &gpUold, const vector<double> &gpV,
                                    const vector<double> &gpVold, const vector<Vector3d> &gpGradU,
                                    const vector<Vector3d> &gpGradUold, const vector<Vector3d> &gpGradV,
                                    const vector<Vector3d> &gpGradVold, const double &test, const Vector3d &grad_test,
                                    const Materials &Mate, const Materials &MateOld, map<string, double> &gpProj) {
    //***********************************************************
    //*** get rid of unused warning
    //***********************************************************
    if(nDim||nNodes||nDofs||t||dt||ctan[0]||gpCoords(1)||gpU[0]||gpUold[0]||gpV[0]||gpVold[0]||
       gpGradU[0](1)||gpGradUold[0](1)||gpGradV[0](1)||gpGradVold[0](1)||
       test||grad_test(1)||gpProj.size()||
       Mate.VectorMaterials.size()||MateOld.ScalarMaterials.size()){}

}