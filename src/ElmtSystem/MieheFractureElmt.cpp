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
//+++ Purpose: implement the residual and jacobian for Miehe's
//+++          phase field fracture model
//+++          1) eta*dD/dt=2(1-D)H-(Gc/L)D+Gc*L*Lap(D) (see Eq. 47)
//+++          2) div(\Sigma)=0
//+++ Reference: A phase field model for rate-independent crack propagation:
//+++            Robust algorithmic implementation based on operator splits
//+++ DOI    : https://doi.org/10.1016/j.cma.2010.04.011
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ElmtSystem/MieheFractureElmt.h"

void MieheFractureElmt::ComputeAll(const FECalcType &calctype, const int &nDim, const int &nNodes, const int &nDofs,
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
        MessagePrinter::PrintErrorTxt("unsupported calculation type in MieheFractureElmt, please check your related code");
        MessagePrinter::AsFem_Exit();
    }
}
//**************************************************
void MieheFractureElmt::ComputeResidual(const int &nDim, const int &nNodes, const int &nDofs, const double &t,
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

    // For R_d
    double viscosity=Mate.ScalarMaterials.at("viscosity");
    double Gc=Mate.ScalarMaterials.at("Gc");
    double L=Mate.ScalarMaterials.at("L");
    double Hist=Mate.ScalarMaterials.at("Hist");
    RankTwoTensor Stress=Mate.Rank2Materials.at("stress");

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
}
//*************************************************************
void MieheFractureElmt::ComputeJacobian(const int &nDim, const int &nNodes, const int &nDofs, const double &t,
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

    //************************************************************
    //*** some intermediate variables
    //************************************************************
    int i;
    double valx,valy,valz;
    double viscosity=Mate.ScalarMaterials.at("viscosity");
    double Gc=Mate.ScalarMaterials.at("Gc");
    double L=Mate.ScalarMaterials.at("L");
    double Hist=Mate.ScalarMaterials.at("Hist");
    RankTwoTensor dStressdD=Mate.Rank2Materials.at("dstressdD");
    RankTwoTensor dHdstrain=Mate.Rank2Materials.at("dHdstrain");

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
    localK(2,2)=Mate.Rank4Materials.at("jacobian").GetIKjlComponent(1,1,grad_test,grad_trial)*ctan[0];
    // K_ux,uy
    localK(2,3)=Mate.Rank4Materials.at("jacobian").GetIKjlComponent(1,2,grad_test,grad_trial)*ctan[0];
    // K_ux,d
    localK(2,1)=dStressdD.IthRow(1)*grad_test*trial*ctan[0];

    // K_uy,uy
    localK(3,3)=Mate.Rank4Materials.at("jacobian").GetIKjlComponent(2,2,grad_test,grad_trial)*ctan[0];
    // K_uy,ux
    localK(3,2)=Mate.Rank4Materials.at("jacobian").GetIKjlComponent(2,1,grad_test,grad_trial)*ctan[0];
    // K_uy,d
    localK(3,1)=dStressdD.IthRow(2)*grad_test*trial*ctan[0];
    if(nDim==3){
        // K_ux,uz
        localK(2,4)=Mate.Rank4Materials.at("jacobian").GetIKjlComponent(1,3,grad_test,grad_trial)*ctan[0];
        // K_uy,uz
        localK(3,4)=Mate.Rank4Materials.at("jacobian").GetIKjlComponent(2,4,grad_test,grad_trial)*ctan[0];

        // K_uz,uz
        localK(4,4)=Mate.Rank4Materials.at("jacobian").GetIKjlComponent(3,3,grad_test,grad_trial)*ctan[0];
        // K_uz,ux
        localK(4,2)=Mate.Rank4Materials.at("jacobian").GetIKjlComponent(3,1,grad_test,grad_trial)*ctan[0];
        // K_uz,uy
        localK(4,3)=Mate.Rank4Materials.at("jacobian").GetIKjlComponent(3,2,grad_test,grad_trial)*ctan[0];
        // K_uz,d
        localK(4,1)=dStressdD.IthRow(3)*grad_test*trial*ctan[0];
    }
}
//**************************************************************************
void MieheFractureElmt::ComputeProjection(const int &nDim, const int &nNodes, const int &nDofs, const double &t,
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