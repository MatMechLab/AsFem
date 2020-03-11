//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "ElmtSystem/ElmtSystem.h"

void ElmtSystem::TensorCahnHilliard(const int &isw,const int &nDim,const int &nNodes,
            const double &t,const double &dt,const double (&ctan)[2],
            const Vector3d &gpCoord,
            const vector<double> &gpU,const vector<double> &gpV,
            const vector<Vector3d> &gpGradU,const vector<Vector3d> &gpGradV,
            const ShapeFun &shp,
            const vector<double> &ScalarMaterials,
            const vector<Vector3d> &VectorMaterials,
            const vector<RankTwoTensor> &Rank2Materials,
            const vector<RankFourTensor> &Rank4Materials,
            vector<double> &Hist,const vector<double> &HistOld,vector<double> &Proj,
            MatrixXd &K,VectorXd &rhs){
    //********************************
    //*** to get rid of warnings
    //********************************
    if(t||dt){}
    if(gpCoord(1)||gpU[0]||gpV[0]||gpGradU[0](1)||gpGradV[0](1)){}
    if(VectorMaterials.size()){}
    if(Rank2Materials.size()||Rank4Materials.size()){}
    //********************************

    //**************************************
    //*** related parameters
    //**************************************
    RankTwoTensor M=Rank2Materials[0];
    RankTwoTensor dMdc=Rank2Materials[1];
    RankTwoTensor Kappa=Rank2Materials[2];
    double F=ScalarMaterials[1];
    double dFdc=ScalarMaterials[2];
    double d2Fdc2=ScalarMaterials[3];
    if(isw==3||isw==6){
        for(int i=1;i<=nNodes;++i){
            //Rc
            rhs(2*i-1)+=gpV[0]*shp.shape_value(i)
                       +M*gpGradU[1]*shp.shape_grad(i);
            //Rmu
            rhs(2*i  )+=gpU[1]*shp.shape_value(i)
                       -dFdc*shp.shape_value(i)
                       -Kappa*gpGradU[0]*shp.shape_grad(i);
            if(isw==6){
                for(int j=1;j<=nNodes;++j){
                    // Kc,cdot
                    K(2*i-1,2*j-1)+=shp.shape_value(j)*shp.shape_value(i)*ctan[1];
                    // Kc,cc
                    K(2*i-1,2*j-1)+=dMdc*shp.shape_value(j)*gpGradU[1]*shp.shape_grad(i)*ctan[0];
                    // Kc,mu
                    K(2*i-1,2*j  )+=M*shp.shape_grad(j)*shp.shape_grad(i)*ctan[0];

                    // Kmu,c
                    K(2*i  ,2*j-1)+=-d2Fdc2*shp.shape_value(j)*shp.shape_value(i)*ctan[0]
                                   -Kappa*shp.shape_grad(j)*shp.shape_grad(i)*ctan[0];
                    // Kmu,mu
                    K(2*i  ,2*j  )+=shp.shape_value(j)*shp.shape_value(i)*ctan[0];
                }
            }
        }
    }
    else if(isw==4){
        // init history variable
        fill(Hist.begin(),Hist.end(),0.0);
    }
    else if(isw==8){
        // update hist
        Hist=HistOld;
    }
    else if(isw==9){
        Proj[0]=F+0.5*Kappa.Trace()*(gpGradU[0]*gpGradU[0]);// toal free energy
        Proj[1]=gpU[1]; // final chemical potential
        if(nDim==2){
            Proj[2]=gpGradU[0](1);
            Proj[3]=gpGradU[0](2);
            Proj[4]=gpGradU[1](1);
            Proj[5]=gpGradU[1](2);
        }
        else if(nDim==3){
            Proj[2]=gpGradU[0](1);
            Proj[3]=gpGradU[0](2);
            Proj[4]=gpGradU[0](3);

            Proj[5]=gpGradU[1](1);
            Proj[6]=gpGradU[1](2);
            Proj[7]=gpGradU[1](3);
        }
    }
}