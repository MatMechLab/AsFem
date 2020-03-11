//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "ElmtSystem/ElmtSystem.h"

void ElmtSystem::Thermal(const int &isw,const int &nDim,const int &nNodes,
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
    if(nDim||t||dt){};
    if(VectorMaterials.size()) {}
    if(Rank2Materials.size()){}
    if(Rank4Materials.size()){}
    
    // In constpoisson material:
    // ScalarMaterials[0]=alpha;     // alpha--> thermal diffusivity=lambda/(rho*c)
    // ScalarMaterials[1]=dalpha/dT;// alpha's derivative
    // ScalarMaterials[2]=F;        // F---> source term
    // ScalarMaterials[2]=dF/dT     //
    double Alpha=ScalarMaterials[0];
    double dAlphadT=ScalarMaterials[1];
    double F=ScalarMaterials[2];
    double dFdT=ScalarMaterials[3];
    if(isw==3||isw==6){
        for(int i=1;i<=nNodes;++i){
            rhs(i)+=gpV[0]*shp.shape_value(i)
                   +Alpha*(gpGradU[0]*shp.shape_grad(i))
                   -F*shp.shape_value(i);
            if(isw==6){
                for(int j=1;j<=nNodes;++j){
                    K(i,j)+=shp.shape_value(j)*shp.shape_value(i)*ctan[1]
                           +dAlphadT*shp.shape_value(j)*(gpGradU[0]*shp.shape_grad(i))*ctan[0]
                           +Alpha*(shp.shape_grad(j)*shp.shape_grad(i))*ctan[0]
                           -dFdT*shp.shape_value(j)*shp.shape_value(i)*ctan[0];
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
        // do projection
        Proj[0]=Alpha*gpGradU[0](1);
        Proj[1]=Alpha*gpGradU[0](2);
        Proj[2]=Alpha*gpGradU[0](3);
        Proj[3]=F;
        Proj[4]=gpU[0];
        Proj[5]=gpCoord(1);
        Proj[6]=gpV[0];
        Proj[7]=gpGradV[0](1);
    }
}