//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "ElmtSystem/ElmtSystem.h"

void ElmtSystem::Mechanics(const int &isw,const int &nDim,const int &nNodes,
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
    //********************************
    if(isw==3||isw==6){
        for(int i=1;i<=nNodes;++i){
            if(nDim==2){
                rhs(2*i-1)+=Rank2Materials[1].IthRow(1)*shp.shape_grad(i);
                rhs(2*i  )+=Rank2Materials[1].IthRow(2)*shp.shape_grad(i);
            }
            else if(nDim==3){
                rhs(3*i-2)+=Rank2Materials[1].IthRow(1)*shp.shape_grad(i);
                rhs(3*i-1)+=Rank2Materials[1].IthRow(2)*shp.shape_grad(i);
                rhs(3*i  )+=Rank2Materials[1].IthRow(3)*shp.shape_grad(i);
            }
            if(isw==6){
                for(int j=1;j<=nNodes;++j){
                    if(nDim==2){
                        // Kux,ux
                        K(2*i-1,2*j-1)+=Rank4Materials[0].GetIKjlComponent(1,1,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kux,uy
                        K(2*i-1,2*j  )+=Rank4Materials[0].GetIKjlComponent(1,2,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        
                        // Kuy,ux
                        K(2*i  ,2*j-1)+=Rank4Materials[0].GetIKjlComponent(2,1,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kuy,uy
                        K(2*i  ,2*j  )+=Rank4Materials[0].GetIKjlComponent(2,2,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];

                    }
                    else if(nDim==3){
                        // Kux,ux
                        K(3*i-2,3*j-2)+=Rank4Materials[0].GetIKjlComponent(1,1,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kux,uy
                        K(3*i-2,3*j-1)+=Rank4Materials[0].GetIKjlComponent(1,2,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kux,uz
                        K(3*i-2,3*j  )+=Rank4Materials[0].GetIKjlComponent(1,3,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];

                        // Kuy,ux
                        K(3*i-1,3*j-2)+=Rank4Materials[0].GetIKjlComponent(2,1,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kuy,uy
                        K(3*i-1,3*j-1)+=Rank4Materials[0].GetIKjlComponent(2,2,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kuy,uz
                        K(3*i-1,3*j  )+=Rank4Materials[0].GetIKjlComponent(2,3,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];

                        // Kuz,ux
                        K(3*i  ,3*j-2)+=Rank4Materials[0].GetIKjlComponent(3,1,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kuz,uy
                        K(3*i  ,3*j-1)+=Rank4Materials[0].GetIKjlComponent(3,2,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kuz,uy
                        K(3*i  ,3*j  )+=Rank4Materials[0].GetIKjlComponent(3,3,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                    }
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
        Proj[0]=ScalarMaterials[0];// vonMises stress
        Proj[1]=ScalarMaterials[1];// hydrostatic stress
        // to get rid of warnings
        if(nDim==2){
            Proj[2]=Rank2Materials[1](1,1);//sigxx
            Proj[3]=Rank2Materials[1](2,2);//sigxy
            Proj[4]=Rank2Materials[1](1,2);//sigyy

            Proj[5]=Rank2Materials[0](1,1);//strxx
            Proj[6]=Rank2Materials[0](2,2);//strxy
            Proj[7]=Rank2Materials[0](1,2);//stryy
        }
        else if(nDim==3){
            Proj[2]=Rank2Materials[1](1,1);//sigxx
            Proj[3]=Rank2Materials[1](2,2);//sigyy
            Proj[4]=Rank2Materials[1](3,3);//sigzz
            Proj[5]=Rank2Materials[1](2,3);//sigyz
            Proj[6]=Rank2Materials[1](1,3);//sigxz
            Proj[7]=Rank2Materials[1](1,2);//sigxy

            Proj[ 8]=Rank2Materials[0](1,1);//strxx
            Proj[ 9]=Rank2Materials[0](2,2);//stryy
            Proj[10]=Rank2Materials[0](3,3);//strzz
            Proj[11]=Rank2Materials[0](2,3);//stryz
            Proj[12]=Rank2Materials[0](1,3);//strxz
            Proj[13]=Rank2Materials[0](1,2);//strxy
        }
    }
}