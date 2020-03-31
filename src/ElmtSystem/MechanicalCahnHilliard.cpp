//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "ElmtSystem/ElmtSystem.h"

void ElmtSystem::MechanicalCahnHilliard(const int &isw,const int &nDim,const int &nNodes,
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
    if(VectorMaterials.size()) {}
    if(Rank2Materials.size()||Rank4Materials.size()){}
    //********************************

    //**************************************
    //*** related parameters
    //**************************************
    double Kappa=ScalarMaterials[0];
    double M=ScalarMaterials[1];
    double dMdc=ScalarMaterials[2];
    double F=ScalarMaterials[5];
    double dFdc=ScalarMaterials[6];
    double d2Fdc2=ScalarMaterials[7];
    if(isw==3||isw==6){
        for(int i=1;i<=nNodes;++i){
            if(nDim==2){
                // R_ux
                rhs(4*i-3)+=Rank2Materials[1].IthRow(1)*shp.shape_grad(i);
                // R_uy
                rhs(4*i-2)+=Rank2Materials[1].IthRow(2)*shp.shape_grad(i);
                //Rc
                rhs(4*i-1)+=gpV[2]*shp.shape_value(i)
                           +M*(gpGradU[3]*shp.shape_grad(i));
                //Rmu
                rhs(4*i  )+=gpU[3]*shp.shape_value(i)
                           -dFdc*shp.shape_value(i)
                           -Kappa*(gpGradU[2]*shp.shape_grad(i));
            }
            else if(nDim==3){
                // R_ux
                rhs(5*i-4)+=Rank2Materials[1].IthRow(1)*shp.shape_grad(i);
                // R_uy
                rhs(5*i-3)+=Rank2Materials[1].IthRow(2)*shp.shape_grad(i);
                // R_uz
                rhs(5*i-2)+=Rank2Materials[1].IthRow(3)*shp.shape_grad(i);
                //Rc
                rhs(5*i-1)+=gpV[3]*shp.shape_value(i)
                           +M*(gpGradU[4]*shp.shape_grad(i));
                //Rmu
                rhs(5*i  )+=gpU[4]*shp.shape_value(i)
                           -dFdc*shp.shape_value(i)
                           -Kappa*(gpGradU[3]*shp.shape_grad(i));
            }
            if(isw==6){
                for(int j=1;j<=nNodes;++j){
                    if(nDim==2){
                        // Kux,ux
                        K(4*i-3,4*j-3)+=Rank4Materials[0].GetIKjlComponent(1,1,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kux,uy
                        K(4*i-3,4*j-2)+=Rank4Materials[0].GetIKjlComponent(1,2,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kux,c
                        K(4*i-3,4*j-1)+=Rank2Materials[2].IthRow(1)*shp.shape_grad(i)*shp.shape_value(j)*ctan[0];

                        // Kuy,ux
                        K(4*i-2,4*j-3)+=Rank4Materials[0].GetIKjlComponent(2,1,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kuy,uy
                        K(4*i-2,4*j-2)+=Rank4Materials[0].GetIKjlComponent(2,2,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kuy,c
                        K(4*i-2,4*j-1)+=Rank2Materials[2].IthRow(2)*shp.shape_grad(i)*shp.shape_value(j)*ctan[0];

                        // Kc,cdot
                        K(4*i-1,4*j-1)+=shp.shape_value(j)*shp.shape_value(i)*ctan[1];
                        // Kc,cc
                        K(4*i-1,4*j-1)+=dMdc*shp.shape_value(j)*(gpGradU[3]*shp.shape_grad(i))*ctan[0];
                        // Kc,mu
                        K(4*i-1,4*j  )+=M*(shp.shape_grad(j)*shp.shape_grad(i))*ctan[0];

                        // Kmu,c
                        K(4*i  ,4*j-1)+=-d2Fdc2*shp.shape_value(j)*shp.shape_value(i)*ctan[0]
                                       -Kappa*(shp.shape_grad(j)*shp.shape_grad(i))*ctan[0];
                        // Kmu,mu
                        K(4*i  ,4*j  )+=shp.shape_value(j)*shp.shape_value(i)*ctan[0];
                    }
                    else if(nDim==3){
                        // Kux,ux
                        K(5*i-4,5*j-4)+=Rank4Materials[0].GetIKjlComponent(1,1,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kux,uy
                        K(5*i-4,5*j-3)+=Rank4Materials[0].GetIKjlComponent(1,2,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kux,uz
                        K(5*i-4,5*j-2)+=Rank4Materials[0].GetIKjlComponent(1,3,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kux,c
                        K(5*i-4,5*j-1)+=Rank2Materials[2].IthRow(1)*shp.shape_grad(i)*shp.shape_value(j)*ctan[0];

                        // Kuy,ux
                        K(5*i-3,5*j-4)+=Rank4Materials[0].GetIKjlComponent(2,1,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kuy,uy
                        K(5*i-3,5*j-3)+=Rank4Materials[0].GetIKjlComponent(2,2,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kuy,uz
                        K(5*i-3,5*j-2)+=Rank4Materials[0].GetIKjlComponent(2,3,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kuy,c
                        K(5*i-3,5*j-1)+=Rank2Materials[2].IthRow(2)*shp.shape_grad(i)*shp.shape_value(j)*ctan[0];

                        // Kuz,ux
                        K(5*i-2,5*j-4)+=Rank4Materials[0].GetIKjlComponent(3,1,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kuz,uy
                        K(5*i-2,5*j-3)+=Rank4Materials[0].GetIKjlComponent(3,2,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kuz,uy
                        K(5*i-2,5*j-2)+=Rank4Materials[0].GetIKjlComponent(3,3,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kuz,c
                        K(5*i-2,5*j-1)+=Rank2Materials[2].IthRow(3)*shp.shape_grad(i)*shp.shape_value(j)*ctan[0];

                        // Kc,cdot
                        K(5*i-1,5*j-1)+=shp.shape_value(j)*shp.shape_value(i)*ctan[1];
                        // Kc,cc
                        K(5*i-1,5*j-1)+=dMdc*shp.shape_value(j)*(gpGradU[4]*shp.shape_grad(i))*ctan[0];
                        // Kc,mu
                        K(5*i-1,5*j  )+=M*(shp.shape_grad(j)*shp.shape_grad(i))*ctan[0];

                        // Kmu,c
                        K(5*i  ,5*j-1)+=-d2Fdc2*shp.shape_value(j)*shp.shape_value(i)*ctan[0]
                                       -Kappa*(shp.shape_grad(j)*shp.shape_grad(i))*ctan[0];
                        // Kmu,mu
                        K(5*i  ,5*j  )+=shp.shape_value(j)*shp.shape_value(i)*ctan[0];
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
        Proj[0]=F+0.5*Kappa*(gpGradU[0]*gpGradU[0]);// toal free energy
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