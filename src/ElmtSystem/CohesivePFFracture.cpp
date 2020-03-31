//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "ElmtSystem/ElmtSystem.h"

void ElmtSystem::CohesivePFFracture(const int &isw,const int &nDim,const int &nNodes,
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
    if(Hist.size()||HistOld.size()){}
    if(VectorMaterials.size()){}
    //*********************************
    //*** related parameters
    //*********************************
    const double viscosity=ScalarMaterials[0];
    const double Gc=ScalarMaterials[1];
    const double L=ScalarMaterials[2];
    const double dG=ScalarMaterials[3];
    const double d2G=ScalarMaterials[4];
    
    //************************************************************
    //*** see the MieheLinearElasticMaterial.cpp
    //*** for the details, one is referred to Miehe's Paper
    //************************************************************
    double H=ScalarMaterials[20];
    RankTwoTensor dHdStrain=Rank2Materials[5];
    

    double valx=0.0,valy=0.0,valz=0.0;

    //***************************
    //*** dofs: ux,uy,d    (2d)
    //***       ux,uy,uz,d (3d)
    //***************************

    if(isw==3||isw==6){
        for(int i=1;i<=nNodes;++i){
            if(nDim==2){
                // R_ux
                rhs(3*i-2)+=Rank2Materials[1].IthRow(1)*shp.shape_grad(i);
                // R_uy
                rhs(3*i-1)+=Rank2Materials[1].IthRow(2)*shp.shape_grad(i);
                // R_d
                rhs(3*i  )+=viscosity*gpV[nDim]*shp.shape_value(i)
                           +(3.0*Gc*L/8.0)*(gpGradU[nDim]*shp.shape_grad(i))
                           +(3.0*Gc/(4.0*L))*gpU[nDim]*shp.shape_value(i)
                           +dG*H*shp.shape_value(i);
            }
            else if(nDim==3){
                // R_ux
                rhs(4*i-3)+=Rank2Materials[1].IthRow(1)*shp.shape_grad(i);
                // R_uy
                rhs(4*i-2)+=Rank2Materials[1].IthRow(2)*shp.shape_grad(i);
                // R_uz
                rhs(4*i-1)+=Rank2Materials[1].IthRow(3)*shp.shape_grad(i);
                // R_d
                rhs(4*i  )+=viscosity*gpV[nDim]*shp.shape_value(i)
                           +(3.0*Gc*L/8.0)*(gpGradU[nDim]*shp.shape_grad(i))
                           +(3.0*Gc/(4.0*L))*gpU[nDim]*shp.shape_value(i)
                           +dG*H*shp.shape_value(i);
            }
            if(isw==6){
                for(int j=1;j<=nNodes;++j){
                    if(nDim==2){
                        // Kux,ux
                        K(3*i-2,3*j-2)+=Rank4Materials[0].GetIKjlComponent(1,1,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kux,uy
                        K(3*i-2,3*j-1)+=Rank4Materials[0].GetIKjlComponent(1,2,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kux,d
                        K(3*i-2,3*j  )+=Rank2Materials[2].IthRow(1)*shp.shape_grad(i)*shp.shape_value(j)*ctan[0];

                        // Kuy,ux
                        K(3*i-1,3*j-2)+=Rank4Materials[0].GetIKjlComponent(2,1,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kuy,uy
                        K(3*i-1,3*j-1)+=Rank4Materials[0].GetIKjlComponent(2,2,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kuy,c
                        K(3*i-1,3*j  )+=Rank2Materials[2].IthRow(2)*shp.shape_grad(i)*shp.shape_value(j)*ctan[0];

                        // Kd,dt
                        K(3*i  ,3*j  )+=viscosity*shp.shape_value(j)*shp.shape_value(i)*ctan[1];
                        // Kd,d
                        K(3*i  ,3*j  )+=(3.0*Gc*L/8.0)*(shp.shape_grad(j)*shp.shape_grad(i))*ctan[0]
                                        +(3.0*Gc/(4.0*L))*shp.shape_value(j)*shp.shape_value(i)*ctan[0]
                                        +d2G*H*shp.shape_value(j)*shp.shape_value(i)*ctan[0];
                        
                        valx=0.0;valy=0.0;
                        for(int k=1;k<=3;++k){
                            valx+=0.5*(dHdStrain(1,k)+dHdStrain(k,1))*shp.shape_grad(j)(k);
                            valy+=0.5*(dHdStrain(2,k)+dHdStrain(k,2))*shp.shape_grad(j)(k);
                        }
                        // Kc,ux
                        K(3*i  ,3*j-2)+=dG*valx*shp.shape_value(i)*ctan[0];
                        // Kc,uy
                        K(3*i  ,3*j-1)+=dG*valy*shp.shape_value(i)*ctan[0];

                    }
                    else if(nDim==3){
                        // Kux,ux
                        K(4*i-3,4*j-3)+=Rank4Materials[0].GetIKjlComponent(1,1,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kux,uy
                        K(4*i-3,4*j-2)+=Rank4Materials[0].GetIKjlComponent(1,2,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kux,uz
                        K(4*i-3,4*j-1)+=Rank4Materials[0].GetIKjlComponent(1,3,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kux,c
                        K(4*i-3,4*j  )+=Rank2Materials[2].IthRow(1)*shp.shape_grad(i)*shp.shape_value(j)*ctan[0];

                        // Kuy,ux
                        K(4*i-2,4*j-3)+=Rank4Materials[0].GetIKjlComponent(2,1,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kuy,uy
                        K(4*i-2,4*j-2)+=Rank4Materials[0].GetIKjlComponent(2,2,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kuy,uz
                        K(4*i-2,4*j-1)+=Rank4Materials[0].GetIKjlComponent(2,3,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kuy,c
                        K(4*i-2,4*j  )+=Rank2Materials[2].IthRow(2)*shp.shape_grad(i)*shp.shape_value(j)*ctan[0];

                        // Kuz,ux
                        K(4*i-1,4*j-3)+=Rank4Materials[0].GetIKjlComponent(3,1,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kuz,uy
                        K(4*i-1,4*j-2)+=Rank4Materials[0].GetIKjlComponent(3,2,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kuz,uy
                        K(4*i-1,4*j-1)+=Rank4Materials[0].GetIKjlComponent(3,3,shp.shape_grad(i),shp.shape_grad(j))*ctan[0];
                        // Kuz,c
                        K(4*i-1,4*j  )+=Rank2Materials[2].IthRow(3)*shp.shape_grad(i)*shp.shape_value(j)*ctan[0];


                        // Kd,dt
                        K(4*i  ,4*j  )+=viscosity*shp.shape_value(j)*shp.shape_value(i)*ctan[1];
                        // Kd,d
                        K(4*i  ,4*j  )+=(3.0*Gc*L/8.0)*(shp.shape_grad(j)*shp.shape_grad(i))*ctan[0]
                                        +(3.0*Gc/(4.0*L))*shp.shape_value(j)*shp.shape_value(i)*ctan[0]
                                        +d2G*H*shp.shape_value(j)*shp.shape_value(i)*ctan[0];
                        
                        valx=0.0;valy=0.0;valz=0.0;
                        for(int k=1;k<=3;++k){
                            valx+=0.5*(dHdStrain(1,k)+dHdStrain(k,1))*shp.shape_grad(j)(k);
                            valy+=0.5*(dHdStrain(2,k)+dHdStrain(k,2))*shp.shape_grad(j)(k);
                            valz+=0.5*(dHdStrain(3,k)+dHdStrain(k,3))*shp.shape_grad(j)(k);
                        }
                        // Kd,ux
                        K(4*i  ,4*j-3)+=dG*valx*shp.shape_value(i)*ctan[0];
                        // Kd,uy
                        K(4*i  ,4*j-2)+=dG*valy*shp.shape_value(i)*ctan[0];
                        // Kd,uz
                        K(4*i  ,4*j-1)+=dG*valz*shp.shape_value(i)*ctan[0];
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
        // TODO: I need a better design for history update
        //       in material, I already calculate the current hist
        //       it seems quite stupid to do it here again?
        // update hist
        //Hist=HistOld;
    }
    else if(isw==9){
        // do projection
        Proj[0]=ScalarMaterials[3];// vonMises stress
        Proj[1]=ScalarMaterials[4];// hydrostatic stress
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