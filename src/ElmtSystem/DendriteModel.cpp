//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "ElmtSystem/ElmtSystem.h"

void ElmtSystem::DendriteModel(const int &isw,const int &nDim,const int &nNodes,
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
    if(Rank2Materials.size()){}
    if(Rank4Materials.size()){}
    //*********************************
    //*** related parameters
    //*********************************
    double L=ScalarMaterials[0];
    double Conduct=ScalarMaterials[1];
    double eta=ScalarMaterials[2];
    double dFdphi=ScalarMaterials[4];
    double d2Fdphi=ScalarMaterials[5];
    double d2FdphidT=ScalarMaterials[6];
    double k=ScalarMaterials[7];
    double dkdtheta=ScalarMaterials[8];

    Vector3d dkdgradphi =VectorMaterials[0];
    Vector3d ddkdgradphi=VectorMaterials[1];
    Vector3d V          =VectorMaterials[2];
    Vector3d dVdphi(0.0);
    

    //***************************
    //*** dofs: ux,uy,d    (2d)
    //***       ux,uy,uz,d (3d)
    //***************************

    if(isw==3||isw==6){
        for(int i=1;i<=nNodes;++i){
            if(nDim==2){
                // R_phi
                rhs(2*i-1)+=gpV[0]*shp.shape_value(i)
                           +L*k*dkdtheta*(V*shp.shape_grad(i))
                           +L*k*k*(gpGradU[0]*shp.shape_grad(i))
                           +L*dFdphi*shp.shape_value(i);
                // cout<<"i="<<i<<":"<<L*k*dkdtheta*(V*shp.shape_grad(i))<<endl;
                // R_T
                rhs(2*i  )+=gpV[1]*shp.shape_value(i)
                           +Conduct*(gpGradU[1]*shp.shape_grad(i))
                           -eta*gpV[0]*shp.shape_value(i);
            }
            else if(nDim==3){
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: current dendrite model only support 2d case          !!!   ***\n");
                Msg_AsFem_Exit();
            }
            if(isw==6){
                for(int j=1;j<=nNodes;++j){
                    // Kphi,phidot
                    K(2*i-1,2*j-1)+=shp.shape_value(j)*shp.shape_value(i)*ctan[1];
                    // Kphi,phi
                    dVdphi(1)=-shp.shape_grad(j)(2);
                    dVdphi(2)= shp.shape_grad(j)(1);
                    dVdphi(3)= 0.0;
                    K(2*i-1,2*j-1)+=L*(dkdgradphi*shp.shape_grad(j))*dkdtheta*(V*shp.shape_grad(i))*ctan[0]
                                   +L*k*(ddkdgradphi*shp.shape_grad(j))*(V*shp.shape_grad(i))*ctan[0]
                                   +L*k*dkdtheta*(dVdphi*shp.shape_grad(i))*ctan[0]
                                   +L*2*k*(dkdgradphi*shp.shape_grad(j))*(gpGradU[0]*shp.shape_grad(i))*ctan[0]
                                   +L*k*k*(shp.shape_grad(j)*shp.shape_grad(i))*ctan[0]
                                   +L*d2Fdphi*shp.shape_value(j)*shp.shape_value(i)*ctan[0];

                    // Kphi,t
                    K(2*i-1,2*j  )+=L*d2FdphidT*shp.shape_value(j)*shp.shape_value(i)*ctan[0];
                    
                    // Kt,tdot
                    K(2*i  ,2*j  )+=shp.shape_value(j)*shp.shape_value(i)*ctan[1];

                    // Kt,t
                    K(2*i  ,2*j  )+=Conduct*(shp.shape_grad(j)*shp.shape_grad(i))*ctan[0];
                    // Kt,phi
                    K(2*i  ,2*j-1)+=-eta*shp.shape_value(j)*shp.shape_value(i)*ctan[1];
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
        Proj[0]=ScalarMaterials[3];// F
        Proj[1]=dFdphi;// dFdphi
        Proj[2]=ScalarMaterials[7];// K
        Proj[3]=ScalarMaterials[8];// dK/dtheta
        Proj[4]=gpGradU[1](1);
        Proj[5]=gpGradU[1](2);
    }
}