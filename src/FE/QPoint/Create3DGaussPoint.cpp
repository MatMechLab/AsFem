//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "FE/QPoint.h"

void QPoint::Create3DGaussPoint(MeshType meshtype){
    switch (meshtype){
        case MeshType::HEX8:
        case MeshType::HEX20:
        case MeshType::HEX27:
        {
            QPoint qp1d(1,GetQpOrder());
            qp1d.Create1DGaussPoint();
            _nQpPoints=qp1d.GetQpPointsNum()*qp1d.GetQpPointsNum()*qp1d.GetQpPointsNum();
            _qp_coords.resize(_nQpPoints*(GetDim()+1),0.0);
            int l=0;
            // for this regular type mesh , it is just simple tensor product
            for(int k=0;k<qp1d.GetQpPointsNum();k++){
                for(int j=0;j<qp1d.GetQpPointsNum();j++){
                    for(int i=0;i<qp1d.GetQpPointsNum();i++){
                        l+=1;
                        (*this)(l,1)=qp1d(i+1,1);
                        (*this)(l,2)=qp1d(j+1,1);
                        (*this)(l,3)=qp1d(k+1,1);
                        (*this)(l,0)=qp1d(i+1,0)*qp1d(j+1,0)*qp1d(k+1,0);
                        // cout<<"gpInd="<<l
                        //     <<": xi="<<(*this)(l,1)
                        //     <<", eta="<<(*this)(l,2)
                        //     <<": zeta="<<(*this)(l,3)<<endl;
                    }
                }
            }
            return;
            break;
        }
        case MeshType::TET4:
        case MeshType::TET10:
        {
            switch (GetQpOrder())
            {
                case 0:
                case 1:
                {
                    _nQpPoints=1;
                    _qp_coords.resize(_nQpPoints*(GetDim()+1),0.0);
                    (*this)(1,1)=0.25;
                    (*this)(1,2)=0.25;
                    (*this)(1,3)=0.25;
                    (*this)(1,0)=1.0/6.0;
                    return;
                }
                case 2:
                {
                    _nQpPoints=4;
                    _qp_coords.resize(_nQpPoints*(GetDim()+1),0.0);
                    const double a=0.585410196624969;
                    const double b=0.138196601125011;

                    (*this)(1,1)=a;
                    (*this)(1,2)=b;
                    (*this)(1,3)=b;
                    (*this)(1,0)=1.0/24.0;

                    (*this)(2,1)=b;
                    (*this)(2,2)=a;
                    (*this)(2,3)=b;
                    (*this)(2,0)=1.0/24.0;
                
                    (*this)(3,1)=b;
                    (*this)(3,2)=b;
                    (*this)(3,3)=a;
                    (*this)(3,0)=1.0/24.0;

                    (*this)(4,1)=b;
                    (*this)(4,2)=b;
                    (*this)(4,3)=b;
                    (*this)(4,0)=1.0/24.0;
                    return;
                }
                case 3:
                {
                    _nQpPoints=5;
                    _qp_coords.resize(_nQpPoints*(GetDim()+1),0.0);

                    (*this)(1,1)=0.25;
                    (*this)(1,2)=0.25;
                    (*this)(1,3)=0.25;
                    (*this)(1,0)=-2.0/15.0;

                    (*this)(2,1)=0.5;
                    (*this)(2,2)=1.0/6.0;
                    (*this)(2,3)=1.0/6.0;
                    (*this)(2,0)=0.075;
                
                    (*this)(3,1)=1.0/6.0;
                    (*this)(3,2)=0.5;
                    (*this)(3,3)=1.0/6.0;
                    (*this)(3,0)=0.075;

                    (*this)(4,1)=1.0/6.0;
                    (*this)(4,2)=1.0/6.0;
                    (*this)(4,3)=0.5;
                    (*this)(4,0)=0.075;

                    (*this)(5,1)=1.0/6.0;
                    (*this)(5,2)=1.0/6.0;
                    (*this)(5,3)=1.0/6.0;
                    (*this)(5,0)=0.075;

                    return;
                }
                default:
                {
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: invalid gauss integration order(=%3d) for 3D case    !!!   ***\n",GetQpOrder());
                    Msg_AsFem_Exit();
                    break;
                }
            }
            return;
        }
        default:
        {
            PetscPrintf(PETSC_COMM_WORLD,"*** Error: unsupported mesh type for 3D gauss integration       !!!   ***\n");
            Msg_AsFem_Exit();
            break;
        }
    }
}

//******************************************************
void QPoint::Create3DGaussLobattoPoint(MeshType meshtype){
    switch (meshtype){
        case MeshType::HEX8:
        case MeshType::HEX20:
        case MeshType::HEX27:
        {
            QPoint qp1d(1,GetQpOrder());
            qp1d.SetQPointType("gausslobatto");
            qp1d.Create1DGaussPoint();
            _nQpPoints=qp1d.GetQpPointsNum()*qp1d.GetQpPointsNum()*qp1d.GetQpPointsNum();
            _qp_coords.resize(_nQpPoints*(GetDim()+1),0.0);
            int l=0;
            // for this regular type mesh , it is just simple tensor product
            for(int i=0;i<qp1d.GetQpPointsNum();i++){
                for(int j=0;j<qp1d.GetQpPointsNum();j++){
                    for(int k=0;k<qp1d.GetQpPointsNum();k++){
                        l+=1;
                        (*this)(l,1)=qp1d(i+1,1);
                        (*this)(l,2)=qp1d(j+1,1);
                        (*this)(l,3)=qp1d(k+1,1);
                        (*this)(l,0)=qp1d(i+1,0)*qp1d(j+1,0)*qp1d(k+1,0);
                    }
                }
            }
            return;
        }
        default:
            PetscPrintf(PETSC_COMM_WORLD,"*** Error: unsupported mesh type for 3D gauss lobatto integration !!! ***\n");
            Msg_AsFem_Exit();
            break;
    }
}