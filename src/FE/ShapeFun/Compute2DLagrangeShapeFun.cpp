//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "FE/ShapeFun.h"

void ShapeFun::Compute2DLagrangeShapeFun(const double &xi,const double &eta,const Nodes &nodes,bool flag){
    _DetJac=0.0;
    switch (GetMeshType()){
        case MeshType::QUAD4:
        {
            (*this)(1,0)=(1.0-xi)*(1.0-eta)/4.0;
            (*this)(2,0)=(1.0+xi)*(1.0-eta)/4.0;
            (*this)(3,0)=(1.0+xi)*(1.0+eta)/4.0;
            (*this)(4,0)=(1.0-xi)*(1.0+eta)/4.0;

            (*this)(1,1)= (eta-1.0)/4.0;
            (*this)(1,2)= (xi -1.0)/4.0;

            (*this)(2,1)= (1.0-eta)/4.0;
            (*this)(2,2)=-(1.0+xi )/4.0;

            (*this)(3,1)= (1.0+eta)/4.0;
            (*this)(3,2)= (1.0+xi )/4.0;

            (*this)(4,1)=-(1.0+eta)/4.0;
            (*this)(4,2)= (1.0-xi )/4.0;
            break;
        }
        case MeshType::QUAD8:
        {
            // 2D-8Nodes rectangle element
            (*this)(1,0)=(1.0-xi)*(1.0-eta)*(-xi-eta-1.0)/4.0;
            (*this)(2,0)=(1.0+xi)*(1.0-eta)*( xi-eta-1.0)/4.0;
            (*this)(3,0)=(1.0+xi)*(1.0+eta)*( xi+eta-1.0)/4.0;
            (*this)(4,0)=(1.0-xi)*(1.0+eta)*(-xi+eta-1.0)/4.0;
            (*this)(5,0)=(1.0-xi*xi)*(1.0-eta    )/2.0;
            (*this)(6,0)=(1.0+xi   )*(1.0-eta*eta)/2.0;
            (*this)(7,0)=(1.0-xi*xi)*(1.0+eta    )/2.0;
            (*this)(8,0)=(1.0-xi   )*(1.0-eta*eta)/2.0;

            // derivatives over xi and eta
            (*this)(1,1)=(1.0-eta)*(2.0*xi+eta)/4.0;
            (*this)(1,2)=(1.0-xi )*(xi+2.0*eta)/4.0;

            (*this)(2,1)=(1.0-eta)*(2.0*xi-eta)/4.0;
            (*this)(2,2)=(1.0+xi )*(2.0*eta-xi)/4.0;

            (*this)(3,1)=(1.0+eta)*(2.0*xi+eta)/4.0;
            (*this)(3,2)=(1.0+xi )*(xi+2.0*eta)/4.0;

            (*this)(4,1)=(1.0+eta)*(2.0*xi-eta)/4.0;
            (*this)(4,2)=(1.0-xi )*(2.0*eta-xi)/4.0;

            (*this)(5,1)=xi*(eta-1.0);
            (*this)(5,2)=(xi*xi-1.0)/2.0;

            (*this)(6,1)=(1.0-eta*eta)/2.0;
            (*this)(6,2)=-(1.0+xi)*eta;

            (*this)(7,1)=-xi*(1.0+eta);
            (*this)(7,2)=(1.0-xi*xi)/2.0;

            (*this)(8,1)=(eta*eta-1.0)/2.0;
            (*this)(8,2)=(xi-1.0)*eta;
            break;
        }
        case MeshType::QUAD9:
        {
            // 2D-9Nodes rectangle element
            (*this)(1,0)=(xi*xi-xi )*(eta*eta-eta)/4.0;
            (*this)(2,0)=(xi*xi+xi )*(eta*eta-eta)/4.0;
            (*this)(3,0)=(xi*xi+xi )*(eta*eta+eta)/4.0;
            (*this)(4,0)=(xi*xi-xi )*(eta*eta+eta)/4.0;
            (*this)(5,0)=(1.0-xi*xi)*(eta*eta-eta)/2.0;
            (*this)(6,0)=(xi*xi+xi )*(1.0-eta*eta)/2.0;
            (*this)(7,0)=(1.0-xi*xi)*(eta*eta+eta)/2.0;
            (*this)(8,0)=(xi*xi-xi )*(1.0-eta*eta)/2.0;
            (*this)(9,0)=(1.0-xi*xi)*(1.0-eta*eta);

            (*this)(1,1)=(2.0*xi-1.0)*(eta*eta-eta)/4.0;
            (*this)(1,2)=(xi*xi-xi  )*(2.0*eta-1.0)/4.0;

            (*this)(2,1)=(2.0*xi+1.0)*(eta*eta-eta)/4.0;
            (*this)(2,2)=(xi*xi+xi  )*(2.0*eta-1.0)/4.0;

            (*this)(3,1)=(2.0*xi+1.0)*(eta*eta+eta)/4.0;
            (*this)(3,2)=(xi*xi+xi  )*(2.0*eta+1.0)/4.0;

            (*this)(4,1)=(2.0*xi-1.0)*(eta*eta+eta)/4.0;
            (*this)(4,2)=(xi*xi-xi  )*(2.0*eta+1.0)/4.0;

            (*this)(5,1)=-xi*(eta*eta-eta);
            (*this)(5,2)=(1.0-xi*xi )*(2.0*eta-1.0)/2.0;

            (*this)(6,1)=(2.0*xi+1.0)*(1.0-eta*eta)/2.0;
            (*this)(6,2)=-(xi*xi+xi )*eta;

            (*this)(7,1)=-xi*(eta*eta+eta);
            (*this)(7,2)=(1.0-xi*xi )*(2.0*eta+1.0)/2.0;

            (*this)(8,1)=(2.0*xi-1.0)*(1.0-eta*eta)/2.0;
            (*this)(8,2)=-(xi*xi-xi )*eta;

            (*this)(9,1)=-2.0*xi*(1.0-eta*eta);
            (*this)(9,2)=-2.0*eta*(1.0-xi*xi);
            break;
        }
        case MeshType::TRI3:
        {
            (*this)(1,0)=1.0-xi-eta;
            (*this)(1,1)=-1.0;
            (*this)(1,2)=-1.0;

            (*this)(2,0)=xi;
            (*this)(2,1)=1.0;
            (*this)(2,2)=0.0;

            (*this)(3,0)=eta;
            (*this)(3,1)=0.0;
            (*this)(3,2)=1.0;
            break;
        }
        case MeshType::TRI6:
        {
            // taken from: http://www.sd.ruhr-uni-bochum.de/downloads/Shape_funct.pdf
            (*this)(1,0)= (1.0-xi-eta)*(1.0-2*xi-2*eta);
            (*this)(1,1)=-3.0+4.0*eta+4.0*xi;
            (*this)(1,2)=-3.0+4.0*eta+4.0*xi;

            (*this)(2,0)= xi*(2.0*xi-1.0);
            (*this)(2,1)= 4.0*xi-1.0;
            (*this)(2,2)= 0.0;

            (*this)(3,0)= eta*(2.0*eta-1.0);
            (*this)(3,1)= 0.0;
            (*this)(3,2)= 4.0*eta-1.0;

            (*this)(4,0)= 4.0*xi*(1.0-xi-eta);
            (*this)(4,1)= 4.0*(1.0-eta-2*xi);
            (*this)(4,2)=-4.0*xi;

            (*this)(5,0)= 4.0*xi*eta;
            (*this)(5,1)= 4.0*eta;
            (*this)(5,2)= 4.0*xi;

            (*this)(6,0)= 4.0*eta*(1.0-xi-eta);
            (*this)(6,1)=-4.0*eta;
            (*this)(6,2)= 4.0*(1-2*eta-xi);
            
            break;
        }
        default:
            break;
    }

    _dxdxi=0.0;_dydxi=0.0;_dzdxi=0.0;
    _dxdeta=0.0;_dydeta=0.0;_dzdeta=0.0;
    for(int i=1;i<=GetShapeFunNums();i++){
        _dxdxi+=(*this)(i,1)*nodes(i,1);
        _dydxi+=(*this)(i,1)*nodes(i,2);
        _dzdxi+=(*this)(i,1)*nodes(i,3);

        _dxdeta+=(*this)(i,2)*nodes(i,1);
        _dydeta+=(*this)(i,2)*nodes(i,2);
        _dzdeta+=(*this)(i,2)*nodes(i,3);
    }
   
    _DetJac=(_dydxi*_dzdeta-_dydeta*_dzdxi)*(_dydxi*_dzdeta-_dydeta*_dzdxi)
           +(_dzdxi*_dxdeta-_dzdeta*_dxdxi)*(_dzdxi*_dxdeta-_dzdeta*_dxdxi)
           +(_dxdxi*_dydeta-_dxdeta*_dydxi)*(_dxdxi*_dydeta-_dxdeta*_dydxi);
    _DetJac=sqrt(_DetJac);

    if(abs(_DetJac)<1.0e-15){
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: singular element in 2D case                          !!!   ***\n");
        Msg_AsFem_Exit();
    }

    

    _Jac[0][0]= _dxdxi;_Jac[0][1]= _dydxi;
    _Jac[1][0]=_dxdeta;_Jac[1][1]=_dydeta;

    _XJac[0][0]= _Jac[1][1]/_DetJac;
    _XJac[0][1]=-_Jac[0][1]/_DetJac;
    _XJac[1][0]=-_Jac[1][0]/_DetJac;
    _XJac[1][1]= _Jac[0][0]/_DetJac;

    double temp;
    for(int i=1;i<=GetShapeFunNums();i++){

        if(flag){
            temp=(*this)(i,1)*_XJac[0][0]+(*this)(i,2)*_XJac[0][1];
            (*this)(i,2)=(*this)(i,1)*_XJac[1][0]+(*this)(i,2)*_XJac[1][1];
            (*this)(i,1)=temp;
        }

       _shape_value[i-1]=(*this)(i,0);
       _shape_grad[i-1].setZero();
       _shape_grad[i-1](1)=(*this)(i,1);
       _shape_grad[i-1](2)=(*this)(i,2);
    }
}