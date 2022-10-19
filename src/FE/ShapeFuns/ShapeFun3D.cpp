//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2022.05.15
//+++ Purpose: implement the 3d shape fun and its local derivatives
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/ShapeFun3D.h"

ShapeFun3D::ShapeFun3D(){
    m_jac.resize(3,3,0.0);
    m_jac.resize(3,3,0.0);
    m_dxdxi=m_dxdeta=m_dxdzeta=0.0;
    m_dydxi=m_dydeta=m_dydzeta=0.0;
    m_dzdxi=m_dzdeta=m_dzdzeta=0.0;
}

void ShapeFun3D::calc3DShapeFun(const MeshType &t_meshtype,
                                const double &xi,const double &eta,const double &zeta,
                                const Nodes &t_nodes,const bool &flag,
                                vector<double> &t_vals,
                                vector<Vector3d> &t_ders,
                                double &jacdet){
    if(t_meshtype==MeshType::TET4){
        ShapeFun3DTet4::calc3DShapeValsAndDerivatives(xi,eta,zeta,t_vals,t_ders);
        m_nodes=4;
    }
    else if(t_meshtype==MeshType::TET10){
        ShapeFun3DTet10::calc3DShapeValsAndDerivatives(xi,eta,zeta,t_vals,t_ders);
        m_nodes=10;
    }
    else if(t_meshtype==MeshType::HEX8){
        ShapeFun3DHex8::calc3DShapeValsAndDerivatives(xi,eta,zeta,t_vals,t_ders);
        m_nodes=8;
    }
    else if(t_meshtype==MeshType::HEX20){
        ShapeFun3DHex20::calc3DShapeValsAndDerivatives(xi,eta,zeta,t_vals,t_ders);
        m_nodes=20;
    }
    else if(t_meshtype==MeshType::HEX27){
        ShapeFun3DHex27::calc3DShapeValsAndDerivatives(xi,eta,zeta,t_vals,t_ders);
        m_nodes=27;
    }
    else{
        MessagePrinter::printErrorTxt("unsupported meshtype in calc3DShapeFun of ShapeFun3D.cpp");
        MessagePrinter::exitAsFem();
    }

    if(flag){
        // calculate the derivatives on global coordinates
        m_dxdxi=0.0;m_dydxi=0.0;m_dzdxi=0.0;
        m_dxdeta=0.0;m_dydeta=0.0;m_dzdeta=0.0;
        m_dxdzeta=0.0;m_dydzeta=0.0;m_dzdzeta=0.0;
        for(int i=1;i<=m_nodes;i++){
            m_dxdxi+=t_ders[i-1](1)*t_nodes(i,1);
            m_dydxi+=t_ders[i-1](1)*t_nodes(i,2);
            m_dzdxi+=t_ders[i-1](1)*t_nodes(i,3);

            m_dxdeta+=t_ders[i-1](2)*t_nodes(i,1);
            m_dydeta+=t_ders[i-1](2)*t_nodes(i,2);
            m_dzdeta+=t_ders[i-1](2)*t_nodes(i,3);

            m_dxdzeta+=t_ders[i-1](3)*t_nodes(i,1);
            m_dydzeta+=t_ders[i-1](3)*t_nodes(i,2);
            m_dzdzeta+=t_ders[i-1](3)*t_nodes(i,3);
        }

        m_jac(1,1)=  m_dxdxi;m_jac(1,2)=  m_dydxi;m_jac(1,3)=  m_dzdxi;
        m_jac(2,1)= m_dxdeta;m_jac(2,2)= m_dydeta;m_jac(2,3)= m_dzdeta;
        m_jac(3,1)=m_dxdzeta;m_jac(3,2)=m_dydzeta;m_jac(3,3)=m_dzdzeta;

        jacdet=m_jac.det();

        if(abs(jacdet)<1.0e-15){
            MessagePrinter::printErrorTxt("singular element in 3d case, error detected in ShapeFun3D");
            MessagePrinter::exitAsFem();
        }

        m_xjac=m_jac.inverse();

        for(int i=1;i<=m_nodes;i++){
            temp1 =t_ders[i-1](1)*m_xjac(1,1)
                  +t_ders[i-1](2)*m_xjac(1,2)
                  +t_ders[i-1](3)*m_xjac(1,3);
            temp2 =t_ders[i-1](1)*m_xjac(2,1)
                  +t_ders[i-1](2)*m_xjac(2,2)
                  +t_ders[i-1](3)*m_xjac(2,3);
            temp3 =t_ders[i-1](1)*m_xjac(3,1)
                  +t_ders[i-1](2)*m_xjac(3,2)
                  +t_ders[i-1](3)*m_xjac(3,3);
            
            t_ders[i-1](1) = temp1;
            t_ders[i-1](2) = temp2;
            t_ders[i-1](3) = temp3;
        }
    }
    else{
        jacdet=1.0;
    }
}