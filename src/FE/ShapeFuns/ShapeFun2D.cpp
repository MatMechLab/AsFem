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
//+++ Date   : 2022.05.14
//+++ Purpose: implement the 2d shape fun and its local derivatives
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/ShapeFun2D.h"

ShapeFun2D::ShapeFun2D(){}

void ShapeFun2D::calc2DShapeFun(const MeshType &t_meshtype,
                                const double &xi,const double &eta,
                                const Nodes &t_nodes,const bool &flag,
                                vector<double> &t_vals,
                                vector<Vector3d> &t_ders,
                                double &jacdet){
    if(t_meshtype==MeshType::TRI3){
        ShapeFun2DTri3::calc2DShapeValsAndDerivatives(xi,eta,t_vals,t_ders);
        m_nodes=3;
    }
    else if(t_meshtype==MeshType::TRI6){
        ShapeFun2DTri6::calc2DShapeValsAndDerivatives(xi,eta,t_vals,t_ders);
        m_nodes=6;
    }
    else if(t_meshtype==MeshType::QUAD4){
        ShapeFun2DQuad4::calc2DShapeValsAndDerivatives(xi,eta,t_vals,t_ders);
        m_nodes=4;
    }
    else if(t_meshtype==MeshType::QUAD8){
        ShapeFun2DQuad8::calc2DShapeValsAndDerivatives(xi,eta,t_vals,t_ders);
        m_nodes=8;
    }
    else if(t_meshtype==MeshType::QUAD9){
        ShapeFun2DQuad9::calc2DShapeValsAndDerivatives(xi,eta,t_vals,t_ders);
        m_nodes=9;
    }
    else{
        MessagePrinter::printErrorTxt("unsupported meshtype in calc2DShapeFun of ShapeFun2D.cpp");
        MessagePrinter::exitAsFem();
    }

    if(flag){
        // calculate the derivatives on global coordinates
        m_dxdxi =0.0;m_dydxi =0.0;m_dzdxi =0.0;
        m_dxdeta=0.0;m_dydeta=0.0;m_dzdeta=0.0;
        for(int i=1;i<=m_nodes;i++){
            m_dxdxi+=t_ders[i-1](1)*t_nodes(i,1);
            m_dydxi+=t_ders[i-1](1)*t_nodes(i,2);
            m_dzdxi+=t_ders[i-1](1)*t_nodes(i,3);

            m_dxdeta+=t_ders[i-1](2)*t_nodes(i,1);
            m_dydeta+=t_ders[i-1](2)*t_nodes(i,2);
            m_dzdeta+=t_ders[i-1](2)*t_nodes(i,3);
        }

        // here the jacobian must contains the contribution from z-axis, even though
        // you only have a "2d" element!!!
        jac11=(m_dxdxi*m_dxdxi+m_dydxi*m_dydxi+m_dzdxi*m_dzdxi);
        jac12=(m_dxdxi*m_dxdeta+m_dydxi*m_dydeta+m_dzdxi*m_dzdeta);
        jac21=jac12;
        jac22=(m_dxdeta*m_dxdeta+m_dydeta*m_dydeta+m_dzdeta*m_dzdeta);
        jacdet=sqrt(jac11*jac22-jac12*jac21);// det=sqrt(det(T*T^t))

        if(abs(jacdet)<1.0e-15){
            MessagePrinter::printErrorTxt("singular element in 2d case, error detected in ShapeFun2D");
            MessagePrinter::exitAsFem();
        }
        // for the inverse of jacobian matrix
        xjac11= jac22/(jacdet*jacdet);
        xjac22= jac11/(jacdet*jacdet);
        xjac12=-jac12/(jacdet*jacdet);
        xjac21=-jac21/(jacdet*jacdet);

        m_dxidx=xjac11*m_dxdxi+xjac12*m_dxdeta;
        m_dxidy=xjac11*m_dydxi+xjac12*m_dydeta;
        m_dxidz=xjac11*m_dzdxi+xjac12*m_dzdeta;

        m_detadx=xjac21*m_dxdxi+xjac22*m_dxdeta;
        m_detady=xjac21*m_dydxi+xjac22*m_dydeta;
        m_detadz=xjac21*m_dzdxi+xjac22*m_dzdeta;

        for(int i=1;i<=m_nodes;i++){
            valx=t_ders[i-1](1);
            valy=t_ders[i-1](2);
            t_ders[i-1](1)=valx*m_dxidx+valy*m_detadx;//dN/dx=dN/dxi*dxi/dx+dN/deta*deta/dx
            t_ders[i-1](2)=valx*m_dxidy+valy*m_detady;//dN/dy=dN/dxi*dxi/dy+dN/deta*deta/dy
            t_ders[i-1](3)=valx*m_dxidz+valy*m_detadz;//dN/dz=dN/dxi*dxi/dz+dN/deta*deta/dz
        }
    }
    else{
        jacdet=1.0;
    }
}