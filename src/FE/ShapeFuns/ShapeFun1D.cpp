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
//+++ Purpose: implement the 1d shape fun and its local derivatives
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/ShapeFun1D.h"

ShapeFun1D::ShapeFun1D(){}

void ShapeFun1D::calc1DShapeFun(const MeshType &t_meshtype,const double &xi,
                                const Nodes &t_nodes,const bool &flag,
                                vector<double> &t_vals,
                                vector<Vector3d> &t_ders,
                                double &jacdet){
    if(t_meshtype==MeshType::EDGE2){
        ShapeFun1DEdge2::calc1DShapeValsAndDerivatives(xi,t_vals,t_ders);
        m_nodes=2;
    }
    else if(t_meshtype==MeshType::EDGE3){
        ShapeFun1DEdge3::calc1DShapeValsAndDerivatives(xi,t_vals,t_ders);
        m_nodes=3;
    }
    else if(t_meshtype==MeshType::EDGE4){
        ShapeFun1DEdge4::calc1DShapeValsAndDerivatives(xi,t_vals,t_ders);
        m_nodes=4;
    }
    if(flag){
        // calculate the derivatives on global coordinates
        m_dxdxi=0.0;m_dydxi=0.0;m_dzdxi=0.0;
        for(int i=1;i<=m_nodes;i++){
            m_dxdxi+=t_ders[i-1](1)*t_nodes(i,1);
            m_dydxi+=t_ders[i-1](1)*t_nodes(i,2);
            m_dzdxi+=t_ders[i-1](1)*t_nodes(i,3);
        }
        jacdet=sqrt(m_dxdxi*m_dxdxi+m_dydxi*m_dydxi+m_dzdxi*m_dzdxi);

        if(abs(jacdet)<1.0e-15){
            MessagePrinter::printErrorTxt("singular element in 1d case, error detected in ShapeFun1D");
            MessagePrinter::exitAsFem();
        }
        for(int i=1;i<=m_nodes;i++){
            // here we use the decomposition for each direction, then you will have the derivatives in 3D!
            val=t_ders[i-1](1);
            t_ders[i-1](1)=m_dxdxi*val/(jacdet*jacdet);
            t_ders[i-1](2)=m_dydxi*val/(jacdet*jacdet);
            t_ders[i-1](3)=m_dzdxi*val/(jacdet*jacdet);
        }
    }
    else{
        jacdet=1.0;
    }
}