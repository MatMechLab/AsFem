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
//+++ Date   : 2022.05.22
//+++ Purpose: defines the user-defined shape functions in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/ShapeFunUser.h"

ShapeFunUser::ShapeFunUser(){
    m_jac.resize(3,3,0.0);
    m_jac.resize(3,3,0.0);
    m_dxdxi=m_dxdeta=m_dxdzeta=0.0;
    m_dydxi=m_dydeta=m_dydzeta=0.0;
    m_dzdxi=m_dzdeta=m_dzdzeta=0.0;
}

void ShapeFunUser::calcUserShapeFun(const ShapeFunType &t_shp_type,
                                    const double &xi,const double &eta,const double &zeta,
                                    const Nodes &t_nodes,const bool &flag,
                                    vector<double> &t_vals,
                                    vector<Vector3d> &t_ders,
                                    double &jacdet){
    if(t_shp_type==ShapeFunType::USER1){
        ShapeFunUser1::calcShapeValsAndDerivatives(xi,eta,zeta,t_vals,t_ders);
        m_nodes=0;// please modify it according to your own shape functions
    }
    else if(t_shp_type==ShapeFunType::USER2){
        ShapeFunUser2::calcShapeValsAndDerivatives(xi,eta,zeta,t_vals,t_ders);
        m_nodes=0;// please modify it according to your own shape functions
    }
    else if(t_shp_type==ShapeFunType::USER3){
        ShapeFunUser3::calcShapeValsAndDerivatives(xi,eta,zeta,t_vals,t_ders);
        m_nodes=0;// please modify it according to your own shape functions
    }
    else if(t_shp_type==ShapeFunType::USER4){
        ShapeFunUser4::calcShapeValsAndDerivatives(xi,eta,zeta,t_vals,t_ders);
        m_nodes=0;// please modify it according to your own shape functions
    }
    else if(t_shp_type==ShapeFunType::USER5){
        ShapeFunUser5::calcShapeValsAndDerivatives(xi,eta,zeta,t_vals,t_ders);
        m_nodes=0;// please modify it according to your own shape functions
    }
    else{
        MessagePrinter::printErrorTxt("unsupported shape fun type in calcUserShapeFun of ShapeFunUser.cpp");
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

        //*****************************************************
        //*** please modifiy the following code according to your 
        //*** shape functions' formula(i.e., 1d, 2d and 3d have different formula!)
        //*** especially the det of your jacobian and the global derivatives!!!
        //*** this is very important to guarantee the correct result!!!
        //*****************************************************
        m_jac(1,1)=  m_dxdxi;m_jac(1,2)=  m_dydxi;m_jac(1,3)=  m_dzdxi;
        m_jac(2,1)= m_dxdeta;m_jac(2,2)= m_dydeta;m_jac(2,3)= m_dzdeta;
        m_jac(3,1)=m_dxdzeta;m_jac(3,2)=m_dydzeta;m_jac(3,3)=m_dzdzeta;

        jacdet=m_jac.det();
        if(abs(jacdet)<1.0e-15){
            MessagePrinter::printErrorTxt("singular element in user-defined case, error detected in ShapeFunUser");
            MessagePrinter::exitAsFem();
        }
        m_xjac=m_jac.inverse();
        for(int i=1;i<=m_nodes;i++){
            // TODO: modify it according to your own formula, for instance, 1d and 2d case definitely use
            // different formula, be careful !!!
            /*
            Hi dude, please modify the following part, your final code should look like below, but not
            exactly this one. Please correct it according to your own formula!!!

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
            */
        }
    }
    else{
        jacdet=1.0;
    }
}