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
//+++ Purpose: implement the general shape function calculation
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/ShapeFun.h"

void ShapeFun::calc(const double &xi,const double &eta,const double &zeta,const Nodes &t_nodes,const bool &flag){
    if(m_shp_type==ShapeFunType::DEFAULT){
        if(m_dim==1){
            calc1DShapeFun(m_mesh_type,xi,t_nodes,flag,m_shpvals,m_shpgrads,m_jacdet);
        }
        else if(m_dim==2){
            calc2DShapeFun(m_mesh_type,xi,eta,t_nodes,flag,m_shpvals,m_shpgrads,m_jacdet);
        }
        else if(m_dim==3){
            calc3DShapeFun(m_mesh_type,xi,eta,zeta,t_nodes,flag,m_shpvals,m_shpgrads,m_jacdet);
        }
        else{
            MessagePrinter::printErrorTxt("dim>3 is not supported for shape function calculation");
            MessagePrinter::exitAsFem();
        }
    }
    else{
        calcUserShapeFun(m_shp_type,xi,eta,zeta,t_nodes,flag,m_shpvals,m_shpgrads,m_jacdet);
    }
}