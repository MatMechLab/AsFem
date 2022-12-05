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
//+++ Purpose: implement the general shape function class
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/ShapeFun.h"


ShapeFun::ShapeFun(){
    m_shp_type=ShapeFunType::DEFAULT;
}

void ShapeFun::init(){
    if(m_shp_type==ShapeFunType::DEFAULT){
        switch (m_mesh_type)
        {
        // for 1d case
        case MeshType::EDGE2:
            m_dim=1;m_funs=2;
            break;
        case MeshType::EDGE3:
            m_dim=1;m_funs=3;
            break;
        case MeshType::EDGE4:
            m_dim=1;m_funs=4;
            break;
        // for 2d case
        case MeshType::TRI3:
            m_dim=2;m_funs=3;
            break;
        case MeshType::TRI6:
            m_dim=2;m_funs=6;
            break;
        case MeshType::QUAD4:
            m_dim=2;m_funs=4;
            break;
        case MeshType::QUAD8:
            m_dim=2;m_funs=8;
            break;
        case MeshType::QUAD9:
            m_dim=2;m_funs=9;
            break;
        // for 3d case
        case MeshType::TET4:
            m_dim=3;m_funs=4;
            break;
        case MeshType::TET10:
            m_dim=3;m_funs=10;
            break;
        case MeshType::HEX8:
            m_dim=3;m_funs=8;
            break;
        case MeshType::HEX20:
            m_dim=3;m_funs=20;
            break;
        case MeshType::HEX27:
            m_dim=3;m_funs=27;
            break;
        default:
            MessagePrinter::printErrorTxt("unsupported mesh type in AsFem, currently we only support "
                                        "1D: edge2, edge3, edge4 shape function, "
                                        "2D: quad4, quad8, quad9 and tri3, tri6 shape function, "
                                        "3D: hex8, hex20, hex27 and tet4, tet10 shape function");
            MessagePrinter::exitAsFem();
            break;
        }
    }
    else{
        MessagePrinter::printWarningTxt("your are trying to use the user-defined shape function, "
                                        "the number of shape funs = "+to_string(m_funs));
    }
    // now we allocate the memory for each array
    m_shpvals.reserve(m_funs);
    m_shpgrads.reserve(m_funs);
    for(int i=0;i<m_funs;i++){
        m_shpvals.push_back(0.0);
        m_shpgrads.push_back(Vector3d(0));
    }
}
//************************************************
void ShapeFun::releaseMemory(){
    m_shpvals.clear();
    m_shpgrads.clear();
}