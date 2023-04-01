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
//+++ Purpose: 3d tet10 shape function
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/ShapeFun3DTet10.h"

ShapeFun3DTet10::ShapeFun3DTet10(){

}

void ShapeFun3DTet10::calc3DShapeValsAndDerivatives(const double &xi,const double &eta,const double &zeta,
                                                    vector<double> &t_shpvals,
                                                    vector<Vector3d> &t_shpders){
    if(t_shpvals.size()<10 || t_shpders.size()<10){
        MessagePrinter::printErrorTxt("your shape val or derivs vector size is smaller than 10, error detected in ShapeFun3DTet10.cpp");
        MessagePrinter::exitAsFem();
    }

    // It should be mentioned that, tet10 mesh has different node ordering, here we use the one defined in Gmsh
    // https://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
    // The ordering of the ten points is:ids 1-4 are the four tetra vertices
    //                                   ids 5  is the midedge node between (1,2), 
    //                                   ids 6  is the midedge node between (2,3), 
    //                                   ids 7  is the midedge node between (3,1), 
    //                                   ids 8  is the midedge node between (1,4), 
    //                                   ids 9  is the midedge node between (2,4),
    //                                   ids 10 is the midedge node between (3,4).
    // please do the reording for node-9 and node-10 in msh file reading!!!
    /**
     * VTK Cell type: vtkTetra
     * bottom layer:
     * 3
     * |\
     * | \
     * 7  6
     * |   \
     * |    \
     * 1--5--2
     * 
     * middle layer
     * (3)10(4)
     *     |\
     *     | \
     *     |  \
     *     |   \
     *     |    \
     *     |     \
     *     +------+
     * (1)8(4)  (2)9(4)  
     * 
     * top layer:
     * 4*
    */

    t_shpvals[1-1]=(1.0-xi-eta-zeta)*(1.0-2.0*xi-2*eta-2.0*zeta);
    t_shpders[1-1](1)= 4.0*xi+4.0*eta+4.0*zeta-3.0;
    t_shpders[1-1](2)= 4.0*xi+4.0*eta+4.0*zeta-3.0;
    t_shpders[1-1](3)= 4.0*xi+4.0*eta+4.0*zeta-3.0;

    t_shpvals[2-1]=xi*(2.0*xi-1.0);
    t_shpders[2-1](1)= 4.0*xi-1.0;
    t_shpders[2-1](2)= 0.0;
    t_shpders[2-1](3)= 0.0;

    t_shpvals[3-1]=eta*(2.0*eta-1.0);
    t_shpders[3-1](1)= 0.0;
    t_shpders[3-1](2)= 4.0*eta-1.0;
    t_shpders[3-1](3)= 0.0;

    t_shpvals[4-1]=zeta*(2.0*zeta-1.0);
    t_shpders[4-1](1)= 0.0;
    t_shpders[4-1](2)= 0.0;
    t_shpders[4-1](3)= 4.0*zeta-1.0;

    t_shpvals[5-1]=4.0*(1.0-xi-eta-zeta)*xi;
    t_shpders[5-1](1)=-4.0*(2.0*xi+eta+zeta-1.0);
    t_shpders[5-1](2)=-4.0*xi;
    t_shpders[5-1](3)=-4.0*xi;

    t_shpvals[6-1]=4.0*xi*eta;
    t_shpders[6-1](1)= 4.0*eta;
    t_shpders[6-1](2)= 4.0*xi;
    t_shpders[6-1](3)= 0.0;

    t_shpvals[7-1]=4.0*eta*(1.0-xi-eta-zeta);
    t_shpders[7-1](1)=-4.0*eta;
    t_shpders[7-1](2)=-4.0*(xi+2.0*eta+zeta-1.0);
    t_shpders[7-1](3)=-4.0*eta;

    t_shpvals[8-1]=4.0*(1.0-xi-eta-zeta)*zeta;
    t_shpders[8-1](1)=-4.0*zeta;
    t_shpders[8-1](2)=-4.0*zeta;
    t_shpders[8-1](3)=-4.0*(xi+eta+2.0*zeta-1.0);

    t_shpvals[9-1]=4.0*xi*zeta;
    t_shpders[9-1](1)= 4.0*zeta;
    t_shpders[9-1](2)= 0.0;
    t_shpders[9-1](3)= 4.0*xi;

    t_shpvals[10-1]=4.0*eta*zeta;
    t_shpders[10-1](1)= 0.0;
    t_shpders[10-1](2)= 4.0*zeta;
    t_shpders[10-1](3)= 4.0*eta;

}