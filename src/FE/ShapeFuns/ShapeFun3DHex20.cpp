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
//+++ Purpose: 3d hex20 shape function
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/ShapeFun3DHex20.h"

ShapeFun3DHex20::ShapeFun3DHex20(){

}

void ShapeFun3DHex20::calc3DShapeValsAndDerivatives(const double &xi,const double &eta,const double &zeta,
                                                    vector<double> &t_shpvals,
                                                    vector<Vector3d> &t_shpders){
    if(t_shpvals.size()<20 || t_shpders.size()<20){
        MessagePrinter::printErrorTxt("your shape val or derivs vector size is smaller than 20, error detected in ShapeFun3DHex20.cpp");
        MessagePrinter::exitAsFem();
    }
    /**
     * VTK cell type: vtkQuadraticHexahedron
     * bottom layer:
     *  4---11---3
     *  |        |
     * 12        10
     *  |        |
     *  1---9----2
     * middle layer:
     *  20------19
     *   |       |
     *   |       |
     *  17------18
     * top layer:
     *   8---15---7
     *   |        |
     *  16       14
     *   |        |
     *   5---13---6
    */

    t_shpvals[1-1]=(1.0-xi)*(1.0-eta)*(1.0-zeta)*(-2.0-xi-eta-zeta)/8.0;
    t_shpders[1-1](1)= (1.0-eta)*(1.0-zeta)*(2*xi+eta+zeta+1.0)/8.0;
    t_shpders[1-1](2)= (1.0-xi)*(1.0-zeta)*(xi+2.0*eta+zeta+1.0)/8.0;
    t_shpders[1-1](3)= (1.0-xi)*(1.0-eta)*(xi+eta+2.0*zeta+1.0)/8.0;

    t_shpvals[2-1]=(1.0+xi)*(1.0-eta)*(1.0-zeta)*(-2.0+xi-eta-zeta)/8.0;
    t_shpders[2-1](1)=-(1.0-eta)*(1.0-zeta)*(-2.0*xi+eta+zeta+1.0)/8.0;
    t_shpders[2-1](2)=-(1.0+xi)*(1.0-zeta)*(xi-2.0*eta-zeta-1.0)/8.0;
    t_shpders[2-1](3)= (1.0+xi)*(1.0-eta)*(-xi+eta+2.0*zeta+1.0)/8.0;

    t_shpvals[3-1]=(1.0+xi)*(1.0+eta)*(1.0-zeta)*(-2.0+xi+eta-zeta)/8.0;
    t_shpders[3-1](1)= (1.0+eta)*(1.0-zeta)*(2.0*xi+eta-zeta-1.0)/8.0;
    t_shpders[3-1](2)= (1.0+xi)*(1.0-zeta)*(xi+2.0*eta-zeta-1.0)/8.0;
    t_shpders[3-1](3)=-(1.0+xi)*(1.0+eta)*(xi+eta-2.0*zeta-1.0)/8.0;

    t_shpvals[4-1]=(1.0-xi)*(1.0+eta)*(1.0-zeta)*(-2.0-xi+eta-zeta)/8.0;
    t_shpders[4-1](1)=-(1.0+eta)*(1.0-zeta)*(-2.0*xi+eta-zeta-1.0)/8.0;
    t_shpders[4-1](2)=-(1.0-xi)*(1.0-zeta)*(xi-2.0*eta+zeta+1.0)/8.0;
    t_shpders[4-1](3)=-(1.0-xi)*(1.0+eta)*(-xi+eta-2.0*zeta-1.0)/8.0;

    t_shpvals[5-1]=(1.0-xi)*(1.0-eta)*(1.0+zeta)*(-2.0-xi-eta+zeta)/8.0;
    t_shpders[5-1](1)= (1.0-eta)*(1.0+zeta)*(2.0*xi+eta-zeta+1.0)/8.0;
    t_shpders[5-1](2)= (1.0-xi)*(1.0+zeta)*(xi+2.0*eta-zeta+1.0)/8.0;
    t_shpders[5-1](3)=-(1.0-xi)*(1.0-eta)*(xi+eta-2.0*zeta+1.0)/8.0;

    t_shpvals[6-1]=(1.0+xi)*(1.0-eta)*(1.0+zeta)*(-2.0+xi-eta+zeta)/8.0;
    t_shpders[6-1](1)=-(1.0-eta)*(1.0+zeta)*(-2.0*xi+eta-zeta+1.0)/8.0;
    t_shpders[6-1](2)=-(1.0+xi)*(1.0+zeta)*(xi-2.0*eta+zeta-1.0)/8.0;
    t_shpders[6-1](3)=-(1.0+xi)*(1.0-eta)*(-xi+eta-2.0*zeta+1.0)/8.0;

    t_shpvals[7-1]=(1.0+xi)*(1.0+eta)*(1.0+zeta)*(-2.0+xi+eta+zeta)/8.0;
    t_shpders[7-1](1)= (1.0+eta)*(1.0+zeta)*(2*xi+eta+zeta-1.0)/8.0;
    t_shpders[7-1](2)= (1.0+xi)*(1.0+zeta)*(xi+2.0*eta+zeta-1.0)/8.0;
    t_shpders[7-1](3)= (1.0+xi)*(1.0+eta)*(xi+eta+2.0*zeta-1.0)/8.0;

    t_shpvals[8-1]=(1.0-xi)*(1.0+eta)*(1.0+zeta)*(-2.0-xi+eta+zeta)/8.0;
    t_shpders[8-1](1)=-(1.0+eta)*(1.0+zeta)*(-2.0*xi+eta+zeta-1.0)/8.0;
    t_shpders[8-1](2)=-(1.0-xi)*(1.0+zeta)*(xi-2.0*eta-zeta+1.0)/8.0;
    t_shpders[8-1](3)= (1.0-xi)*(1.0+eta)*(-xi+eta+2.0*zeta-1.0)/8.0;

    t_shpvals[9-1]=(1.0-xi*xi)*(1.0-eta)*(1.0-zeta)/4.0;
    t_shpders[9-1](1)=-xi*(1.0-eta)*(1.0-zeta)/2.0;
    t_shpders[9-1](2)=-(1.0-xi*xi)*(1.0-zeta)/4.0;
    t_shpders[9-1](3)=-(1.0-xi*xi)*(1.0-eta)/4.0;

    t_shpvals[10-1]=(1.0+xi)*(1.0-eta*eta)*(1.0-zeta)/4.0;
    t_shpders[10-1](1)= (1.0-eta*eta)*(1.0-zeta)/4.0;
    t_shpders[10-1](2)=-(1.0+xi)*eta*(1.0-zeta)/2.0;
    t_shpders[10-1](3)=-(1.0+xi)*(1.0-eta*eta)/4.0;

    t_shpvals[11-1]=(1.0-xi*xi)*(1.0+eta)*(1.0-zeta)/4.0;
    t_shpders[11-1](1)=-xi*(1.0+eta)*(1.0-zeta)/2.0;
    t_shpders[11-1](2)= (1.0-xi*xi)*(1.0-zeta)/4.0;
    t_shpders[11-1](3)=-(1.0-xi*xi)*(1.0+eta)/4.0;

    t_shpvals[12-1]=(1.0-xi)*(1.0-eta*eta)*(1.0-zeta)/4.0;
    t_shpders[12-1](1)=-(1.0-eta*eta)*(1.0-zeta)/4.0;
    t_shpders[12-1](2)=-(1.0-xi)*eta*(1.0-zeta)/2.0;
    t_shpders[12-1](3)=-(1.0-xi)*(1.0-eta*eta)/4.0;

    t_shpvals[13-1]=(1.0-xi*xi)*(1.0-eta)*(1.0+zeta)/4.0;
    t_shpders[13-1](1)=-xi*(1.0-eta)*(1.0+zeta)/2.0;
    t_shpders[13-1](2)=-(1.0-xi*xi)*(1.0+zeta)/4.0;
    t_shpders[13-1](3)= (1.0-xi*xi)*(1.0-eta)/4.0;

    t_shpvals[14-1]=(1.0+xi)*(1.0-eta*eta)*(1.0+zeta)/4.0;
    t_shpders[14-1](1)= (1.0-eta*eta)*(1.0+zeta)/4.0;
    t_shpders[14-1](2)=-(1.0+xi)*eta*(1.0+zeta)/2.0;
    t_shpders[14-1](3)= (1.0+xi)*(1.0-eta*eta)/4.0;

    t_shpvals[15-1]=(1.0-xi*xi)*(1.0+eta)*(1.0+zeta)/4.0;
    t_shpders[15-1](1)=-xi*(1.0+eta)*(1.0+zeta)/2.0;
    t_shpders[15-1](2)= (1.0-xi*xi)*(1.0+zeta)/4.0;
    t_shpders[15-1](3)= (1.0-xi*xi)*(1.0+eta)/4.0;

    t_shpvals[16-1]=(1.0-xi)*(1.0-eta*eta)*(1.0+zeta)/4.0;
    t_shpders[16-1](1)=-(1.0-eta*eta)*(1.0+zeta)/4.0;
    t_shpders[16-1](2)=-(1.0-xi)*eta*(1.0+zeta)/2.0;
    t_shpders[16-1](3)= (1.0-xi)*(1.0-eta*eta)/4.0;

    t_shpvals[17-1]=(1.0-xi)*(1.0-eta)*(1.0-zeta*zeta)/4.0;
    t_shpders[17-1](1)=-(1.0-eta)*(1.0-zeta*zeta)/4.0;
    t_shpders[17-1](2)=-(1.0-xi)*(1.0-zeta*zeta)/4.0;
    t_shpders[17-1](3)=-(1.0-xi)*(1.0-eta)*zeta/2.0;

    t_shpvals[18-1]=(1.0+xi)*(1.0-eta)*(1.0-zeta*zeta)/4.0;
    t_shpders[18-1](1)= (1.0-eta)*(1.0-zeta*zeta)/4.0;
    t_shpders[18-1](2)=-(1.0+xi)*(1.0-zeta*zeta)/4.0;
    t_shpders[18-1](3)=-(1.0+xi)*(1.0-eta)*zeta/2.0;

    t_shpvals[19-1]=(1.0+xi)*(1.0+eta)*(1.0-zeta*zeta)/4.0;
    t_shpders[19-1](1)= (1.0+eta)*(1.0-zeta*zeta)/4.0;
    t_shpders[19-1](2)= (1.0+xi)*(1.0-zeta*zeta)/4.0;
    t_shpders[19-1](3)=-(1.0+xi)*(1.0+eta)*zeta/2.0;

    t_shpvals[20-1]=(1.0-xi)*(1.0+eta)*(1.0-zeta*zeta)/4.0;
    t_shpders[20-1](1)=-(1.0+eta)*(1.0-zeta*zeta)/4.0;
    t_shpders[20-1](2)= (1.0-xi)*(1.0-zeta*zeta)/4.0;
    t_shpders[20-1](3)=-(1.0-xi)*(1.0+eta)*zeta/2.0;

}