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
//+++ Purpose: 3d hex27 shape function
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/ShapeFun3DHex27.h"

ShapeFun3DHex27::ShapeFun3DHex27(){

}

void ShapeFun3DHex27::calc3DShapeValsAndDerivatives(const double &xi,const double &eta,const double &zeta,
                                                    vector<double> &t_shpvals,
                                                    vector<Vector3d> &t_shpders){
    if(t_shpvals.size()<27 || t_shpders.size()<27){
        MessagePrinter::printErrorTxt("your shape val or derivs vector size is smaller than 27, error detected in ShapeFun3DHex27.cpp");
        MessagePrinter::exitAsFem();
    }
    /**
     * VTK cell type: vtkTriQuadraticHexahedron
     * bottom layer:
     *  4--10---3
     *  |   |   |
     * 12--25--10
     *  |   |   |
     *  1---9---2
     * middle layer:
     *  20---24---19
     *   |    |   |
     *  21---27---22
     *   |    |   |
     *  17---23---18
     * top layer:
     *   8--15--7
     *   |   |  |
     *  16--26--14
     *   |   |  |
     *   5--13--6
    */

    t_shpvals[1-1] = xi * (xi - 1) * eta * (eta - 1) * zeta * (zeta - 1) / 8.0;
    t_shpders[1-1](1) = (2 * xi - 1) * eta * (eta - 1) * zeta * (zeta - 1) / 8.0;
    t_shpders[1-1](2) = xi * (xi - 1) * (2 * eta - 1) * zeta * (zeta - 1) / 8.0;
    t_shpders[1-1](3) = xi * (xi - 1) * eta * (eta - 1) * (2 * zeta - 1) / 8.0;

    t_shpvals[2-1] = xi * (xi + 1) * eta * (eta - 1) * zeta * (zeta - 1) / 8.0;
    t_shpders[2-1](1) = (2 * xi + 1) * eta * (eta - 1) * zeta * (zeta - 1) / 8.0;
    t_shpders[2-1](2) = xi * (xi + 1) * (2 * eta - 1) * zeta * (zeta - 1) / 8.0;
    t_shpders[2-1](3) = xi * (xi + 1) * eta * (eta - 1) * (2 * zeta - 1) / 8.0;

    t_shpvals[3-1] = xi * (xi + 1) * eta * (eta + 1) * zeta * (zeta - 1) / 8.0;
    t_shpders[3-1](1) = (2 * xi + 1) * eta * (eta + 1) * zeta * (zeta - 1) / 8.0;
    t_shpders[3-1](2) = xi * (xi + 1) * (2 * eta + 1) * zeta * (zeta - 1) / 8.0;
    t_shpders[3-1](3) = xi * (xi + 1) * eta * (eta + 1) * (2 * zeta - 1) / 8.0;

    t_shpvals[4-1] = xi * (xi - 1) * eta * (eta + 1) * zeta * (zeta - 1) / 8.0;
    t_shpders[4-1](1) = (2 * xi - 1) * eta * (eta + 1) * zeta * (zeta - 1) / 8.0;
    t_shpders[4-1](2) = xi * (xi - 1) * (2 * eta + 1) * zeta * (zeta - 1) / 8.0;
    t_shpders[4-1](3) = xi * (xi - 1) * eta * (eta + 1) * (2 * zeta - 1) / 8.0;

    t_shpvals[5-1] = xi * (xi - 1) * eta * (eta - 1) * zeta * (zeta + 1) / 8.0;
    t_shpders[5-1](1) = (2 * xi - 1) * eta * (eta - 1) * zeta * (zeta + 1) / 8.0;
    t_shpders[5-1](2) = xi * (xi - 1) * (2 * eta - 1) * zeta * (zeta + 1) / 8.0;
    t_shpders[5-1](3) = xi * (xi - 1) * eta * (eta - 1) * (2 * zeta + 1) / 8.0;

    t_shpvals[6-1] = xi * (xi + 1) * eta * (eta - 1) * zeta * (zeta + 1) / 8.0;
    t_shpders[6-1](1) = (2 * xi + 1) * eta * (eta - 1) * zeta * (zeta + 1) / 8.0;
    t_shpders[6-1](2) = xi * (xi + 1) * (2 * eta - 1) * zeta * (zeta + 1) / 8.0;
    t_shpders[6-1](3) = xi * (xi + 1) * eta * (eta - 1) * (2 * zeta + 1) / 8.0;

    t_shpvals[7-1] = xi * (xi + 1) * eta * (eta + 1) * zeta * (zeta + 1) / 8.0;
    t_shpders[7-1](1) = (2 * xi + 1) * eta * (eta + 1) * zeta * (zeta + 1) / 8.0;
    t_shpders[7-1](2) = xi * (xi + 1) * (2 * eta + 1) * zeta * (zeta + 1) / 8.0;
    t_shpders[7-1](3) = xi * (xi + 1) * eta * (eta + 1) * (2 * zeta + 1) / 8.0;

    t_shpvals[8-1] = xi * (xi - 1) * eta * (eta + 1) * zeta * (zeta + 1) / 8.0;
    t_shpders[8-1](1) = (2 * xi - 1) * eta * (eta + 1) * zeta * (zeta + 1) / 8.0;
    t_shpders[8-1](2) = xi * (xi - 1) * (2 * eta + 1) * zeta * (zeta + 1) / 8.0;
    t_shpders[8-1](3) = xi * (xi - 1) * eta * (eta + 1) * (2 * zeta + 1) / 8.0;

    t_shpvals[9-1] = (1 - xi * xi) * eta * (eta - 1) * zeta * (zeta - 1) / 4.0;
    t_shpders[9-1](1) = -xi * eta * (eta - 1) * zeta * (zeta - 1) / 2.0;
    t_shpders[9-1](2) = (1 - xi * xi) * (2 * eta - 1) * zeta * (zeta - 1) / 4.0;
    t_shpders[9-1](3) = (1 - xi * xi) * eta * (eta - 1) * (2 * zeta - 1) / 4.0;

    t_shpvals[10-1] = xi * (xi + 1) * (1 - eta * eta) * zeta * (zeta - 1) / 4.0;
    t_shpders[10-1](1) = (2 * xi + 1) * (1 - eta * eta) * zeta * (zeta - 1) / 4.0;
    t_shpders[10-1](2) = -xi * (xi + 1) * eta * zeta * (zeta - 1) / 2.0;
    t_shpders[10-1](3) = xi * (xi + 1) * (1 - eta * eta) * (2 * zeta - 1) / 4.0;

    t_shpvals[11-1] = (1 - xi * xi) * eta * (eta + 1) * zeta * (zeta - 1) / 4.0;
    t_shpders[11-1](1) = -xi * eta * (eta + 1) * zeta * (zeta - 1) / 2.0;
    t_shpders[11-1](2) = (1 - xi * xi) * (2 * eta + 1) * zeta * (zeta - 1) / 4.0;
    t_shpders[11-1](3) = (1 - xi * xi) * eta * (eta + 1) * (2 * zeta - 1) / 4.0;

    t_shpvals[12-1] = xi * (xi - 1) * (1 - eta * eta) * zeta * (zeta - 1) / 4.0;
    t_shpders[12-1](1) = (2 * xi - 1) * (1 - eta * eta) * zeta * (zeta - 1) / 4.0;
    t_shpders[12-1](2) = -xi * (xi - 1) * eta * zeta * (zeta - 1) / 2.0;
    t_shpders[12-1](3) = xi * (xi - 1) * (1 - eta * eta) * (2 * zeta - 1) / 4.0;

    t_shpvals[13-1] = (1 - xi * xi) * eta * (eta - 1) * zeta * (zeta + 1) / 4.0;
    t_shpders[13-1](1) = -xi * eta * (eta - 1) * zeta * (zeta + 1) / 2.0;
    t_shpders[13-1](2) = (1 - xi * xi) * (2 * eta - 1) * zeta * (zeta + 1) / 4.0;
    t_shpders[13-1](3) = (1 - xi * xi) * eta * (eta - 1) * (2 * zeta + 1) / 4.0;

    t_shpvals[14-1] = xi * (xi + 1) * (1 - eta * eta) * zeta * (zeta + 1) / 4.0;
    t_shpders[14-1](1) = (2 * xi + 1) * (1 - eta * eta) * zeta * (zeta + 1) / 4.0;
    t_shpders[14-1](2) =-xi * (xi + 1) * eta * zeta * (zeta + 1) / 2.0;
    t_shpders[14-1](3) = xi * (xi + 1) * (1 - eta * eta) * (2 * zeta + 1) / 4.0;

    t_shpvals[15-1] = (1 - xi * xi) * eta * (eta + 1) * zeta * (zeta + 1) / 4.0;
    t_shpders[15-1](1) = -xi * eta * (eta + 1) * zeta * (zeta + 1) / 2.0;
    t_shpders[15-1](2) = (1 - xi * xi) * (2 * eta + 1) * zeta * (zeta + 1) / 4.0;
    t_shpders[15-1](3) = (1 - xi * xi) * eta * (eta + 1) * (2 * zeta + 1) / 4.0;

    t_shpvals[16-1] = xi * (xi - 1) * (1 - eta * eta) * zeta * (zeta + 1) / 4.0;
    t_shpders[16-1](1) = (2 * xi - 1) * (1 - eta * eta) * zeta * (zeta + 1) / 4.0;
    t_shpders[16-1](2) = -xi * (xi - 1) * eta * zeta * (zeta + 1) / 2.0;
    t_shpders[16-1](3) = xi * (xi - 1) * (1 - eta * eta) * (2 * zeta + 1) / 4.0;

    t_shpvals[17-1] = xi * (xi - 1) * eta * (eta - 1) * (1 - zeta * zeta) / 4.0;
    t_shpders[17-1](1) = (2 * xi - 1) * eta * (eta - 1) * (1 - zeta * zeta) / 4.0;
    t_shpders[17-1](2) = xi * (xi - 1) * (2 * eta - 1) * (1 - zeta * zeta) / 4.0;
    t_shpders[17-1](3) = -xi * (xi - 1) * eta * (eta - 1) * zeta / 2.0;

    t_shpvals[18-1] = xi * (xi + 1) * eta * (eta - 1) * (1 - zeta * zeta) / 4.0;
    t_shpders[18-1](1) = (2 * xi + 1) * eta * (eta - 1) * (1 - zeta * zeta) / 4.0;
    t_shpders[18-1](2) = xi * (xi + 1) * (2 * eta - 1) * (1 - zeta * zeta) / 4.0;
    t_shpders[18-1](3) = -xi * (xi + 1) * eta * (eta - 1) * zeta / 2.0;

    t_shpvals[19-1] = xi * (xi + 1) * eta * (eta + 1) * (1 - zeta * zeta) / 4.0;
    t_shpders[19-1](1) = (2 * xi + 1) * eta * (eta + 1) * (1 - zeta * zeta) / 4.0;
    t_shpders[19-1](2) = xi * (xi + 1) * (2 * eta + 1) * (1 - zeta * zeta) / 4.0;
    t_shpders[19-1](3) = -xi * (xi + 1) * eta * (eta + 1) * zeta / 2.0;

    t_shpvals[20-1] = xi * (xi - 1) * eta * (eta + 1) * (1 - zeta * zeta) / 4.0;
    t_shpders[20-1](1) = (2 * xi - 1) * eta * (eta + 1) * (1 - zeta * zeta) / 4.0;
    t_shpders[20-1](2) = xi * (xi - 1) * (2 * eta + 1) * (1 - zeta * zeta) / 4.0;
    t_shpders[20-1](3) = -xi * (xi - 1) * eta * (eta + 1) * zeta / 2.0;

    t_shpvals[21-1] = xi * (xi - 1) * (1 - eta * eta) * (1 - zeta * zeta) / 2.0;
    t_shpders[21-1](1) = (2 * xi - 1) * (1 - eta * eta) * (1 - zeta * zeta) / 2.0;
    t_shpders[21-1](2) = -xi * (xi - 1) * eta * (1 - zeta * zeta);
    t_shpders[21-1](3) = -xi * (xi - 1) * (1 - eta * eta) * zeta;

    t_shpvals[22-1] = xi * (xi + 1) * (1 - eta * eta) * (1 - zeta * zeta) / 2.0;
    t_shpders[22-1](1) = (2 * xi + 1) * (1 - eta * eta) * (1 - zeta * zeta) / 2.0;
    t_shpders[22-1](2) = -xi * (xi + 1) * eta * (1 - zeta * zeta);
    t_shpders[22-1](3) = -xi * (xi + 1) * (1 - eta * eta) * zeta;

    t_shpvals[23-1] = (1 - xi * xi) * eta * (eta - 1) * (1 - zeta * zeta) / 2.0;
    t_shpders[23-1](1) = -xi * eta * (eta - 1) * (1 - zeta * zeta);
    t_shpders[23-1](2) = (1 - xi * xi) * (2 * eta - 1) * (1 - zeta * zeta) / 2.0;
    t_shpders[23-1](3) = -(1 - xi * xi) * eta * (eta - 1) * zeta;

    t_shpvals[24-1] = (1 - xi * xi) * eta * (eta + 1) * (1 - zeta * zeta) / 2.0;
    t_shpders[24-1](1) = -xi * eta * (eta + 1) * (1 - zeta * zeta);
    t_shpders[24-1](2) = (1 - xi * xi) * (2 * eta + 1) * (1 - zeta * zeta) / 2.0;
    t_shpders[24-1](3) = -(1 - xi * xi) * eta * (eta + 1) * zeta;

    t_shpvals[25-1] = (1 - xi * xi) * (1 - eta * eta) * zeta * (zeta - 1) / 2.0;
    t_shpders[25-1](1) = -xi * (1 - eta * eta) * zeta * (zeta - 1);
    t_shpders[25-1](2) = -(1 - xi * xi) * eta * zeta * (zeta - 1);
    t_shpders[25-1](3) = (1 - xi * xi) * (1 - eta * eta) * (2 * zeta - 1) / 2.0;

    t_shpvals[26-1] = (1 - xi * xi) * (1 - eta * eta) * zeta * (zeta + 1) / 2.0;
    t_shpders[26-1](1) = -xi * (1 - eta * eta) * zeta * (zeta + 1);
    t_shpders[26-1](2) = -(1 - xi * xi) * eta * zeta * (zeta + 1);
    t_shpders[26-1](3) = (1 - xi * xi) * (1 - eta * eta) * (2 * zeta + 1) / 2.0;

    t_shpvals[27-1] = (1 - xi * xi) * (1 - eta * eta) * (1 - zeta * zeta);
    t_shpders[27-1](1) = -2 * xi * (1 - eta * eta) * (1 - zeta * zeta);
    t_shpders[27-1](2) = -2 * (1 - xi * xi) * eta * (1 - zeta * zeta);
    t_shpders[27-1](3) = -2 * (1 - xi * xi) * (1 - eta * eta) * zeta;

}