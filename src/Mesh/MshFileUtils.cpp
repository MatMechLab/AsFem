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
//+++ Date   : 2022.08.16
//+++ Purpose: This class offers the utils function for msh file
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Mesh/MshFileUtils.h"

MshFileUtils::MshFileUtils(){

}

int MshFileUtils::getElmtNodesNumFromElmtType(const int &elmttype){
    switch(elmttype){
        case 1:
            // 2-node line
            return 2;
        case 2:
            // 3-node triangle
            return 3;
        case 3:
            // 4-node quadrangle
            return 4;
        case 4:
            // 4-node tetrahedron
            return 4;
        case 5:
            // 8-node hexahedron
            return 8;
        case 6:
            // 6-node prism
            return 6;
        case 7:
            // 5-node pyramid
            return 5;
        case 8:
            //3-node second order line
            return 3;
        case 9:
            // 6-node second order triangle
            return 6;
        case 10:
            // 9-node second order quadrangle
            return 9;
        case 11:
            // 10-node second order tetrahedron
            return 10;
        case 12:
            // 27-node second order hexahedron
            return 27;
        case 13:
            // 18-node second order prism
            return 18;
        case 14:
            // 14-node second order pyramid
            return 14;
        case 15:
            // 1-node point
            return 1;
        case 16:
            // 8-node second order quadrangle
            return 8;
        case 17:
            // 20-node second order hexahedron
            return 20;
        case 18:
            // 15-node second order prism
            return 15;
        case 19:
            // 13-node second order pyramid
            return 13;
        case 20:
            // 9-node third order incomplete triangle
            return 9;
        case 21:
            // 10-node third order triangle
            return 10;
        case 22:
            // 12-node fourth order incomplete triangle
            return 12;
        case 23:
            // 15-node fourth order triangle
            return 15;
        case 24:
            // 15-node fifth order incomplete triangle
            return 15;
        case 25:
            // 21-node fifth order complete triangle
            return 21;
        case 26:
            // 4-node third order edge
            return 4;
        case 27:
            // 5-node fourth order edge
            return 5;
        case 28:
            // 6-node fifth order edge
            return 6;
        case 29:
            // 20-node third order tetrahedron
            return 20;
        case 30:
            // 35-node fourth order tetrahedron
            return 35;
        case 31:
            // 56-node fifth order tetrahedron
            return 56;
        case 92:
            // 64-node third order hexahedron
            return 64;
        case 93:
            // 125-node fourth order hexahedron
            return 125;
        default:
            MessagePrinter::printErrorTxt("Unsupported gmsh element type in getElmtNodesNumFromElmtType, please check your code or your mesh-importer");
            MessagePrinter::exitAsFem();
            return -1;
    }
    return -1;
}
//************************************************************
int MshFileUtils::getElmtDimFromElmtType(const int &elmttype){
    switch(elmttype){
        case 1:
            // 2-node line
            return 1;
        case 2:
            // 3-node triangle
            return 2;
        case 3:
            // 4-node quadrangle
            return 2;
        case 4:
            // 4-node tetrahedron
            return 3;
        case 5:
            // 8-node hexahedron
            return 3;
        case 6:
            // 6-node prism
            return 3;
        case 7:
            // 5-node pyramid
            return 3;
        case 8:
            //3-node second order line
            return 1;
        case 9:
            // 6-node second order triangle
            return 2;
        case 10:
            // 9-node second order quadrangle
            return 2;
        case 11:
            // 10-node second order tetrahedron
            return 3;
        case 12:
            // 27-node second order hexahedron
            return 3;
        case 13:
            // 18-node second order prism
            return 3;
        case 14:
            // 14-node second order pyramid
            return 3;
        case 15:
            // 1-node point
            return 0;
        case 16:
            // 8-node second order quadrangle
            return 2;
        case 17:
            // 20-node second order hexahedron
            return 3;
        case 18:
            // 15-node second order prism
            return 3;
        case 19:
            // 13-node second order pyramid
            return 3;
        case 20:
            // 9-node third order incomplete triangle
            return 2;
        case 21:
            // 10-node third order triangle
            return 2;
        case 22:
            // 12-node fourth order incomplete triangle
            return 2;
        case 23:
            // 15-node fourth order triangle
            return 2;
        case 24:
            // 15-node fifth order incomplete triangle
            return 2;
        case 25:
            // 21-node fifth order complete triangle
            return 2;
        case 26:
            // 4-node third order edge
            return 1;
        case 27:
            // 5-node fourth order edge
            return 1;
        case 28:
            // 6-node fifth order edge
            return 1;
        case 29:
            // 20-node third order tetrahedron
            return 3;
        case 30:
            // 35-node fourth order tetrahedron
            return 3;
        case 31:
            // 56-node fifth order tetrahedron
            return 3;
        case 92:
            // 64-node third order hexahedron
            return 3;
        case 93:
            // 125-node fourth order hexahedron
            return 3;
        default:
            MessagePrinter::printErrorTxt("Unsupported gmsh element type in getElmtDimFromElmtType, please check your code or your mesh-importer");
            MessagePrinter::exitAsFem();
            return -1;
    }
    return -1;
}
//*********************************************************
int MshFileUtils::getElmtOrderFromElmtType(const int &elmttype){
    switch(elmttype){
        case 1:
            // 2-node line
            return 1;
        case 2:
            // 3-node triangle
            return 1;
        case 3:
            // 4-node quadrangle
            return 1;
        case 4:
            // 4-node tetrahedron
            return 1;
        case 5:
            // 8-node hexahedron
            return 1;
        case 6:
            // 6-node prism
            return 1;
        case 7:
            // 5-node pyramid
            return 1;
        case 8:
            //3-node second order line
            return 2;
        case 9:
            // 6-node second order triangle
            return 2;
        case 10:
            // 9-node second order quadrangle
            return 2;
        case 11:
            // 10-node second order tetrahedron
            return 2;
        case 12:
            // 27-node second order hexahedron
            return 2;
        case 13:
            // 18-node second order prism
            return 2;
        case 14:
            // 14-node second order pyramid
            return 2;
        case 15:
            // 1-node point
            return 0;
        case 16:
            // 8-node second order quadrangle
            return 2;
        case 17:
            // 20-node second order hexahedron
            return 2;
        case 18:
            // 15-node second order prism
            return 2;
        case 19:
            // 13-node second order pyramid
            return 2;
        case 20:
            // 9-node third order incomplete triangle
            return 3;
        case 21:
            // 10-node third order triangle
            return 3;
        case 22:
            // 12-node fourth order incomplete triangle
            return 4;
        case 23:
            // 15-node fourth order triangle
            return 4;
        case 24:
            // 15-node fifth order incomplete triangle
            return 5;
        case 25:
            // 21-node fifth order complete triangle
            return 5;
        case 26:
            // 4-node third order edge
            return 3;
        case 27:
            // 5-node fourth order edge
            return 4;
        case 28:
            // 6-node fifth order edge
            return 5;
        case 29:
            // 20-node third order tetrahedron
            return 3;
        case 30:
            // 35-node fourth order tetrahedron
            return 4;
        case 31:
            // 56-node fifth order tetrahedron
            return 5;
        case 92:
            // 64-node third order hexahedron
            return 3;
        case 93:
            // 125-node fourth order hexahedron
            return 4;
        default:
            MessagePrinter::printErrorTxt("Unsupported gmsh element type in getElmtOrderFromElmtType, please check your code or your mesh-importer");
            MessagePrinter::exitAsFem();
            return -1;
    }
    return -1;
}
//*********************************************************
int MshFileUtils::getElmtVTKCellTypeFromElmtType(const int &elmttype){
    switch(elmttype){
        case 1:
            // 2-node line
            return 3;
        case 2:
            // 3-node triangle
            return 5;
        case 3:
            // 4-node quadrangle
            return 9;
        case 4:
            // 4-node tetrahedron
            return 10;
        case 5:
            // 8-node hexahedron
            return 12;
        case 6:
            // 6-node prism
            return 13;
        case 7:
            // 5-node pyramid
            return 14;
        case 8:
            //3-node second order line
            return 4;
        case 9:
            // 6-node second order triangle
            return 22;
        case 10:
            // 9-node second order quadrangle
            return 28;
        case 11:
            // 10-node second order tetrahedron
            return 24;
        case 12:
            // 27-node second order hexahedron
            return 29;
        case 13:
            // 18-node second order prism
            return 6;
        case 14:
            // 14-node second order pyramid
            return 6;
        case 15:
            // 1-node point
            return 0;
        case 16:
            // 8-node second order quadrangle
            return 23;
        case 17:
            // 20-node second order hexahedron
            return 25;
        case 18:
            // 15-node second order prism
            return 5;
        case 19:
            // 13-node second order pyramid
            return 6;
        case 20:
            // 9-node third order incomplete triangle
            return 0;
        case 21:
            // 10-node third order triangle
            return 0;
        case 22:
            // 12-node fourth order incomplete triangle
            return 0;
        case 23:
            // 15-node fourth order triangle
            return 0;
        case 24:
            // 15-node fifth order incomplete triangle
            return 0;
        case 25:
            // 21-node fifth order complete triangle
            return 0;
        case 26:
            // 4-node third order edge
            return 4;
        case 27:
            // 5-node fourth order edge
            return 4;
        case 28:
            // 6-node fifth order edge
            return 4;
        case 29:
            // 20-node third order tetrahedron
            return 4;
        case 30:
            // 35-node fourth order tetrahedron
            return 4;
        case 31:
            // 56-node fifth order tetrahedron
            return 4;
        case 92:
            // 64-node third order hexahedron
            return 6;
        case 93:
            // 125-node fourth order hexahedron
            return 6;
        default:
            MessagePrinter::printErrorTxt("Unsupported gmsh element type in getElmtVTKCellTypeFromElmtType, please check your code or your mesh-importer");
            MessagePrinter::exitAsFem();
            return 0;
    }
    return -1;
}
//******************************************************************
MeshType MshFileUtils::getElmtMeshTypeFromElmtType(const int &elmttype){
    switch(elmttype){
        case 1:
            // 2-node line
            return MeshType::EDGE2;
        case 2:
            // 3-node triangle
            return MeshType::TRI3;
        case 3:
            // 4-node quadrangle
            return MeshType::QUAD4;
        case 4:
            // 4-node tetrahedron
            return MeshType::TET4;
        case 5:
            // 8-node hexahedron
            return MeshType::HEX8;
        case 6:
            // 6-node prism
            return MeshType::NULLTYPE;
        case 7:
            // 5-node pyramid
            return MeshType::NULLTYPE;
        case 8:
            //3-node second order line
            return MeshType::EDGE3;
        case 9:
            // 6-node second order triangle
            return MeshType::TRI6;
        case 10:
            // 9-node second order quadrangle
            return MeshType::QUAD9;
        case 11:
            // 10-node second order tetrahedron
            return MeshType::TET10;
        case 12:
            // 27-node second order hexahedron
            return MeshType::HEX27;
        case 13:
            // 18-node second order prism
            return MeshType::NULLTYPE;
        case 14:
            // 14-node second order pyramid
            return MeshType::NULLTYPE;
        case 15:
            // 1-node point
            return MeshType::POINT;
        case 16:
            // 8-node second order quadrangle
            return MeshType::QUAD8;
        case 17:
            // 20-node second order hexahedron
            return MeshType::HEX20;
        case 18:
            // 15-node second order prism
            return MeshType::NULLTYPE;
        case 19:
            // 13-node second order pyramid
            return MeshType::NULLTYPE;
        case 20:
            // 9-node third order incomplete triangle
            return MeshType::NULLTYPE;
        case 21:
            // 10-node third order triangle
            return MeshType::NULLTYPE;
        case 22:
            // 12-node fourth order incomplete triangle
            return MeshType::NULLTYPE;
        case 23:
            // 15-node fourth order triangle
            return MeshType::NULLTYPE;
        case 24:
            // 15-node fifth order incomplete triangle
            return MeshType::NULLTYPE;
        case 25:
            // 21-node fifth order complete triangle
            return MeshType::NULLTYPE;
        case 26:
            // 4-node third order edge
            return MeshType::EDGE4;
        case 27:
            // 5-node fourth order edge
            return MeshType::EDGE5;
        case 28:
            // 6-node fifth order edge
            return MeshType::NULLTYPE;
        case 29:
            // 20-node third order tetrahedron
            return MeshType::NULLTYPE;
        case 30:
            // 35-node fourth order tetrahedron
            return MeshType::NULLTYPE;
        case 31:
            // 56-node fifth order tetrahedron
            return MeshType::NULLTYPE;
        case 92:
            // 64-node third order hexahedron
            return MeshType::NULLTYPE;
        case 93:
            // 125-node fourth order hexahedron
            return MeshType::NULLTYPE;
        default:
            MessagePrinter::printErrorTxt("Unsupported gmsh element type in getElmtMeshTypeFromElmtType, please check your code or your mesh-importer");
            MessagePrinter::exitAsFem();
            return MeshType::NULLTYPE;
    }
    return MeshType::NULLTYPE;
}
//*******************************************************************
string MshFileUtils::getElmtMeshTypeNameFromElmtType(const int &elmttype){
    switch(elmttype){
        case 1:
            // 2-node line
            return "edge2";
        case 2:
            // 3-node triangle
            return "tri3";
        case 3:
            // 4-node quadrangle
            return "quad4";
        case 4:
            // 4-node tetrahedron
            return "tet4";
        case 5:
            // 8-node hexahedron
            return "hex8";
        case 6:
            // 6-node prism
            return "prims6";
        case 7:
            // 5-node pyramid
            return "pyramid5";
        case 8:
            //3-node second order line
            return "edge3";
        case 9:
            // 6-node second order triangle
            return "tri6";
        case 10:
            // 9-node second order quadrangle
            return "quad9";
        case 11:
            // 10-node second order tetrahedron
            return "tet10";
        case 12:
            // 27-node second order hexahedron
            return "hex27";
        case 13:
            // 18-node second order prism
            return "prism18";
        case 14:
            // 14-node second order pyramid
            return "pyramid14";
        case 15:
            // 1-node point
            return "point";
        case 16:
            // 8-node second order quadrangle
            return "quad8";
        case 17:
            // 20-node second order hexahedron
            return "hex20";
        case 18:
            // 15-node second order prism
            return "prism15";
        case 19:
            // 13-node second order pyramid
            return "pyramid13";
        case 20:
            // 9-node third order incomplete triangle
            return "tri9";
        case 21:
            // 10-node third order triangle
            return "tri10";
        case 22:
            // 12-node fourth order incomplete triangle
            return "tri12";
        case 23:
            // 15-node fourth order triangle
            return "tri15";
        case 24:
            // 15-node fifth order incomplete triangle
            return "tri15";
        case 25:
            // 21-node fifth order complete triangle
            return "tri21";
        case 26:
            // 4-node third order edge
            return "edge4";
        case 27:
            // 5-node fourth order edge
            return "edge5";
        case 28:
            // 6-node fifth order edge
            return "edge6";
        case 29:
            // 20-node third order tetrahedron
            return "tet20";
        case 30:
            // 35-node fourth order tetrahedron
            return "tet35";
        case 31:
            // 56-node fifth order tetrahedron
            return "tet56";
        case 92:
            // 64-node third order hexahedron
            return "hex64";
        case 93:
            // 125-node fourth order hexahedron
            return "hex125";
        default:
            MessagePrinter::printErrorTxt("Unsupported gmsh element type in getElmtMeshTypeNameFromElmtType, please check your code or your mesh-importer");
            MessagePrinter::exitAsFem();
            return "unknown";
    }
    return "unknown";
}
//*****************************************************************
//*** for lower dimension element (mesh)
//*****************************************************************
int MshFileUtils::getSubElmtNodesNumFromElmtType(const int &elmttype){
    switch(elmttype){
        case 1:
            // 2-node line
            return 1;
        case 2:
            // 3-node triangle
            return 2;
        case 3:
            // 4-node quadrangle
            return 2;
        case 4:
            // 4-node tetrahedron
            return 3;
        case 5:
            // 8-node hexahedron
            return 4;
        case 6:
            // 6-node prism
            return -1;
        case 7:
            // 5-node pyramid
            return -1;
        case 8:
            //3-node second order line
            return 1;
        case 9:
            // 6-node second order triangle
            return 3;
        case 10:
            // 9-node second order quadrangle
            return 3;
        case 11:
            // 10-node second order tetrahedron
            return 6;
        case 12:
            // 27-node second order hexahedron
            return 9;
        case 13:
            // 18-node second order prism
            return -1;
        case 14:
            // 14-node second order pyramid
            return -1;
        case 15:
            // 1-node point
            return 0;
        case 16:
            // 8-node second order quadrangle
            return 3;
        case 17:
            // 20-node second order hexahedron
            return 8;
        case 18:
            // 15-node second order prism
            return -1;
        case 19:
            // 13-node second order pyramid
            return -1;
        case 20:
            // 9-node third order incomplete triangle
            return -1;
        case 21:
            // 10-node third order triangle
            return -1;
        case 22:
            // 12-node fourth order incomplete triangle
            return -1;
        case 23:
            // 15-node fourth order triangle
            return -1;
        case 24:
            // 15-node fifth order incomplete triangle
            return -1;
        case 25:
            // 21-node fifth order complete triangle
            return -1;
        case 26:
            // 4-node third order edge
            return 1;
        case 27:
            // 5-node fourth order edge
            return 1;
        case 28:
            // 6-node fifth order edge
            return 1;
        case 29:
            // 20-node third order tetrahedron
            return -1;
        case 30:
            // 35-node fourth order tetrahedron
            return -1;
        case 31:
            // 56-node fifth order tetrahedron
            return -1;
        case 92:
            // 64-node third order hexahedron
            return -1;
        case 93:
            // 125-node fourth order hexahedron
            return -1;
        default:
            MessagePrinter::printErrorTxt("Unsupported gmsh element type in getElmtMeshTypeNameFromElmtType, please check your code or your mesh-importer");
            MessagePrinter::exitAsFem();
            return -1;
    }
}
//********************************************************
int MshFileUtils::getSubElmtDimFromElmtType(const int &elmttype){
    switch(elmttype){
        case 1:
            // 2-node line
            return 1-1;
        case 2:
            // 3-node triangle
            return 2-1;
        case 3:
            // 4-node quadrangle
            return 2-1;
        case 4:
            // 4-node tetrahedron
            return 3-1;
        case 5:
            // 8-node hexahedron
            return 3-1;
        case 6:
            // 6-node prism
            return 3-1;
        case 7:
            // 5-node pyramid
            return 3-1;
        case 8:
            //3-node second order line
            return 1-1;
        case 9:
            // 6-node second order triangle
            return 2-1;
        case 10:
            // 9-node second order quadrangle
            return 2-1;
        case 11:
            // 10-node second order tetrahedron
            return 3-1;
        case 12:
            // 27-node second order hexahedron
            return 3-1;
        case 13:
            // 18-node second order prism
            return 3-1;
        case 14:
            // 14-node second order pyramid
            return 3-1;
        case 15:
            // 1-node point
            return 0-1;
        case 16:
            // 8-node second order quadrangle
            return 2-1;
        case 17:
            // 20-node second order hexahedron
            return 3-1;
        case 18:
            // 15-node second order prism
            return 3-1;
        case 19:
            // 13-node second order pyramid
            return 3-1;
        case 20:
            // 9-node third order incomplete triangle
            return 2-1;
        case 21:
            // 10-node third order triangle
            return 2-1;
        case 22:
            // 12-node fourth order incomplete triangle
            return 2-1;
        case 23:
            // 15-node fourth order triangle
            return 2-1;
        case 24:
            // 15-node fifth order incomplete triangle
            return 2-1;
        case 25:
            // 21-node fifth order complete triangle
            return 2-1;
        case 26:
            // 4-node third order edge
            return 1-1;
        case 27:
            // 5-node fourth order edge
            return 1-1;
        case 28:
            // 6-node fifth order edge
            return 1-1;
        case 29:
            // 20-node third order tetrahedron
            return 3-1;
        case 30:
            // 35-node fourth order tetrahedron
            return 3-1;
        case 31:
            // 56-node fifth order tetrahedron
            return 3-1;
        case 92:
            // 64-node third order hexahedron
            return 3-1;
        case 93:
            // 125-node fourth order hexahedron
            return 3-1;
        default:
            MessagePrinter::printErrorTxt("Unsupported gmsh element type in getSubElmtDimFromElmtType, please check your code or your mesh-importer");
            MessagePrinter::exitAsFem();
            return -1;
    }
    return -1;
}
//**************************************************************
int MshFileUtils::getSubElmtOrderFromElmtType(const int &elmttype){
    switch(elmttype){
        case 1:
            // 2-node line
            return 1;
        case 2:
            // 3-node triangle
            return 1;
        case 3:
            // 4-node quadrangle
            return 1;
        case 4:
            // 4-node tetrahedron
            return 1;
        case 5:
            // 8-node hexahedron
            return 1;
        case 6:
            // 6-node prism
            return 1;
        case 7:
            // 5-node pyramid
            return 1;
        case 8:
            //3-node second order line
            return 1;
        case 9:
            // 6-node second order triangle
            return 2;
        case 10:
            // 9-node second order quadrangle
            return 2;
        case 11:
            // 10-node second order tetrahedron
            return 2;
        case 12:
            // 27-node second order hexahedron
            return 2;
        case 13:
            // 18-node second order prism
            return 2;
        case 14:
            // 14-node second order pyramid
            return 2;
        case 15:
            // 1-node point
            return 0;
        case 16:
            // 8-node second order quadrangle
            return 2;
        case 17:
            // 20-node second order hexahedron
            return 2;
        case 18:
            // 15-node second order prism
            return 2;
        case 19:
            // 13-node second order pyramid
            return 2;
        case 20:
            // 9-node third order incomplete triangle
            return 3;
        case 21:
            // 10-node third order triangle
            return 3;
        case 22:
            // 12-node fourth order incomplete triangle
            return 4;
        case 23:
            // 15-node fourth order triangle
            return 4;
        case 24:
            // 15-node fifth order incomplete triangle
            return 5;
        case 25:
            // 21-node fifth order complete triangle
            return 5;
        case 26:
            // 4-node third order edge
            return 1;
        case 27:
            // 5-node fourth order edge
            return 1;
        case 28:
            // 6-node fifth order edge
            return 1;
        case 29:
            // 20-node third order tetrahedron
            return 3;
        case 30:
            // 35-node fourth order tetrahedron
            return 4;
        case 31:
            // 56-node fifth order tetrahedron
            return 5;
        case 92:
            // 64-node third order hexahedron
            return 3;
        case 93:
            // 125-node fourth order hexahedron
            return 4;
        default:
            MessagePrinter::printErrorTxt("Unsupported gmsh element type in getSubElmtOrderFromElmtType, please check your code or your mesh-importer");
            MessagePrinter::exitAsFem();
            return -1;
    }
    return -1;
}
//***************************************************************
MeshType MshFileUtils::getSubElmtMeshTypeFromElmtType(const int &elmttype){
    switch(elmttype){
        case 1:
            // 2-node line
            return MeshType::POINT;
        case 2:
            // 3-node triangle
            return MeshType::EDGE2;
        case 3:
            // 4-node quadrangle
            return MeshType::EDGE2;
        case 4:
            // 4-node tetrahedron
            return MeshType::TRI3;
        case 5:
            // 8-node hexahedron
            return MeshType::QUAD4;
        case 6:
            // 6-node prism
            return MeshType::NULLTYPE;
        case 7:
            // 5-node pyramid
            return MeshType::NULLTYPE;
        case 8:
            //3-node second order line
            return MeshType::POINT;
        case 9:
            // 6-node second order triangle
            return MeshType::EDGE3;
        case 10:
            // 9-node second order quadrangle
            return MeshType::EDGE3;
        case 11:
            // 10-node second order tetrahedron
            return MeshType::TRI6;
        case 12:
            // 27-node second order hexahedron
            return MeshType::QUAD9;
        case 13:
            // 18-node second order prism
            return MeshType::NULLTYPE;
        case 14:
            // 14-node second order pyramid
            return MeshType::NULLTYPE;
        case 15:
            // 1-node point
            return MeshType::NULLTYPE;
        case 16:
            // 8-node second order quadrangle
            return MeshType::EDGE3;
        case 17:
            // 20-node second order hexahedron
            return MeshType::QUAD8;
        case 18:
            // 15-node second order prism
            return MeshType::NULLTYPE;
        case 19:
            // 13-node second order pyramid
            return MeshType::NULLTYPE;
        case 20:
            // 9-node third order incomplete triangle
            return MeshType::NULLTYPE;
        case 21:
            // 10-node third order triangle
            return MeshType::NULLTYPE;
        case 22:
            // 12-node fourth order incomplete triangle
            return MeshType::NULLTYPE;
        case 23:
            // 15-node fourth order triangle
            return MeshType::NULLTYPE;
        case 24:
            // 15-node fifth order incomplete triangle
            return MeshType::NULLTYPE;
        case 25:
            // 21-node fifth order complete triangle
            return MeshType::NULLTYPE;
        case 26:
            // 4-node third order edge
            return MeshType::POINT;
        case 27:
            // 5-node fourth order edge
            return MeshType::POINT;
        case 28:
            // 6-node fifth order edge
            return MeshType::POINT;
        case 29:
            // 20-node third order tetrahedron
            return MeshType::NULLTYPE;
        case 30:
            // 35-node fourth order tetrahedron
            return MeshType::NULLTYPE;
        case 31:
            // 56-node fifth order tetrahedron
            return MeshType::NULLTYPE;
        case 92:
            // 64-node third order hexahedron
            return MeshType::NULLTYPE;
        case 93:
            // 125-node fourth order hexahedron
            return MeshType::NULLTYPE;
        default:
            MessagePrinter::printErrorTxt("Unsupported gmsh element type in getSubElmtMeshTypeFromElmtType, please check your code or your mesh-importer");
            MessagePrinter::exitAsFem();
            return MeshType::NULLTYPE;
    }
    return MeshType::NULLTYPE;
}
//*********************************************************
//*** for node reordering
//*********************************************************
void MshFileUtils::reorderNodesIndex(const int &elmttype,vector<int> &elmtconn){
    switch(elmttype){
        case 1:
            // 2-node line
            return;
        case 2:
            // 3-node triangle
            return;
        case 3:
            // 4-node quadrangle
            return;
        case 4:
            // 4-node tetrahedron
            return;
        case 5:
            // 8-node hexahedron
            return;
        case 6:
            // 6-node prism
            return;
        case 7:
            // 5-node pyramid
            return;
        case 8:{
            //3-node second order line
            // change 0--2--1 to 0--1--2
            int i=elmtconn[1];
            elmtconn[1]=elmtconn[2];
            elmtconn[2]=i;
            }
            return;
        case 9:
            // 6-node second order triangle
            return;
        case 10:
            // 9-node second order quadrangle
            return;
        case 11:{
            // 10-node second order tetrahedron
            int i=elmtconn[9-1];
            elmtconn[9-1]=elmtconn[10-1];
            elmtconn[10-1]=i;
            }
            return;
        case 12:
            // 27-node second order hexahedron
            return;
        case 13:
            // 18-node second order prism
            return;
        case 14:
            // 14-node second order pyramid
            return;
        case 15:
            // 1-node point
            return;
        case 16:
            // 8-node second order quadrangle
            return;
        case 17:
            // 20-node second order hexahedron
            return;
        case 18:
            // 15-node second order prism
            return;
        case 19:
            // 13-node second order pyramid
            return;
        case 20:
            // 9-node third order incomplete triangle
            return;
        case 21:
            // 10-node third order triangle
            return;
        case 22:
            // 12-node fourth order incomplete triangle
            return;
        case 23:
            // 15-node fourth order triangle
            return;
        case 24:
            // 15-node fifth order incomplete triangle
            return;
        case 25:
            // 21-node fifth order complete triangle
            return;
        case 26:{
            // 4-node third order edge
            // 0---2---3---1 to 0---1---2---3
            int i2=elmtconn[4-1];
            int i3=elmtconn[2-1];
            int i4=elmtconn[3-1];
            elmtconn[2-1]=i2;
            elmtconn[3-1]=i3;
            elmtconn[4-1]=i4;
            }
            return;
        case 27:
            // 5-node fourth order edge
            return;
        case 28:
            // 6-node fifth order edge
            return;
        case 29:
            // 20-node third order tetrahedron
            return;
        case 30:
            // 35-node fourth order tetrahedron
            return;
        case 31:
            // 56-node fifth order tetrahedron
            return;
        case 92:
            // 64-node third order hexahedron
            return;
        case 93:
            // 125-node fourth order hexahedron
            return;
        default:
            MessagePrinter::printErrorTxt("Unsupported gmsh element type in reorderNodesIndex, please check your code or your mesh-importer");
            MessagePrinter::exitAsFem();
            return;
    }
    
}
void MshFileUtils::reorderGmsh2NodesIndex(const int &elmttype,vector<int> &elmtconn){
    switch(elmttype){
        case 1:
            // 2-node line
            return;
        case 2:
            // 3-node triangle
            return;
        case 3:
            // 4-node quadrangle
            return;
        case 4:{
            // 4-node tetrahedron
            int i=elmtconn[1];
            elmtconn[1]=elmtconn[2];
            elmtconn[2]=i;
            }
            return;
        case 5:
            // 8-node hexahedron
            return;
        case 6:
            // 6-node prism
            return;
        case 7:
            // 5-node pyramid
            return;
        case 8:{
            //3-node second order line
            // change 0--2--1 to 0--1--2
            int i=elmtconn[1];
            elmtconn[1]=elmtconn[2];
            elmtconn[2]=i;
            }
            return;
        case 9:
            // 6-node second order triangle
            return;
        case 10:
            // 9-node second order quadrangle
            return;
        case 11:{
            // 10-node second order tetrahedron
            int i2=elmtconn[4-1];
            int i4=elmtconn[2-1];
            int i5=elmtconn[8-1];
            int i6=elmtconn[9-1];
            int i8=elmtconn[5-1];
            int i9=elmtconn[10-1];
            int i10=elmtconn[6-1];
            elmtconn[2-1]=i2;
            elmtconn[4-1]=i4;
            elmtconn[5-1]=i5;
            elmtconn[6-1]=i6;
            elmtconn[8-1]=i8;
            elmtconn[9-1]=i9;
            elmtconn[10-1]=i10;
            }
            return;
        case 12:
            // 27-node second order hexahedron
            return;
        case 13:
            // 18-node second order prism
            return;
        case 14:
            // 14-node second order pyramid
            return;
        case 15:
            // 1-node point
            return;
        case 16:
            // 8-node second order quadrangle
            return;
        case 17:
            // 20-node second order hexahedron
            return;
        case 18:
            // 15-node second order prism
            return;
        case 19:
            // 13-node second order pyramid
            return;
        case 20:
            // 9-node third order incomplete triangle
            return;
        case 21:
            // 10-node third order triangle
            return;
        case 22:
            // 12-node fourth order incomplete triangle
            return;
        case 23:
            // 15-node fourth order triangle
            return;
        case 24:
            // 15-node fifth order incomplete triangle
            return;
        case 25:
            // 21-node fifth order complete triangle
            return;
        case 26:{
            // 4-node third order edge
            // 0---2---3---1 to 0---1---2---3
            int i2=elmtconn[4-1];
            int i3=elmtconn[2-1];
            int i4=elmtconn[3-1];
            elmtconn[2-1]=i2;
            elmtconn[3-1]=i3;
            elmtconn[4-1]=i4;
            }
            return;
        case 27:
            // 5-node fourth order edge
            return;
        case 28:
            // 6-node fifth order edge
            return;
        case 29:
            // 20-node third order tetrahedron
            return;
        case 30:
            // 35-node fourth order tetrahedron
            return;
        case 31:
            // 56-node fifth order tetrahedron
            return;
        case 92:
            // 64-node third order hexahedron
            return;
        case 93:
            // 125-node fourth order hexahedron
            return;
        default:
            MessagePrinter::printErrorTxt("Unsupported gmsh2 element type in reorderNodesIndex, please check your code or your mesh-importer");
            MessagePrinter::exitAsFem();
            return;
    }
    
}