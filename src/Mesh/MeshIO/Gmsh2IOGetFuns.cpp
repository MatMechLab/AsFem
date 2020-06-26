//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.06.26
//+++ Purpose: add getting functions for gmsh2 io importor
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Mesh/Gmsh2IO.h"

int Gmsh2IO::GetElmtNodesNumFromGmshElmtType(const int &elmttype) const{
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
            return -1;
    }
}
//******************************************************
int Gmsh2IO::GetElmtDimFromGmshElmtType(const int &elmttype) const{
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
            return -1;
    }
}
//**************************************************
int Gmsh2IO::GetElmtOrderFromGmshElmtType(const int &elmttype)const{
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
            return -1;
    }
}
//***************************************************************
MeshType Gmsh2IO::GetElmtMeshTypeFromGmshElmtType(const int &elmttype)const{
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
            return MeshType::NULLTYPE;
    }
}
string Gmsh2IO::GetElmtTypeNameFromGmshElmtType(const int &elmttype)const{
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
            return "unknown";
    }
}

//****************************************************************
int Gmsh2IO::GetSubElmtNodesNumFromGmshElmtType(const int &elmttype) const{
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
            return -1;
    }
}
//***************************************************************
int Gmsh2IO::GetSubElmtDimFromGmshElmtType(const int &elmttype) const{
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
            return -1;
    }
}
//*************************************************************
int Gmsh2IO::GetSubElmtOrderFromGmshElmtType(const int &elmttype)const{
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
            return -1;
    }
}
//*************************************************************
MeshType Gmsh2IO::GetSubElmtMeshTypeFromGmshElmtType(const int &elmttype)const{
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
            return MeshType::NULLTYPE;
    }
}
