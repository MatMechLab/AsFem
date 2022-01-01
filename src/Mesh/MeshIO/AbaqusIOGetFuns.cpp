//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author  : Yang Bai
//+++ Date    : 2021.02.18
//+++ Reviewer: Xiaoyuan @2021.08.28
//+++ Purpose : add getting functions for abaqus io importor
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Mesh/AbaqusIO.h"

//********************************************************
int AbaqusIO::GetElmtNodesNumFromInp() const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get element nodes number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int i;
    substr.clear();
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Element, type=")!=string::npos){
            i=str.find("=");
            substr=str.substr(i+1,string::npos);
            break;
        }
    }
    in.close();
    if(substr.find("CPS4")!=string::npos){
        return 4;
    }
    else if(substr.find("CPS8")!=string::npos){
        return 8;
    }
    else if(substr.find("CPS3")!=string::npos){
        return 3;
    }
    else if(substr.find("CPS6")!=string::npos){
        return 6;
    }
    else if(substr.find("C3D4")!=string::npos){
        return 4;
    }
    else if(substr.find("C3D8")!=string::npos){
        return 8;
    }
    else if(substr.find("C3D10")!=string::npos){
        return 10;
    }
    else if(substr.find("C3D20")!=string::npos){
        return 20;
    }
    else if(substr.find("C3D27")!=string::npos){
        return 27;
    }
    else{
        MessagePrinter::PrintErrorTxt("can not call GetElmtNodesNumFromInp, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return -1;
}
//******************************************************************
int AbaqusIO::GetElmtDimFromInp() const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get element dimension number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int i;
    substr.clear();
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Element, type=")!=string::npos){
            i=str.find("=");
            substr=str.substr(i+1,string::npos);
            break;
        }
    }
    in.close();
    if(substr.find("CPS4")!=string::npos){
        return 2;
    }
    else if(substr.find("CPS8")!=string::npos){
        return 2;
    }
    else if(substr.find("CPS3")!=string::npos){
        return 2;
    }
    else if(substr.find("CPS6")!=string::npos){
        return 2;
    }
    else if(substr.find("C3D4")!=string::npos){
        return 3;
    }
    else if(substr.find("C3D8")!=string::npos){
        return 3;
    }
    else if(substr.find("C3D10")!=string::npos){
        return 3;
    }
    else if(substr.find("C3D20")!=string::npos){
        return 3;
    }
    else if(substr.find("C3D27")!=string::npos){
        return 3;
    }
    else{
        MessagePrinter::PrintErrorTxt("can not call GetElmtDimFromInp, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return -1;
}
//***************************************************************
int AbaqusIO::GetElmtOrderFromInp()const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get element dimension number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int i;
    substr.clear();
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Element, type=")!=string::npos){
            i=str.find("=");
            substr=str.substr(i+1,string::npos);
            break;
        }
    }
    in.close();
    if(substr.find("CPS4")!=string::npos){
        return 1;
    }
    else if(substr.find("CPS8")!=string::npos){
        return 2;
    }
    else if(substr.find("CPS3")!=string::npos){
        return 1;
    }
    else if(substr.find("CPS6")!=string::npos){
        return 2;
    }
    else if(substr.find("C3D4")!=string::npos){
        return 1;
    }
    else if(substr.find("C3D8")!=string::npos){
        return 1;
    }
    else if(substr.find("C3D10")!=string::npos){
        return 2;
    }
    else if(substr.find("C3D20")!=string::npos){
        return 2;
    }
    else if(substr.find("C3D27")!=string::npos){
        return 2;
    }
    else{
        MessagePrinter::PrintErrorTxt("can not call GetElmtOrderFromInp, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return -1;
}
//*****************************************************
int AbaqusIO::GetElmtVTKCellTypeFromInp()const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get element dimension number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int i;
    substr.clear();
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Element, type=")!=string::npos){
            i=str.find("=");
            substr=str.substr(i+1,string::npos);
            break;
        }
    }
    in.close();
    if(substr.find("CPS4")!=string::npos){
        // 4-node quadrangle
        return 9;
    }
    else if(substr.find("CPS8")!=string::npos){
        // 8-node second order quadrangle
        return 23;
    }
    else if(substr.find("CPS3")!=string::npos){
        // 3-node triangle
        return 5;
    }
    else if(substr.find("CPS6")!=string::npos){
        // 6-node second order triangle
        return 22;
    }
    else if(substr.find("C3D4")!=string::npos){
        // 4-node tetrahedron
        return 10;
    }
    else if(substr.find("C3D8")!=string::npos){
        // 8-node hexahedron
        return 12;
    }
    else if(substr.find("C3D10")!=string::npos){
        // 10-node second order tetrahedron
        return 24;
    }
    else if(substr.find("C3D20")!=string::npos){
        // 20-node second order hexahedron
        return 25;
    }
    else if(substr.find("C3D27")!=string::npos){
        // 27-node second order hexahedron
        return 29;
    }
    else{
        MessagePrinter::PrintErrorTxt("can not call GetElmtVTKCellTypeFromInp, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return -1;
}
//**************************************************
MeshType AbaqusIO::GetElmtMeshTypeFromInp()const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get element dimension number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int i;
    substr.clear();
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Element, type=")!=string::npos){
            i=str.find("=");
            substr=str.substr(i+1,string::npos);
            break;
        }
    }
    in.close();
    if(substr.find("CPS4")!=string::npos){
        // 4-node quadrangle
        return MeshType::QUAD4;
    }
    else if(substr.find("CPS8")!=string::npos){
        // 8-node second order quadrangle
        return MeshType::QUAD8;
    }
    else if(substr.find("CPS3")!=string::npos){
        // 3-node triangle
        return MeshType::TRI3;
    }
    else if(substr.find("CPS6")!=string::npos){
        // 6-node second order triangle
        return MeshType::TRI6;
    }
    else if(substr.find("C3D4")!=string::npos){
        // 4-node tetrahedron
        return MeshType::TET4;
    }
    else if(substr.find("C3D8")!=string::npos){
        // 8-node hexahedron
        return MeshType::HEX8;
    }
    else if(substr.find("C3D10")!=string::npos){
        // 10-node second order tetrahedron
        return MeshType::TET10;
    }
    else if(substr.find("C3D20")!=string::npos){
        // 20-node second order hexahedron
        return MeshType::HEX20;
    }
    else if(substr.find("C3D27")!=string::npos){
        // 27-node second order hexahedron
        return MeshType::HEX27;
    }
    else{
        MessagePrinter::PrintErrorTxt("can not call GetElmtMeshTypeFromInp, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return MeshType::NULLTYPE;
}
//************************************************************
string AbaqusIO::GetElmtMeshTypeNameFromInp()const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get element dimension number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int i;
    substr.clear();
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Element, type=")!=string::npos){
            i=str.find("=");
            substr=str.substr(i+1,string::npos);
            break;
        }
    }
    in.close();
    if(substr.find("CPS4")!=string::npos){
        // 4-node quadrangle
        return "quad4";
    }
    else if(substr.find("CPS8")!=string::npos){
        // 8-node second order quadrangle
        return "quad8";
    }
    else if(substr.find("CPS3")!=string::npos){
        // 3-node triangle
        return "tri3";
    }
    else if(substr.find("CPS6")!=string::npos){
        // 6-node second order triangle
        return "tri6";
    }
    else if(substr.find("C3D4")!=string::npos){
        // 4-node tetrahedron
        return "tet4";
    }
    else if(substr.find("C3D8")!=string::npos){
        // 8-node hexahedron
        return "hex8";
    }
    else if(substr.find("C3D10")!=string::npos){
        // 10-node second order tetrahedron
        return "tet10";
    }
    else if(substr.find("C3D20")!=string::npos){
        // 20-node second order hexahedron
        return "hex20";
    }
    else if(substr.find("C3D27")!=string::npos){
        // 27-node second order hexahedron
        return "hex27";
    }
    else{
        MessagePrinter::PrintErrorTxt("can not call GetElmtMeshTypeFromInp, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return "unknowm-mesh";
}
//************************************************************
int AbaqusIO::GetSubElmtNodesNumFromInp() const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get sub element nodes number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int i;
    substr.clear();
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Element, type=")!=string::npos){
            i=str.find("=");
            substr=str.substr(i+1,string::npos);
            break;
        }
    }
    in.close();
    if(substr.find("CPS4")!=string::npos){
        return 2;
    }
    else if(substr.find("CPS8")!=string::npos){
        return 3;
    }
    else if(substr.find("CPS3")!=string::npos){
        return 2;
    }
    else if(substr.find("CPS6")!=string::npos){
        return 3;
    }
    else if(substr.find("C3D4")!=string::npos){
        return 3;
    }
    else if(substr.find("C3D8")!=string::npos){
        return 4;
    }
    else if(substr.find("C3D10")!=string::npos){
        // 10-node second order tetrahedron
        return 6;
    }
    else if(substr.find("C3D20")!=string::npos){
        return 8;
    }
    else if(substr.find("C3D27")!=string::npos){
        return 9;
    }
    else{
        MessagePrinter::PrintErrorTxt("can not call GetSubElmtNodesNumFromInp, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return -1;
}
//*****************************************************************
int AbaqusIO::GetSubSubElmtNodesNumFromInp() const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get sub sub element nodes number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int i;
    substr.clear();
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Element, type=")!=string::npos){
            i=str.find("=");
            substr=str.substr(i+1,string::npos);
            break;
        }
    }
    in.close();
    if(substr.find("CPS4")!=string::npos){
        // quad4-edge2-point
        return 1;
    }
    else if(substr.find("CPS8")!=string::npos){
        // quad8-edge3-point
        return 1;
    }
    else if(substr.find("CPS3")!=string::npos){
        // tri3-edge2-point
        return 1;
    }
    else if(substr.find("CPS6")!=string::npos){
        // tri6-edge3-point
        return 1;
    }
    else if(substr.find("C3D4")!=string::npos){
        // tet4-tri3-edge2
        return 2;
    }
    else if(substr.find("C3D8")!=string::npos){
        // hex8-quad4-edge2
        return 2;
    }
    else if(substr.find("C3D10")!=string::npos){
        // 10-node second order tetrahedron-tri6-edge3
        return 3;
    }
    else if(substr.find("C3D20")!=string::npos){
        // hex20-quad8-edge3
        return 3;
    }
    else if(substr.find("C3D27")!=string::npos){
        // hex27-quad9-edge3
        return 3;
    }
    else{
        MessagePrinter::PrintErrorTxt("can not call GetSubSubElmtNodesNumFromInp, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return -1;
}
//*****************************************************************
int AbaqusIO::GetSubElmtDimFromInp() const{
    return GetElmtDimFromInp()-1;
}
int AbaqusIO::GetSubElmtOrderFromInpElmt()const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get element order number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int i;
    substr.clear();
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Element, type=")!=string::npos){
            i=str.find("=");
            substr=str.substr(i+1,string::npos);
            break;
        }
    }
    in.close();
    if(substr.find("CPS4")!=string::npos){
        return 1;
    }
    else if(substr.find("CPS8")!=string::npos){
        return 2;
    }
    else if(substr.find("CPS3")!=string::npos){
        return 1;
    }
    else if(substr.find("CPS6")!=string::npos){
        return 2;
    }
    else if(substr.find("C3D4")!=string::npos){
        return 1;
    }
    else if(substr.find("C3D8")!=string::npos){
        return 1;
    }
    else if(substr.find("C3D10")!=string::npos){
        // 10-node second order tetrahedron
        return 2;
    }
    else if(substr.find("C3D20")!=string::npos){
        return 2;
    }
    else if(substr.find("C3D27")!=string::npos){
        return 2;
    }
    else{
        MessagePrinter::PrintErrorTxt("can not call GetSubElmtOrderFromInpElmt, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return -1;
}
//*****************************************************************
MeshType AbaqusIO::GetSubElmtMeshTypeFromInp()const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get sub element type information, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int i;
    substr.clear();
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Element, type=")!=string::npos){
            i=str.find("=");
            substr=str.substr(i+1,string::npos);
            break;
        }
    }
    in.close();
    if(substr.find("CPS4")!=string::npos){
        // 4-node quadrangle
        return MeshType::EDGE2;
    }
    else if(substr.find("CPS8")!=string::npos){
        // 8-node second order quadrangle
        return MeshType::EDGE3;
    }
    else if(substr.find("CPS3")!=string::npos){
        // 3-node triangle
        return MeshType::EDGE2;
    }
    else if(substr.find("CPS6")!=string::npos){
        // 6-node second order triangle
        return MeshType::EDGE3;
    }
    else if(substr.find("C3D4")!=string::npos){
        // 4-node tetrahedron
        return MeshType::TRI3;
    }
    else if(substr.find("C3D8")!=string::npos){
        // 8-node hexahedron
        return MeshType::QUAD4;
    }
    else if(substr.find("C3D10")!=string::npos){
        // 10-node second order tetrahedron
        return MeshType::TRI6;
    }
    else if(substr.find("C3D20")!=string::npos){
        // 20-node second order hexahedron
        return MeshType::QUAD8;
    }
    else if(substr.find("C3D27")!=string::npos){
        // 27-node second order hexahedron
        return MeshType::QUAD9;
    }
    else{
        MessagePrinter::PrintErrorTxt("can not call GetSubElmtMeshTypeFromInp, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return MeshType::NULLTYPE;
}
//************************************************************
MeshType AbaqusIO::GetSubSubElmtMeshTypeFromInp()const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get sub sub element type information, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int i;
    substr.clear();
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Element, type=")!=string::npos){
            i=str.find("=");
            substr=str.substr(i+1,string::npos);
            break;
        }
    }
    in.close();
    if(substr.find("CPS4")!=string::npos){
        // 4-node quadrangle-egde2-point
        return MeshType::POINT;
    }
    else if(substr.find("CPS8")!=string::npos){
        // 8-node second order quadrangle-edge2-point
        return MeshType::POINT;
    }
    else if(substr.find("CPS3")!=string::npos){
        // 3-node triangle-edge2-point
        return MeshType::POINT;
    }
    else if(substr.find("CPS6")!=string::npos){
        // 6-node second order triangle-edge3-point
        return MeshType::POINT;
    }
    else if(substr.find("C3D4")!=string::npos){
        // 4-node tetrahedron-tri3-edge2
        return MeshType::EDGE2;
    }
    else if(substr.find("C3D8")!=string::npos){
        // 8-node hexahedron-quad4-edge2
        return MeshType::EDGE2;
    }
    else if(substr.find("C3D10")!=string::npos){
        // 10-node second order tetrahedron
        return MeshType::EDGE3;
    }
    else if(substr.find("C3D20")!=string::npos){
        // 20-node second order hexahedron-quad8-edge3
        return MeshType::EDGE3;
    }
    else if(substr.find("C3D27")!=string::npos){
        // 27-node second order hexahedron-quad9-edge3
        return MeshType::EDGE3;
    }
    else{
        MessagePrinter::PrintErrorTxt("can not call GetSubSubElmtMeshTypeFromInp, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return MeshType::NULLTYPE;
}
//************************************************************
int AbaqusIO::GetNodesNumFromInp()const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get nodes number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int nNodes;
    substr.clear();
    nNodes=0;
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Node")!=string::npos){
            nNodes=0;
            getline(in,str);// this line already contains 1st node's x,y,z
            while(str.find("*")==string::npos){
                nNodes+=1;
                getline(in,str);
            }
        }
    }
    in.close();
    return nNodes;
}
//*******************************************************
int AbaqusIO::GetElmtsNumFromInp()const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get elements number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int nElmts;
    substr.clear();
    nElmts=0;
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Element,")!=string::npos){
            nElmts=0;
            getline(in,str);// this line already contains 1st node's x,y,z
            while(str.find("*")==string::npos){
                nElmts+=1;
                getline(in,str);
            }
        }
    }
    in.close();
    return nElmts;
}
//***********************************
int AbaqusIO::GetNsetsNumFromInp()const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get Nset number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int nNsets;
    substr.clear();
    nNsets=0;
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Nset,")!=string::npos){
            nNsets+=1;
        }
    }
    in.close();
    return nNsets;
}
//***********************************
int AbaqusIO::GetElsetsNumFromInp()const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get Elset number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int nElsets;
    substr.clear();
    nElsets=0;
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Elset,")!=string::npos){
            nElsets+=1;
        }
    }
    in.close();
    return nElsets;
}
//***********************************
int AbaqusIO::GetSurfacesNumFromInp()const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get Surface number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,substr;
    int nSurfaces;
    substr.clear();
    nSurfaces=0;
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Surface,")!=string::npos){
            nSurfaces+=1;
        }
    }
    in.close();
    return nSurfaces;
}
//*********************************************************
int AbaqusIO::GetSurfaceElmtsNumFromInp()const{
    ifstream in,insub;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get Surface elements number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,str1,substr;
    string SurfaceSetName;
    vector<double> numbers;
    int nElmts;
    substr.clear();
    nElmts=0;
    int i;
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Surface,")!=string::npos){
            // *Surface, type=ELEMENT, name=Surf-Left
            // _Surf-Left_S4, S4
            // here we only read '_Surf-Left_S4'
            getline(in,str);
            i=str.find(",");
            SurfaceSetName=str.substr(0,i-1);// '_Surf-Left_S4'
            insub.open(_MeshFileName,ios::in);
            while(!insub.eof()){
                getline(insub,str1);
                if(str1.find("*Elset, elset=")!=string::npos&&
                   str1.find(SurfaceSetName)!=string::npos){
                    // now we find the related Elset information
                    if(str1.find("generate")!=string::npos){
                        getline(insub,substr);// read the element index information
                        numbers=StringUtils::SplitStrNum(substr);
                        if(static_cast<int>(numbers.size())!=3){
                            MessagePrinter::PrintErrorTxt(SurfaceSetName+" is using generate way, however, your element index number is not equal to 3");
                            MessagePrinter::AsFem_Exit();
                        }
                        nElmts+=1+(static_cast<int>(numbers[1])-static_cast<int>(numbers[0]))/static_cast<int>(numbers[2]);
                    }
                    else{
                        for(i=0;i<numeric_limits<int>::max();i++){
                            getline(insub,substr);
                            numbers=StringUtils::SplitStrNum(substr);
                            nElmts+=static_cast<int>(numbers.size());
                            if(substr.find("*")!=string::npos){
                                break;
                            }
                        }
                    }
                } // end-of-element-set-information
            }// end-of-while-insub
            insub.close();
        }
    }
    in.close();
    return nElmts;
}
//*********************************************************
int AbaqusIO::GetSurfaceEdgeIDViaSurfaceNameFromInp(string surfacesetname)const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can not get Surface elements number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,str1,substr;
    string SurfaceSetName;
    vector<double> numbers;
    vector<int> ElmtIDs;
    substr.clear();
    ElmtIDs.clear();
    int i,edge;
    edge=0;
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Surface,")!=string::npos && str.find(surfacesetname)!=string::npos){
            // *Surface, type=ELEMENT, name=Surf-Left
            // _Surf-Left_S4, S4
            getline(in,str);// read _Surf-Left_S4, S4
            i=str.find(",");
            substr=str.substr(i+1,string::npos);// '_Surf-Left_S4'
            numbers=StringUtils::SplitStrNum(substr);// split s4 into 4
            if(numbers.size()<1){
                MessagePrinter::PrintErrorTxt("can not find any numbers in "+str+" ,please check your inp file");
                MessagePrinter::AsFem_Exit();
            }
            edge=static_cast<int>(numbers[0]);
        }
    }
    in.close();
    return edge;
}
//*********************************************************
vector<int> AbaqusIO::GetSurfaceElmtIDViaSurfaceNameFromInp(string surfacesetname)const{
    ifstream in,insub;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can nott get Surface elements number, we cant open "+_MeshFileName);
        MessagePrinter::AsFem_Exit();
    }
    string str,str1,substr;
    string SurfaceSetName;
    vector<double> numbers;
    vector<int> ElmtIDs;
    substr.clear();
    int i;
    ElmtIDs.clear();
    while(!in.eof()){
        getline(in,str);
        if(str.find("*Surface,")!=string::npos && str.find(surfacesetname)!=string::npos){
            // *Surface, type=ELEMENT, name=Surf-Left
            // _Surf-Left_S4, S4
            getline(in,str);// read _Surf-Left_S4, S4
            i=str.find(",");
            SurfaceSetName=str.substr(0,i-1);// '_Surf-Left_S4'
            insub.open(_MeshFileName,ios::in);
            while(!insub.eof()){
                getline(insub,str1);
                if(str1.find("*Elset, elset=")!=string::npos&&
                   str1.find(SurfaceSetName)!=string::npos){
                    // now we find the related Elset information
                    if(str1.find("generate")!=string::npos){
                        getline(insub,substr);// read the element index information
                        numbers=StringUtils::SplitStrNum(substr);
                        if(static_cast<int>(numbers.size())!=3){
                            MessagePrinter::PrintErrorTxt(SurfaceSetName+" is using generate way, however, your element index number is not equal to 3");
                            MessagePrinter::AsFem_Exit();
                        }
                        for(int e=static_cast<int>(numbers[0]);e<=static_cast<int>(numbers[1]);e+=static_cast<int>(numbers[2])){
                            ElmtIDs.push_back(e);
                        }
                    }
                    else{
                        for(i=0;i<numeric_limits<int>::max();i++){
                            getline(insub,substr);
                            numbers=StringUtils::SplitStrNum(substr);
                            for(i=0;i<static_cast<int>(numbers.size());i++){
                                ElmtIDs.push_back(static_cast<int>(numbers[i]));
                            }
                            if(substr.find("*")!=string::npos){
                                break;
                            }
                        }
                    }
                } // end-of-element-set-information
            }// end-of-while-insub
            insub.close();
        }
    }
    in.close();
    return ElmtIDs;// please keep in mind, this element id is just the global bulk id, not the boundary element id !!!
}
