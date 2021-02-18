//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.02.18
//+++ Purpose: add getting functions for abaqus io importor
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Mesh/AbaqusIO.h"

//********************************************************
int AbaqusIO::GetElmtNodesNumFromInp() const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can\'t get element nodes number, we can\'t open "+_MeshFileName);
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
    else if(substr.find("C3D20")!=string::npos){
        return 20;
    }
    else if(substr.find("C3D27")!=string::npos){
        return 27;
    }
    else{
        MessagePrinter::PrintErrorTxt("can\'t call GetElmtNodesNumFromInp, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return -1;
}
//******************************************************************
int AbaqusIO::GetElmtDimFromInp() const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can\'t get element dimension number, we can\'t open "+_MeshFileName);
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
    else if(substr.find("C3D20")!=string::npos){
        return 3;
    }
    else if(substr.find("C3D27")!=string::npos){
        return 3;
    }
    else{
        MessagePrinter::PrintErrorTxt("can\'t call GetElmtDimFromInp, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return -1;
}
//***************************************************************
int AbaqusIO::GetElmtOrderFromInp()const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can\'t get element dimension number, we can\'t open "+_MeshFileName);
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
    else if(substr.find("C3D20")!=string::npos){
        return 2;
    }
    else if(substr.find("C3D27")!=string::npos){
        return 2;
    }
    else{
        MessagePrinter::PrintErrorTxt("can\'t call GetElmtOrderFromInp, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return -1;
}
//*****************************************************
int AbaqusIO::GetElmtVTKCellTypeFromInp()const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can\'t get element dimension number, we can\'t open "+_MeshFileName);
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
    else if(substr.find("C3D20")!=string::npos){
        // 20-node second order hexahedron
        return 25;
    }
    else if(substr.find("C3D27")!=string::npos){
        // 27-node second order hexahedron
        return 29;
    }
    else{
        MessagePrinter::PrintErrorTxt("can\'t call GetElmtVTKCellTypeFromInp, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return -1;
}
//**************************************************
MeshType AbaqusIO::GetElmtMeshTypeFromInp()const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can\'t get element dimension number, we can\'t open "+_MeshFileName);
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
    else if(substr.find("C3D20")!=string::npos){
        // 20-node second order hexahedron
        return MeshType::HEX20;
    }
    else if(substr.find("C3D27")!=string::npos){
        // 27-node second order hexahedron
        return MeshType::HEX27;
    }
    else{
        MessagePrinter::PrintErrorTxt("can\'t call GetElmtMeshTypeFromInp, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return MeshType::NULLTYPE;
}

//************************************************************
int AbaqusIO::GetSubElmtNodesNumFromInp() const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can\'t get sub element nodes number, we can\'t open "+_MeshFileName);
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
    else if(substr.find("C3D20")!=string::npos){
        return 8;
    }
    else if(substr.find("C3D27")!=string::npos){
        return 9;
    }
    else{
        MessagePrinter::PrintErrorTxt("can\'t call GetSubElmtNodesNumFromInp, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return -1;
}
int AbaqusIO::GetSubElmtDimFromGmshInp() const{
    return GetElmtDimFromInp()-1;
}
int AbaqusIO::GetSubElmtOrderFromInpElmt()const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can\'t get element order number, we can\'t open "+_MeshFileName);
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
    else if(substr.find("C3D20")!=string::npos){
        return 2;
    }
    else if(substr.find("C3D27")!=string::npos){
        return 2;
    }
    else{
        MessagePrinter::PrintErrorTxt("can\'t call GetSubElmtOrderFromInpElmt, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return -1;
}
MeshType AbaqusIO::GetSubElmtMeshTypeFromInp()const{
    ifstream in;
    in.open(_MeshFileName,ios::in);
    if(!in.is_open()){
        MessagePrinter::PrintErrorTxt("can\'t get sub element type information, we can\'t open "+_MeshFileName);
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
    else if(substr.find("C3D20")!=string::npos){
        // 20-node second order hexahedron
        return MeshType::QUAD8;
    }
    else if(substr.find("C3D27")!=string::npos){
        // 27-node second order hexahedron
        return MeshType::QUAD9;
    }
    else{
        MessagePrinter::PrintErrorTxt("can\'t call GetSubElmtMeshTypeFromInp, unsupported element type");
        MessagePrinter::AsFem_Exit();
    }
    return MeshType::NULLTYPE;
}