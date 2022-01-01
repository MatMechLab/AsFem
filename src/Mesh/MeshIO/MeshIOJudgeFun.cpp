//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.06.26
//+++ Purpose: implement several judge functions to detect which 
//+++          kind of mesh file we are handling
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Mesh/MeshIO.h"

bool MeshIO::IsGmsh2MeshFile(string filename){
    _in.open(filename.c_str(),ios::in);
    if(!_in.is_open()){
        string str="can\'t read the .msh file(="+filename+"), please make sure your mesh file name is correct";
        MessagePrinter::PrintErrorTxt(str);
        MessagePrinter::AsFem_Exit();
        return false;
    }
    double version;
    bool HasVersion=false;
    string str;
    int format,size;
    version=0.0;HasVersion=false;
    while(!_in.eof()){
        getline(_in,str);
        if(str.find("$MeshFormat")!=string::npos){
            _in>>version>>format>>size;
            HasVersion=true;
            if((version!=2.0)&&(version!=2.1)&&(version!=2.2)){
                _in.close();
                return false;
            }
        }
    }
    _in.close();
    if(!HasVersion){
        return false;
    }
    else{
        return true;
    }
}
//***********************************************
bool MeshIO::IsGmsh4MeshFile(string filename){
    _in.open(filename.c_str(),ios::in);
    if(!_in.is_open()){
        string str="can\'t read the .msh file(="+filename+"), please make sure your mesh file name is correct";
        MessagePrinter::PrintErrorTxt(str);
        MessagePrinter::AsFem_Exit();
        return false;
    }
    double version;
    bool HasVersion=false;
    string str;
    int format,size;
    version=0.0;HasVersion=false;
    while(!_in.eof()){
        getline(_in,str);
        if(str.find("$MeshFormat")!=string::npos){
            _in>>version>>format>>size;
            HasVersion=true;
            if((version!=4.0)&&(version!=4.1)&&(version!=4.2)){
                _in.close();
                return false;
            }
        }
    }
    _in.close();
    if(!HasVersion){
        return false;
    }
    else{
        return true;
    }
}
bool MeshIO::IsNetgenMeshFile(string filename){
    // the extension of netgen mesh file is: xxx.gmsh2
    int iInd;
    string str;
    iInd=filename.find_last_of(".");
    str=filename.substr(iInd+1);
    if(str.find("gmsh2")!=string::npos){
        return true;
    }
    else{
        return false;
    }
}
bool MeshIO::IsAbaqusMeshFile(string filename){
    // the extension of abaqus mesh file is: xxx.inp
    int iInd;
    string str;
    iInd=filename.find_last_of(".");
    str=filename.substr(iInd+1);
    if(str.find("inp")!=string::npos){
        return true;
    }
    else{
        return false;
    }
}