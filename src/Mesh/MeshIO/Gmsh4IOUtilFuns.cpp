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
//+++ Date   : 2021.01.24
//+++ Purpose: add ifstream utils for msh4 file reading
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Mesh/Gmsh4IO.h"


bool Gmsh4IO::NextLineIsSingleNumber(ifstream &in){
    string str;
    vector<double> numbers;

    streampos oldpos=in.tellg();
    getline(in,str);
    numbers=StringUtils::SplitStrNum(str);
    in.seekg(oldpos);
    if(numbers.size()==1){
        return true;
    }
    return false;
}
bool Gmsh4IO::NextLineIsEntities(ifstream &in){
    string str;
    vector<double> numbers;
    streampos oldpos=in.tellg();
    getline(in,str);
    numbers=StringUtils::SplitStrNum(str);
    in.seekg(oldpos);
    if(numbers.size()==4){
        return true;
    }
    return false;
}
bool Gmsh4IO::NextLineIsCoordinate(ifstream &in){
    string str;
    vector<double> numbers;
    streampos oldpos=in.tellg();
    getline(in,str);
    numbers=StringUtils::SplitStrNum(str);
    in.seekg(oldpos);
    if(numbers.size()==4){
        return true;
    }
    return false;
}
