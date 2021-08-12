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
//+++ Date   : 2021.08.12
//+++ Purpose: Implement the materials class for AsFem, this class
//+++          contains the scalar, vector, rank-2, and rank4 type
//+++          materials
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/Materials.h"

Materials::Materials(){
    _ScalarMaterials.clear();
    _VectorMaterials.clear();
    _Rank2Materials.clear();
    _Rank4Materials.clear();
    _VectorNull.setZero();
    _Rank2Null.SetToZeros();
    _Rank4Null.SetToZeros();
}
Materials::Materials(const Materials &newmate){
    Clean();
    for(const auto &it:newmate._ScalarMaterials){
        _ScalarMaterials[it.first]=it.second;
    }
    for(const auto &it:newmate._VectorMaterials){
        _VectorMaterials[it.first]=it.second;
    }
    for(const auto &it:newmate._Rank2Materials){
        _Rank2Materials[it.first]=it.second;
    }
    for(const auto &it:newmate._Rank4Materials){
        _Rank4Materials[it.first]=it.second;
    }
}
//*************************************************
void Materials::Clean(){
    _ScalarMaterials.clear();
    _VectorMaterials.clear();
    _Rank2Materials.clear();
    _Rank4Materials.clear();
}
//*************************************************
double& Materials::ScalarMaterials(string matename){
    return _ScalarMaterials[matename];
}
double Materials::ScalarMaterials(string matename)const{
    try {
        return _ScalarMaterials.at(matename);
    } catch (...) {
        MessagePrinter::PrintErrorTxt("cant find "+matename+" in your scalar materials, please check your material definiation in the material system");
        MessagePrinter::AsFem_Exit();
        return -1;
    }
}
//************************************************
Vector3d& Materials::VectorMaterials(string matename){
    return _VectorMaterials[matename];
}
Vector3d Materials::VectorMaterials(string matename)const{
    try {
        return _VectorMaterials.at(matename);
    } catch (...) {
        MessagePrinter::PrintErrorTxt("cant find "+matename+" in your vector materials, please check your material definiation in the material system");
        MessagePrinter::AsFem_Exit();
        return _VectorNull;
    }
}
//*************************************************
RankTwoTensor& Materials::Rank2Materials(string matename){
    return _Rank2Materials[matename];
}
RankTwoTensor Materials::Rank2Materials(string matename)const{
    try {
        return _Rank2Materials.at(matename);
    } catch (...) {
        MessagePrinter::PrintErrorTxt("cant find "+matename+" in your rank-2 tensor materials, please check your material definiation in the material system");
        MessagePrinter::AsFem_Exit();
        return _Rank2Null;
    }
}
//**************************************************
RankFourTensor& Materials::Rank4Materials(string matename){
    return _Rank4Materials[matename];
}
RankFourTensor Materials::Rank4Materials(string matename)const{
    try {
        return _Rank4Materials.at(matename);
    } catch (...) {
        MessagePrinter::PrintErrorTxt("cant find "+matename+" in your rank-4 tensor materials, please check your material definiation in the material system");
        MessagePrinter::AsFem_Exit();
        return _Rank4Null;
    }
}
