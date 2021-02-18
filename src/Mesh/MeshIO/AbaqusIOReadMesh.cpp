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

bool AbaqusIO::ReadMeshFromFile(Mesh &mesh){
    if(mesh.GetDim()){
        return false;
    }
    else{
        return true;
    }
}