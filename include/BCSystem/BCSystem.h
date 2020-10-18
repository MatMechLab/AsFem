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
//+++ Date   : 2020.07.10
//+++ Purpose: Define the boundary condition system in AsFem
//+++          Here we can apply:
//+++               1) dirichlet bc, i.e. displacement, temperature ...
//+++               2) neuman bc, i.e. flux, force
//+++               3) robin bc as well as user-defined-bc (ubc)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <string>
#include <vector>


//******************************************
//*** for AsFem own header
//******************************************
#include "Utils/MessagePrinter.h"

#include "BCSystem/BCBlock.h"
#include "BCSystem/BCType.h"

using namespace std;

class BCSystem{
public:
    BCSystem();
    void AddBCBlock2List(BCBlock &bcblock);

    inline int GetBCBlockNums()const{return _nBCBlocks;}
    inline BCBlock GetIthBCBlock(const int &i)const{return _BCBlockList[i-1];}

private:
    int _nBCBlocks;
    vector<BCBlock> _BCBlockList;
};