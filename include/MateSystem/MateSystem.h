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
//+++ Purpose: Implement the materials system for AsFem
//+++          this class offer some built-in material models
//+++          as well as the User-Defined-Material (umat) models
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <iomanip>
#include <cmath>


//**********************************
//*** For AsFem's own header file
//**********************************
#include "Utils/MessagePrinter.h"

#include "MateSystem/MateBlock.h"

class MateSystem{
public:
    MateSystem();

    void AddMateBlock2List(MateBlock &mateblock);
    //********************************************
    //*** for some basic getting functions
    //********************************************
    inline int GetMateBlockNums()const{return _nMateBlocks;}



private:
    int _nMateBlocks;
    vector<MateBlock> _MateBlockList;
};