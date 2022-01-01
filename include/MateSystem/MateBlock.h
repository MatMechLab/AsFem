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
//+++ Date   : 2020.07.10
//+++ Purpose: Define the [mate] block
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once


#include <iostream>
#include <iomanip>

#include "Utils/MessagePrinter.h"

#include "MateSystem/MateType.h"

using namespace std;

/**
 * This class defines the basic information for [mate] block 
 * in your input file
 */
class MateBlock{
public:
    MateBlock(){
        _MateBlockName.clear();
        _MateTypeName.clear();
        _MateType=MateType::NULLMATE;
        _Parameters.clear();
        _MateBlockIndex=0;
    }

    string         _MateBlockName;
    string         _MateTypeName;
    MateType       _MateType;
    vector<double> _Parameters;
    int            _MateBlockIndex;

    void Init(){
        _MateBlockName.clear();
        _MateTypeName.clear();
        _MateType=MateType::NULLMATE;
        _Parameters.clear();
        _MateBlockIndex=0;
    }
};
