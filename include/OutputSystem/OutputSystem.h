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
//+++ Date   : 2020.07.12
//+++ Purpose: define the output system for AsFem, where all the 
//+++          results should be written out to the result file
//+++          by this class
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once 

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

#include "OutputSystem/OutputBlock.h"

using namespace std;


class OutputSystem{
public:
    OutputSystem();

    void InitFromOutputBlock(OutputBlock &outputblock);

private:
    int _Interval;
    OutputType _OutputType;
    string _OutputTypeName;
};