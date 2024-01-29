//****************************************************************
//* This file is part of the AsFem framework
//* Advanced Simulation kit based on Finite Element Method (AsFem)
//* All rights reserved, Yang Bai/MM-Lab@CopyRight 2020-present
//* https://github.com/MatMechLab/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.01.17
//+++ Purpose: Implement some commonly used mathematic functions
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MathUtils/MathFuns.h"

MathFuns::MathFuns(){}

double MathFuns::bracketPos(const double &x){
    return 0.5*(x+abs(x));
}
double MathFuns::bracketNeg(const double &x){
    return 0.5*(x-abs(x));
}
//***************************************
double MathFuns::sign(const double &x){
    return x >= 0.0 ? 1.0 : -1.0;
}