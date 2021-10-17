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
//+++ Date   : 2020.12.26
//+++ Purpose: define some basic information to control our FEM 
//+++          calculation
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>

#include "petsc.h"

using namespace std;

class FEControlInfo{
public:
    FEControlInfo(){
        ctan[0]=1.0;ctan[1]=0.0;ctan[2]=0.0;
        dt=1.0e-6;
        t=1.0;
        CurrentStep=0;
        FinalStep=0;
        IsDebug=true;
        IsDepDebug=false;
        IsProjection=false;
    }

    void Init(){
        ctan[0]=1.0;ctan[1]=0.0;ctan[2]=0.0;
        dt=1.0e-6;
        t=1.0;
        CurrentStep=0;
        FinalStep=0;
        IsDebug=true;
        IsDepDebug=false;
        IsProjection=false;
    }

    double ctan[3];
    double dt=1.0e-6;
    double t=1.0;
    int CurrentStep=0;
    int FinalStep=0;
    bool IsDebug=true;
    bool IsDepDebug=false;
    bool IsProjection=false;

};
