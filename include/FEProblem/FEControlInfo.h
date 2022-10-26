//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
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

#include "TimeStepping/TimeSteppingType.h"
#include "petsc.h"


/**
 * This class defines the basic information to control to FE calculation
 */
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

        m_timesteppingtype=TimeSteppingType::BACKWARDEULER;
    }

    void init(){
        ctan[0]=1.0;ctan[1]=0.0;ctan[2]=0.0;
        dt=1.0e-6;
        t=1.0;
        CurrentStep=0;
        FinalStep=0;
        IsDebug=true;
        IsDepDebug=false;
        IsProjection=false;

        m_timesteppingtype=TimeSteppingType::BACKWARDEULER;
    }

    double ctan[3];
    double dt=1.0e-6;
    double t=1.0;
    int CurrentStep=0;
    int FinalStep=0;
    bool IsDebug=true;
    bool IsDepDebug=false;
    bool IsProjection=false;

    // for time stepping
    TimeSteppingType m_timesteppingtype=TimeSteppingType::BACKWARDEULER;

};