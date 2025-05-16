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
        Ctan[0]=1.0;Ctan[1]=0.0;Ctan[2]=0.0;
        Dt=1.0e-6;
        T=1.0;
        CurrentStep=0;
        FinalStep=0;
        IsDebug=true;
        IsDepDebug=false;
        IsProjection=false;

        m_TimesteppingType=TimeSteppingType::BACKWARDEULER;
    }

    void init(){
        Ctan[0]=1.0;Ctan[1]=0.0;Ctan[2]=0.0;
        Dt=1.0e-6;
        T=1.0;
        CurrentStep=0;
        FinalStep=0;
        IsDebug=true;
        IsDepDebug=false;
        IsProjection=false;

        m_TimesteppingType=TimeSteppingType::BACKWARDEULER;
    }

    double Ctan[3];
    double Dt=1.0e-6;
    double T=1.0;
    int CurrentStep=0;
    int FinalStep=0;
    bool IsDebug=true;
    bool IsDepDebug=false;
    bool IsProjection=false;

    // for time stepping
    TimeSteppingType m_TimesteppingType=TimeSteppingType::BACKWARDEULER;

};