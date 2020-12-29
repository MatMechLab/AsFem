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
//+++ Date   : 2020.12.29
//+++ Purpose: Define time stepping system for transient analysis
//+++          in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "TimeStepping/TimeStepping.h"

TimeStepping::TimeStepping(){
    _Dt=1.0e-5;
    _FinalT=1.0e-3;
    _Adaptive=false;
    _TotalSteps=-1;
    _TimeSteppingType=TimeSteppingType::BACKWARDEULER;
    _TimeSteppingTypeName="backward-euler";
}