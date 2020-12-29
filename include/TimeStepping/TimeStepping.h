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

#pragma once

#include <iostream>
#include <string>

#include "Utils/MessagePrinter.h"

#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"

#include "BCSystem/BCSystem.h"
#include "ICSystem/ICSystem.h"

#include "MateSystem/MateSystem.h"
#include "ElmtSystem/ElmtSystem.h"

#include "FE/FE.h"
#include "FESystem/FESystem.h"

#include "SolutionSystem/SolutionSystem.h"
#include "EquationSystem/EquationSystem.h"

#include "OutputSystem/OutputSystem.h"

#include "TimeStepping/TimeSteppingBlock.h"
#include "TimeStepping/TimeSteppingType.h"

using namespace std;


class TimeStepping{
public:
    TimeStepping();

    void InitTimeStepping(TimeSteppingBlock &timeSteppingBlock);

private:
    //*****************************************************************
    //*** basic variables for time stepping
    //*****************************************************************
    double _Dt=1.0e-5;
    double _FinalT=1.0e-3;
    bool _Adaptive=false;
    long int _TotalSteps=-1;
    TimeSteppingType _TimeSteppingType=TimeSteppingType::BACKWARDEULER;
    string _TimeSteppingTypeName="backward-euler";

};