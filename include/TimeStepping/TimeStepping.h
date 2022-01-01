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
#include "Postprocess/Postprocess.h"

#include "NonlinearSolver/NonlinearSolver.h"

#include "TimeStepping/TimeSteppingBlock.h"
#include "TimeStepping/TimeSteppingType.h"

#include "FEProblem/FEControlInfo.h"

using namespace std;


/**
 * This class responsible for the time stepping of transient analysis
 */
class TimeStepping{
public:
    TimeStepping();

    /**
     * Setup the timeStepping solver from the [timestepping] block
     */
    void SetOpitonsFromTimeSteppingBlock(TimeSteppingBlock &timeSteppingBlock);

    /**
     * Get the current time stepping method
     */
    TimeSteppingType GetCurrentSteppingMethod()const{return _TimeSteppingType;}

    /**
     * Check wether the adaptive is enabled
     */
    bool IsAdaptive()const{return _Adaptive;}

    /**
     * Do the time stepping until the maximum step is arrived
     */
    bool Solve(Mesh &mesh,DofHandler &dofHandler,
            ElmtSystem &elmtSystem,MateSystem &mateSystem,
            BCSystem &bcSystem,ICSystem &icSystem,
            SolutionSystem &solutionSystem,EquationSystem &equationSystem,
            FE &fe,FESystem &feSystem,
            OutputSystem &outputSystem,
            Postprocess &postprocessSystem,
            FEControlInfo &fectrlinfo,
            NonlinearSolver &nonlinearSolver);

    /**
     * Print out the basic information of time stepping class to your terminal
     */
    void PrintTimeSteppingInfo()const;
private:
    //*****************************************************************
    //*** basic variables for time stepping
    //*****************************************************************
    double _Dt=1.0e-5;
    double _FinalT=1.0e-3;
    bool _Adaptive=false;
    long int _TotalSteps=-1;
    long int _CurrentStep=-1;
    TimeSteppingType _TimeSteppingType=TimeSteppingType::BACKWARDEULER;
    string _TimeSteppingTypeName="backward-euler";
    double _GrowthFactor=1.1,_CutBackFactor=0.85;
    int _OptIters;
    double _DtMin,_DtMax;

};
