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
//+++ Date   : 2022.05.09
//+++ Purpose: the FEProblem class of AsFem, the top level of the
//+++          whole program
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "InputSystem/InputSystem.h"
#include "FECell/FECell.h"
#include "DofHandler/DofHandler.h"
#include "ElmtSystem/ElmtSystem.h"
#include "MateSystem/MateSystem.h"
#include "FE/FE.h"
#include "FESystem/FESystem.h"
#include "BCSystem/BCSystem.h"
#include "ICSystem/ICSystem.h"
#include "ProjectionSystem/ProjectionSystem.h"
#include "EquationSystem/EquationSystem.h"
#include "SolutionSystem/SolutionSystem.h"
#include "LinearSolver/LinearSolver.h"
#include "NonlinearSolver/NonlinearSolver.h"
#include "TimeStepping/TimeStepping.h"
#include "OutputSystem/OutputSystem.h"
#include "Postprocess/Postprocessor.h"
#include "FEProblem/FEJobBlock.h"
#include "FEProblem/FEControlInfo.h"

#include "Utils/Timer.h"
#include "Utils/MessagePrinter.h"


/**
 * the class implement the FEProblem for the whole procedure of FEM calculation
 */
class FEProblem{
public:
    /**
     * constructor
     */
    FEProblem();

    /**
     * initialize the fe problem,this will read the input file, allocate the memory, but do not
     * execute the 'real' FEM simulation.
     */
    void initFEProblem(int args,char *argv[]);

    /**
     * execute the FEM analysis
     */
    void run();

    /**
     * finalize the FEM simulation
     */
    void finalize();

private:
    /**
     * run the FEM simulation for static problem
     */
    void runStaticAnalysis();
    /**
     * run the FEM simulation for transient problem
     */
    void runTransientAnalysis();
    
private:
    InputSystem m_InputSystem;/**< input system for input file reading */
    FECell m_FECell;/**< fe cell class for mesh I/O, mesh generation, etc. */
    DofHandler m_DofHandler;/**< dof class for Dof map generation and management */
    ElmtSystem m_ElmtSystem;/**< the element system class */
    MateSystem m_MateSystem;/**< the material system class */
    FE m_FE;/**< the fe space class */
    FESystem m_FESystem;/**< the FE system class */
    Timer m_Timer;/**< timer class for time counting */
    BCSystem m_BCSystem;/**< boundary condition system */
    ICSystem m_ICSystem;/**< initial condition system */
    ProjectionSystem m_ProjSystem;/**< projection system */
    EquationSystem m_EqSystem;/**< equation system */
    SolutionSystem m_SolnSystem;/**< solution system */
    LinearSolver m_LinearSolver;/**< the linear solver system */
    NonlinearSolver m_NLSolver;/**< the nonlinear solver system */
    TimeStepping m_TimeStepping;/**< the time stepping system */
    OutputSystem m_Output;/**< the output system */
    Postprocessor m_PostProcessor;/**< the post processor system */

    FEJobBlock m_JobBlock;/**< the fe job block */
    FEControlInfo m_FECtrlInfo;/**< the fe control block */

};