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
//+++ Date   : 2022.05.09
//+++ Purpose: the FEProblem class of AsFem, the top level of the
//+++          whole program
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "InputSystem/InputSystem.h"
#include "Mesh/Mesh.h"
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
    InputSystem m_inputSystem;/**< input system for input file reading */
    Mesh m_mesh;/**< mesh class for mesh I/O, mesh generation, etc. */
    DofHandler m_dofhandler;/**< dof class for Dof map generation and management */
    ElmtSystem m_elmtsystem;/**< the element system class */
    MateSystem m_matesystem;/**< the material system class */
    FE m_fe;/**< the fe space class */
    FESystem m_fesystem;/**< the FE system class */
    Timer m_timer;/**< timer class for time counting */
    BCSystem m_bcsystem;/**< boundary condition system */
    ICSystem m_icsystem;/**< initial condition system */
    ProjectionSystem m_projsystem;/**< projection system */
    EquationSystem m_equationsystem;/**< equation system */
    SolutionSystem m_solutionsystem;/**< solution system */
    NonlinearSolver m_nlsolver;/**< the nonlinear solver system */
    TimeStepping m_timestepping;/**< the time stepping system */
    OutputSystem m_output;/**< the output system */
    Postprocessor m_postprocessor;/**< the post processor system */

    FEJobBlock m_jobblock;/**< the fe job block */
    FEControlInfo m_fectrlinfo;/**< the fe control block */

};