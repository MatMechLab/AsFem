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
//+++ Date   : 2020.12.29
//+++ Purpose: Define time stepping system for transient analysis
//+++          in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "FECell/FECell.h"
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
#include "NonlinearSolver/NonlinearSolver.h"
#include "Postprocess/Postprocessor.h"
#include "FEProblem/FEControlInfo.h"

#include "TimeStepping/TimeSteppingType.h"
#include "TimeStepping/TimeSteppingData.h"


/**
 * This class resposible for the time stepping in transient analysis
 */
class TimeStepping{
public:
    /**
     * constructor
     */
    TimeStepping();

    //**************************************************
    //*** general settings
    //**************************************************
    /**
     * set up the time stepping method
     * @param method the time stepping method
     */
    void setTimeSteppingMethod(const TimeSteppingType &method){m_Data.m_SteppingType=method;}
    /**
     * set up the initial delta t
     * @param dt0 the initial delta t
     */
    void setInitialDt(const double &dt0){m_Data.m_Dt0=dt0;}
    /**
     * setup the maximum delta t
     * @param dtmax the maximum delta t
     */
    void setMaxDt(const double &dtmax){m_Data.m_DtMax=dtmax;}
    /**
     * setup the minimum delta t
     * @param dtmin the minimum delta t
     */
    void setMinDt(const double &dtmin){m_Data.m_DtMin=dtmin;}
    /**
     * setup the final time
     * @param t the final t
     */
    void setFinalTime(const double &t){m_Data.m_FinalTime=t;}
    /**
     * setup the optimize nonlinear iterations
     * @param iters the optimize iteration
     */
    void setOptimizeIters(const int &iters){m_Data.m_OptimizeIters=iters;}
    /**
     * setup the cut back factor
     * @param factor the cut back factor value
     */
    void setCutbackFactor(const double &factor){m_Data.m_CutFactor=factor;}
    /**
     * setup the growth factor
     * @param factor the cut back factor value
     */
    void setGrowthFactor(const double &factor){m_Data.m_GrowthFactor=factor;}
    /**
     * setup the adaptive status
     * @param flag true to enable the adaptive time stepping
     */
    void setAdaptiveFlag(const bool &flag){m_Data.m_IsAdaptive=flag;}

    /**
     * apply the default time stepping settings
     */
    void applyDefaultSettings();

    /**
     * print out the time stepping system information
     */
    void printInfo()const;

    /**
     * get the current delta t
     */
    inline double getDt()const{return m_Data.m_Dt;}
    /**
     * get the inital delta t
     */
    inline double getDt0()const{return m_Data.m_Dt0;}
    /**
     * get max delta t
     */
    inline double getMaxDt()const{return m_Data.m_DtMax;}
    /**
     * get min delta t
     */
    inline double getMinDt()const{return m_Data.m_DtMin;}
    /**
     * get final time
     */
    inline double getFinalTime()const{return m_Data.m_FinalTime;}
    /**
     * get the cut back factor
     */
    inline double getCutbackFactor()const{return m_Data.m_CutFactor;}
    /**
     * get the growth factor
     */
    inline double getGrowthFactor()const{return m_Data.m_GrowthFactor;}
    /**
     * get the optimal iterations
     */
    inline int getOptimizeIters()const{return m_Data.m_OptimizeIters;}
    /**
     * get the adaptive status
     */
    inline bool isAdaptive()const{return m_Data.m_IsAdaptive;}
    /**
     * get the time integration method
     */
    inline TimeSteppingType getTimeSteppingType()const{return m_Data.m_SteppingType;}

    /**
     * solve the transient equation
     * @param t_FECell the fe cell class
     * @param t_DofHandler the dof class
     * @param t_FE the fe class
     * @param t_ElmtSystem the element system class
     * @param t_MateSystem the material system class
     * @param t_FESystem the fe system class
     * @param t_BCSystem the boundary condition system
     * @param t_ICSystem the initial condition system
     * @param t_SolnSystem the solution system class
     * @param t_EqSystem the equation system class
     * @param t_ProjSystem the projection system
     * @param t_FECtrlInfo the fe control info
     * @param t_LinearSolver the linear solver system
     * @param t_NLSolver the nonlinear solver
     * @param t_Output  the output system
     * @param t_PostProcess the postprocess system
     */
    bool solve(FECell &t_FECell,
               DofHandler &t_DofHandler,
               FE &t_FE,
               ElmtSystem &t_ElmtSystem,
               MateSystem &t_MateSystem,
               FESystem &t_FESystem,
               BCSystem &t_BCSystem,
               ICSystem &t_ICSystem,
               SolutionSystem &t_SolnSystem,
               EquationSystem &t_EqSystem,
               ProjectionSystem &t_ProjSystem,
               FEControlInfo &t_FECtrlInfo,
               LinearSolver &t_LinearSolver,
               NonlinearSolver &t_NLSolver,
               OutputSystem &t_Output,
               Postprocessor &t_PostProcess);



private:
    TimeSteppingData m_Data;/**< the time stepping data */

};