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
//+++ Date   : 2020.12.29
//+++ Purpose: Define time stepping system for transient analysis
//+++          in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

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
    void setTimeSteppingMethod(const TimeSteppingType &method){m_data.m_stepping_type=method;}
    /**
     * set up the initial delta t
     * @param dt0 the initial delta t
     */
    void setInitialDt(const double &dt0){m_data.m_dt0=dt0;}
    /**
     * setup the maximum delta t
     * @param dtmax the maximum delta t
     */
    void setMaxDt(const double &dtmax){m_data.m_dtmax=dtmax;}
    /**
     * setup the minimum delta t
     * @param dtmin the minimum delta t
     */
    void setMinDt(const double &dtmin){m_data.m_dtmin=dtmin;}
    /**
     * setup the final time
     * @param t the final t
     */
    void setFinalTime(const double &t){m_data.m_finaltime=t;}
    /**
     * setup the optimize nonlinear iterations
     * @param iters the optimize iteration
     */
    void setOptimizeIters(const int &iters){m_data.m_optimize_iters=iters;}
    /**
     * setup the cut back factor
     * @param factor the cut back factor value
     */
    void setCutbackFactor(const double &factor){m_data.m_cutfactor=factor;}
    /**
     * setup the growth factor
     * @param factor the cut back factor value
     */
    void setGrowthFactor(const double &factor){m_data.m_growthfactor=factor;}
    /**
     * setup the adaptive status
     * @param flag true to enable the adaptive time stepping
     */
    void setAdaptiveFlag(const bool &flag){m_data.m_isadaptive=flag;}

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
    inline double getDt()const{return m_data.m_dt;}
    /**
     * get the inital delta t
     */
    inline double getDt0()const{return m_data.m_dt0;}
    /**
     * get max delta t
     */
    inline double getMaxDt()const{return m_data.m_dtmax;}
    /**
     * get min delta t
     */
    inline double getMinDt()const{return m_data.m_dtmin;}
    /**
     * get final time
     */
    inline double getFinalTime()const{return m_data.m_finaltime;}
    /**
     * get the cut back factor
     */
    inline double getCutbackFactor()const{return m_data.m_cutfactor;}
    /**
     * get the growth factor
     */
    inline double getGrowthFactor()const{return m_data.m_growthfactor;}
    /**
     * get the optimal iterations
     */
    inline int getOptimizeIters()const{return m_data.m_optimize_iters;}
    /**
     * get the adaptive status
     */
    inline bool isAdaptive()const{return m_data.m_isadaptive;}
    /**
     * get the time integration method
     */
    inline TimeSteppingType getTimeSteppingType()const{return m_data.m_stepping_type;}

    /**
     * solve the transient equation
     * @param t_mesh the mesh class
     * @param t_dofhandler the dof class
     * @param t_fe the fe class
     * @param t_elmtsystem the element system class
     * @param t_matesystem the material system class
     * @param t_fesystem the fe system class
     * @param t_bcsystem the boundary condition system
     * @param t_icsystem the initial condition system
     * @param t_solutionsystem the solution system class
     * @param t_equationsystem the equation system class
     * @param t_projection the projection system
     * @param t_fectrlinfo the fe control info
     * @param t_nlsolver the nonlinear solver
     * @param t_output  the output system
     * @param t_postprocess the postprocess system
     */
    bool solve(Mesh &t_mesh,DofHandler &t_dofhandler,FE &t_fe,
               ElmtSystem &t_elmtsystem,MateSystem &t_matesystem,
               FESystem &t_fesystem,
               BCSystem &t_bcsystem,
               ICSystem &t_icsystem,
               SolutionSystem &t_solutionsystem,
               EquationSystem &t_equationsystem,
               ProjectionSystem &t_projection,
               FEControlInfo &t_fectrlinfo,
               NonlinearSolver &t_nlsolver,
               OutputSystem &t_output,
               Postprocessor &t_postprocess);



private:
    TimeSteppingData m_data;/**< the time stepping data */

};