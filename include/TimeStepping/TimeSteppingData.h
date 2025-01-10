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
//+++ Date   : 2022.08.26
//+++ Purpose: Defines the data for time stepping
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

/**
 * time stepping data
 */
struct TimeSteppingData{
    int m_TotalStep;/**< the total time step */
    double m_Dt0=1.0e-6;/**< the initial delta t */
    double m_Dt=1.0e-6;/**< the current delta t */
    double m_DtMax=1.0e1;/**< the maximum time increment */
    double m_DtMin=1.0e-13;/**< the minimum time increment */
    double m_FinalTime;/**< the final simulation time */
    double m_CutFactor;/**< the cut back factor time time adaptive */
    double m_GrowthFactor;/**< the growth factor for time adaptive */
    bool m_IsAdaptive;/**< boolean flag for adaptive */
    int m_OptimizeIters=3;/**< optimize nonlinear iterations for time adaptive */

    TimeSteppingType m_SteppingType;/**< the time stepping type */
};
