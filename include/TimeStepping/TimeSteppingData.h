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
//+++ Date   : 2022.08.26
//+++ Purpose: Defines the data for time stepping
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

/**
 * time stepping data
 */
struct TimeSteppingData{
    int m_totalstep;/**< the total time step */
    double m_dt0=1.0e-6;/**< the initial delta t */
    double m_dt=1.0e-6;/**< the current delta t */
    double m_dtmax=1.0e1;/**< the maximum time increment */
    double m_dtmin=1.0e-13;/**< the minimum time increment */
    double m_finaltime;/**< the final simulation time */
    double m_cutfactor;/**< the cut back factor time time adaptive */
    double m_growthfactor;/**< the growth factor for time adaptive */
    bool m_isadaptive;/**< boolean flag for adaptive */
    int m_optimize_iters=3;/**< optimize nonlinear iterations for time adaptive */

    TimeSteppingType m_stepping_type;/**< the time stepping type */
};
