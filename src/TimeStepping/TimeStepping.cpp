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

#include "TimeStepping/TimeStepping.h"

TimeStepping::TimeStepping(){
    m_Data.m_TotalStep=1;/**< the total time step */
    m_Data.m_Dt0=1.0e-6;/**< the initial delta t */
    m_Data.m_DtMax=1.0e-3;/**< the maximum time increment */
    m_Data.m_DtMin=1.0e-13;/**< the minimum time increment */
    m_Data.m_FinalTime=1.0e-2;/**< the final simulation time */
    m_Data.m_CutFactor=0.85;/**< the cut back factor time time adaptive */
    m_Data.m_GrowthFactor=1.1;/**< the growth factor for time adaptive */
    m_Data.m_IsAdaptive=false;/**< boolean flag for adaptive */
    m_Data.m_OptimizeIters=4;/**< optimize nonlinear iterations for time adaptive */

    m_Data.m_SteppingType=TimeSteppingType::BACKWARDEULER;/**< the time stepping type */
}

void TimeStepping::applyDefaultSettings(){
    m_Data.m_TotalStep=1;/**< the total time step */
    m_Data.m_Dt0=1.0e-6;/**< the initial delta t */
    m_Data.m_DtMax=1.0e-3;/**< the maximum time increment */
    m_Data.m_DtMin=1.0e-13;/**< the minimum time increment */
    m_Data.m_FinalTime=1.0e-2;/**< the final simulation time */
    m_Data.m_CutFactor=0.85;/**< the cut back factor time time adaptive */
    m_Data.m_GrowthFactor=1.1;/**< the growth factor for time adaptive */
    m_Data.m_IsAdaptive=false;/**< boolean flag for adaptive */
    m_Data.m_OptimizeIters=4;/**< optimize nonlinear iterations for time adaptive */

    m_Data.m_SteppingType=TimeSteppingType::BACKWARDEULER;/**< the time stepping type */
}

void TimeStepping::printInfo()const{
    string str;
    char buff[69];
    
    MessagePrinter::printNormalTxt("Time stepping system information summary");

    snprintf(buff,69,"  final t=%13.5e, optimal iters=%2d",
                     getFinalTime(),getOptimizeIters());
    str=buff;
    MessagePrinter::printNormalTxt(str);

    snprintf(buff,69,"  delta t=%13.5e, max dt=%13.5e, min dt=%13.5e",
                     getDt0(),getMaxDt(),getMinDt());
    str=buff;
    MessagePrinter::printNormalTxt(str);

    snprintf(buff,69,"  growth factor=%13.5e, cut back factor=%13.5e",
                     getGrowthFactor(),getCutbackFactor());
    str=buff;
    MessagePrinter::printNormalTxt(str);
    
    if(isAdaptive()){
        MessagePrinter::printNormalTxt("  adaptive = true");
    }
    else{
        MessagePrinter::printNormalTxt("  adaptive = false");
    }

    if(getTimeSteppingType()==TimeSteppingType::BACKWARDEULER){
        MessagePrinter::printNormalTxt("  stepping method = backward euler(BE)");
    }
    else if(getTimeSteppingType()==TimeSteppingType::CRANCKNICOLSON){
        MessagePrinter::printNormalTxt("  stepping method = cranck nicolson (CN)");
    }
    else if(getTimeSteppingType()==TimeSteppingType::STATIC){
        MessagePrinter::printNormalTxt("  stepping method = static");
    }
    else if(getTimeSteppingType()==TimeSteppingType::BDF2){
        MessagePrinter::printNormalTxt("  stepping method = BDF2");
    }

    MessagePrinter::printStars();
}