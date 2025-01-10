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
//+++ Date   : 2022.08.25
//+++ Purpose: Implement aux functions for time stepping
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "SolutionSystem/SolutionSystem.h"
#include "TimeStepping/TimeSteppingType.h"
#include "FEProblem/FEControlInfo.h"

/**
 * compute the different orders of time derivatives according to choosen time stepping method
 * @param fectrlinfo fe control information
 * @param U the intermediate PETSc Vec for current solution during iteration
 * @param soln the solution system
 */
void computeTimeDerivatives(FEControlInfo &fectrlinfo,
                            const Vec &U,
                            SolutionSystem &soln){
    if(fectrlinfo.m_TimesteppingType==TimeSteppingType::STATIC){
        fectrlinfo.Ctan[0]=1.0;
        fectrlinfo.Ctan[1]=0.0;
        fectrlinfo.Ctan[2]=0.0;
        soln.m_Utemp=U;
        soln.m_V.setToZero();
        soln.m_A.setToZero();
    }
    else if(fectrlinfo.m_TimesteppingType==TimeSteppingType::BACKWARDEULER){
        fectrlinfo.Ctan[0]=1.0;
        fectrlinfo.Ctan[1]=1.0/fectrlinfo.Dt;
        fectrlinfo.Ctan[2]=0.0;

        soln.m_Utemp=U;
        soln.m_V=soln.m_Utemp;
        soln.m_V-=soln.m_Uold;
        soln.m_V*=1.0/fectrlinfo.Dt;

        soln.m_A.setToZero();
    }
    else if(fectrlinfo.m_TimesteppingType==TimeSteppingType::CRANCKNICOLSON){
        fectrlinfo.Ctan[0]=0.5;
        fectrlinfo.Ctan[1]=1.0/fectrlinfo.Dt;
        fectrlinfo.Ctan[2]=0.0;

        soln.m_Utemp=U;
        soln.m_Utemp+=soln.m_Uold;
        soln.m_Utemp*=0.5;

        soln.m_V=soln.m_Utemp;
        soln.m_V-=soln.m_Uold;
        soln.m_V*=1.0/fectrlinfo.Dt;

        soln.m_A.setToZero();
    }
    else if(fectrlinfo.m_TimesteppingType==TimeSteppingType::BDF2){
        // the expression for BDF2 is:
        // [Un+2-(4/3)Un+1+(1/3)Un]/dt=(2/3)f(Un+2)

        if(fectrlinfo.CurrentStep<2){
            // for the first step, it is backward euler
            fectrlinfo.Ctan[0]=1.0;
            fectrlinfo.Ctan[1]=1.0/fectrlinfo.Dt;
            fectrlinfo.Ctan[2]=0.0;
            // calculate the current velocity
            soln.m_V =U;
            soln.m_V-=soln.m_Uold;
            soln.m_V*=fectrlinfo.Ctan[1];
        }
        else{
            fectrlinfo.Ctan[0]=2.0/3.0;
            fectrlinfo.Ctan[1]=1.0/fectrlinfo.Dt;
            fectrlinfo.Ctan[2]=0.0;
            // calculate the current velocity
            soln.m_V=U;
            soln.m_V+=soln.m_Uold*(-4.0/3.0);
            soln.m_V+=soln.m_Uolder*(1.0/3.0);
            soln.m_V*=fectrlinfo.Ctan[1];
        }
        soln.m_A.setToZero();
    }
    else{
        MessagePrinter::printErrorTxt("Unsupported time stepping method for time derivates (TimeSteppingTool)");
        MessagePrinter::exitAsFem();
    }
}