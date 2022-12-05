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
    if(fectrlinfo.m_timesteppingtype==TimeSteppingType::STATIC){
        fectrlinfo.ctan[0]=1.0;
        fectrlinfo.ctan[1]=0.0;
        fectrlinfo.ctan[2]=0.0;
        soln.m_u_temp=U;
        soln.m_v.setToZero();
        soln.m_a.setToZero();
    }
    else if(fectrlinfo.m_timesteppingtype==TimeSteppingType::BACKWARDEULER){
        fectrlinfo.ctan[0]=1.0;
        fectrlinfo.ctan[1]=1.0/fectrlinfo.dt;
        fectrlinfo.ctan[2]=0.0;

        soln.m_u_temp=U;
        soln.m_v=soln.m_u_temp;
        soln.m_v-=soln.m_u_old;
        soln.m_v*=1.0/fectrlinfo.dt;

        soln.m_a.setToZero();
    }
    else if(fectrlinfo.m_timesteppingtype==TimeSteppingType::CRANCKNICOLSON){
        fectrlinfo.ctan[0]=0.5;
        fectrlinfo.ctan[1]=1.0/fectrlinfo.dt;
        fectrlinfo.ctan[2]=0.0;

        soln.m_u_temp=U;
        soln.m_u_temp+=soln.m_u_old;
        soln.m_u_temp*=0.5;

        soln.m_v=soln.m_u_temp;
        soln.m_v-=soln.m_u_old;
        soln.m_v*=1.0/fectrlinfo.dt;

        soln.m_a.setToZero();
    }
    else if(fectrlinfo.m_timesteppingtype==TimeSteppingType::BDF2){
        // the expression for BDF2 is:
        // [Un+2-(4/3)Un+1+(1/3)Un]/dt=(2/3)f(Un+2)

        if(fectrlinfo.CurrentStep<2){
            // for the first step, it is backward euler
            fectrlinfo.ctan[0]=1.0;
            fectrlinfo.ctan[1]=1.0/fectrlinfo.dt;
            fectrlinfo.ctan[2]=0.0;
            // calculate the current velocity
            soln.m_v =U;
            soln.m_v-=soln.m_u_old;
            soln.m_v*=fectrlinfo.ctan[1];
        }
        else{
            fectrlinfo.ctan[0]=2.0/3.0;
            fectrlinfo.ctan[1]=1.0/fectrlinfo.dt;
            fectrlinfo.ctan[2]=0.0;
            // calculate the current velocity
            soln.m_v=U;
            soln.m_v+=soln.m_u_old*(-4.0/3.0);
            soln.m_v+=soln.m_u_older*(1.0/3.0);
            soln.m_v*=fectrlinfo.ctan[1];
        }
        soln.m_a.setToZero();
    }
    else{
        MessagePrinter::printErrorTxt("Unsupported time stepping method for time derivates (TimeSteppingTool)");
        MessagePrinter::exitAsFem();
    }
}