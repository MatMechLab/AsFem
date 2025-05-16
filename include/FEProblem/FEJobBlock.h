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
//+++ Purpose: Define the job block for FEM analysis
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "Utils/MessagePrinter.h"
#include "FEProblem/FEJobType.h"

/**
 * This class defines the basic info for a job block
 */
class FEJobBlock{
public:
    FEJobType m_JobType=FEJobType::STATIC;/**< for the job type, i.e., static, transient. */
    string   m_JobTypeName="static";/**< the job type name */
    bool m_IsDebug=true;/**< message print level */
    bool m_IsDepDebug=false;/**< for the dep message print */

    /**
     * init the job block
     */
    void init(){
        m_JobType=FEJobType::STATIC;
        m_JobTypeName="static";
        m_IsDebug=true;
        m_IsDepDebug=false;
    }
    /**
     * print out the job block information
     */
    void printJobInfo(){
        MessagePrinter::printNormalTxt("Job information summary:");
        MessagePrinter::printNormalTxt("  job type="+m_JobTypeName);
        if(m_IsDebug){
            if(m_IsDepDebug){
                MessagePrinter::printNormalTxt("  dep message print is enabled");
            }
            else{
                MessagePrinter::printNormalTxt("  message print is enabled");
            }
        }
        MessagePrinter::printStars();
    }
};