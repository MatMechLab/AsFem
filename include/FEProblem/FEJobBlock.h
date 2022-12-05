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
    FEJobType m_jobtype=FEJobType::STATIC;/**< for the job type, i.e., static, transient. */
    string   m_jobtypename="static";/**< the job type name */
    bool m_isdebug=true;/**< message print level */
    bool m_isdepdebug=false;/**< for the dep message print */

    /**
     * init the job block
     */
    void init(){
        m_jobtype=FEJobType::STATIC;
        m_jobtypename="static";
        m_isdebug=true;
        m_isdepdebug=false;
    }
    /**
     * print out the job block information
     */
    void printJobInfo(){
        MessagePrinter::printNormalTxt("Job information summary:");
        MessagePrinter::printNormalTxt("  job type="+m_jobtypename);
        if(m_isdebug){
            if(m_isdepdebug){
                MessagePrinter::printNormalTxt("  dep message print is enabled");
            }
            else{
                MessagePrinter::printNormalTxt("  message print is enabled");
            }
        }
        MessagePrinter::printStars();
    }
};