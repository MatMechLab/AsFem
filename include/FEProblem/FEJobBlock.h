//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
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

#include <iostream>
#include <string>

#include "Utils/MessagePrinter.h"
#include "FEProblem/FEJobType.h"

using namespace std;


class FEJobBlock{
public:
    FEJobType _jobType=FEJobType::STATIC;
    string   _jobTypeName="static";
    bool _IsDebug=true,_IsDepDebug=false;


    void Init(){
        _jobType=FEJobType::STATIC;
        _jobTypeName="static";
        _IsDebug=true;
        _IsDepDebug=false;
    }

    void PrintJobInfo(){
        MessagePrinter::PrintNormalTxt("Job information summary:");
        MessagePrinter::PrintNormalTxt("  job type="+_jobTypeName);
        if(_IsDebug){
            if(_IsDepDebug){
                MessagePrinter::PrintNormalTxt("  debug dep print is enabled");
            }
            else{
                MessagePrinter::PrintNormalTxt("  debug print is enabled");
            }
        }
        MessagePrinter::PrintDashLine();
    }
};