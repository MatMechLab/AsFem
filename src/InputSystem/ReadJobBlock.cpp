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
//+++ Date   : 2022.08.07
//+++ Purpose: read the job block 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/InputSystem.h"

bool InputSystem::readJobBlock(nlohmann::json &t_json,FEJobBlock &t_jobblock){
    // the json already contains "job" block
    bool HasType;
    if(t_json.contains("type")){
        if(!t_json.at("type").is_string()){
            MessagePrinter::printErrorTxt("the type name of your job block is not a valid string");
            return false;
        }
        string jobtypename=t_json.at("type");
        if(jobtypename=="static"){
            t_jobblock.m_jobtypename="static";
            t_jobblock.m_jobtype=FEJobType::STATIC;
        }
        else if(jobtypename=="transient"){
            t_jobblock.m_jobtypename="transient";
            t_jobblock.m_jobtype=FEJobType::TRANSIENT;
        }
        else{
            MessagePrinter::printErrorTxt("type="+jobtypename+" is invalid in [job] block, please check your input file");
            MessagePrinter::exitAsFem();
        }
        HasType=true;
    }
    else{
        MessagePrinter::printErrorTxt("can\'t find 'type'= in your job block, please check your input file");
        MessagePrinter::exitAsFem();
    }

    if(t_json.contains("print")){
        if(!t_json.at("print").is_string()){
            MessagePrinter::printErrorTxt("the 'print' option in your job block is not a valid string");
            return false;
        }
        string print=t_json.at("print");
        if(print=="on"){
            t_jobblock.m_isdebug=true;
            t_jobblock.m_isdepdebug=false;
        }
        else if(print=="dep"){
            t_jobblock.m_isdebug=true;
            t_jobblock.m_isdepdebug=true;
        }
        else if(print=="off"){
            t_jobblock.m_isdebug=false;
            t_jobblock.m_isdepdebug=false;
        }
        else{
            MessagePrinter::printErrorTxt("unsupported option for 'print' in your job block, please check your input file" );
            MessagePrinter::exitAsFem();
        }
    }
    else{
        MessagePrinter::printWarningTxt("can\'t find 'print' in your job block, then the default linear solver (gmres) will be used");
        t_jobblock.m_isdebug=true;
        t_jobblock.m_isdepdebug=false;
    }


    return HasType;
}