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
//+++ Purpose: read the nonlinear solver block 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/InputSystem.h"

bool InputSystem::readOutputBlock(nlohmann::json &t_json,OutputSystem &t_output){
    // the json already contains "output" block
    bool HasType;
    if(t_json.contains("type")){
        if(!t_json.at("type").is_string()){
            MessagePrinter::printErrorTxt("the type name of your output block is not a valid string");
            return false;
        }
        string filetypename=t_json.at("type");

        if(filetypename=="default"){
            t_output.setFileFormat(ResultFileFormat::VTU);
        }
        else if(filetypename=="vtu"){
            t_output.setFileFormat(ResultFileFormat::VTU);
        }
        else if(filetypename=="vtk"){
            t_output.setFileFormat(ResultFileFormat::VTK);
        }
        else if(filetypename=="csv"){
            t_output.setFileFormat(ResultFileFormat::CSV);
        }
        else{
            MessagePrinter::printErrorTxt("Unsupported 'type' (file format) in the output block, please check your input file");
            return false;
        }
        HasType=true;
    }
    else{
        MessagePrinter::printWarningTxt("no 'type' (file format) found in the output block, the default one (VTU) will be used");
        t_output.setFileFormat(ResultFileFormat::VTU);
        HasType=true;
    }

    if(t_json.contains("interval")){
        if(!t_json.at("interval").is_number_integer()){
            MessagePrinter::printErrorTxt("the interval value is not valid, please check your output block");
            return false;
        }
        int step=t_json.at("interval");
        if(step<1){
            MessagePrinter::printErrorTxt("interval="+to_string(step)+" is invalid, please check your output block");
            return false;
        }
        t_output.setIntervalNum(step);
    }
    else{
        MessagePrinter::printWarningTxt("no 'interval' found in the output block, the default value (1) will be used");
        t_output.setIntervalNum(1);
    }

    return HasType;
}