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
//+++ Date   : 2022.05.09
//+++ Purpose: read the dofs block from json
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/InputSystem.h"

bool InputSystem::readDofsBlock(nlohmann::json &t_json,DofHandler &t_dofhandler){
    // now the json already read "dofs"
    t_dofhandler.init();
    if(t_json.contains("names")){
        if(t_json.at("names").size()<1){
            MessagePrinter::printErrorTxt("can\'t find dof names or empty dof names in 'names' of your 'dofs' block, please check your input file");
            MessagePrinter::exitAsFem();
        }
        string dofname;
        for(int i=0;i<static_cast<int>(t_json.at("names").size());i++){
            if(!t_json.at("names").at(i).is_string()){
                MessagePrinter::printErrorTxt(to_string(i+1)+"-th dofs name is invalid in 'names' of your 'dofs' block, please check your input file");
                MessagePrinter::exitAsFem();
            }
            dofname=t_json.at("names").at(i);
            if(dofname.size()<1){
                MessagePrinter::printErrorTxt(to_string(i+1)+"-th dofs name is invalid or empty in 'names' of your 'dofs' block, please check your input file");
                MessagePrinter::exitAsFem();
            }
            t_dofhandler.addDofName2List(dofname);
        }
    }
    else{
        MessagePrinter::printErrorTxt("can\'t find 'names' in your 'dofs' block, please check your input file");
        MessagePrinter::exitAsFem();
    }
    return true;
}