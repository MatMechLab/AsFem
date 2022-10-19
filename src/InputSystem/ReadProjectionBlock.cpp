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
//+++ Date   : 2022.05.07
//+++ Purpose: read the "projection" block from the input file
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/InputSystem.h"

bool InputSystem::readProjectionBlock(nlohmann::json &t_json,const DofHandler &t_dofhandler,ProjectionSystem &t_projsystem){
    // here json already contains "projection"


    t_projsystem.setScalarMaterialNum(0);
    t_projsystem.setVectorMaterialNum(0);
    t_projsystem.setRank2TensorMaterialNum(0);
    t_projsystem.setRank4TensorMaterialNum(0);

    if(t_json.contains("type")){
        string methodname=t_json.at("type");
        if(methodname=="default"){
            t_projsystem.setProjectionType(ProjectionType::DEFAULT);
        }
        else if(methodname=="leastsquare"){
            t_projsystem.setProjectionType(ProjectionType::LEASTSQUARE);
        }
        else{
            MessagePrinter::printErrorTxt("type="+methodname+" is unsupported in your 'projection' block, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    else{
        MessagePrinter::printWarningTxt("can\'t find 'type=' option in your 'projection' block, then the default one will be used");
    }

    if(t_json.contains("scalarmate")){
        string matename;
        for(int i=0;i<static_cast<int>(t_json.at("scalarmate").size());i++){
            if(!t_json.at("scalarmate").at(i).is_string()){
                MessagePrinter::printErrorTxt(to_string(i+1)+"-th mate name is invalid in 'scalarmate' of your 'projection' block, please check your input file");
                MessagePrinter::exitAsFem();
            }
            matename=t_json.at("scalarmate").at(i);
            if(matename.size()<1){
                MessagePrinter::printErrorTxt(to_string(i+1)+"-th mate name is invalid or empty in 'scalarmate' of your 'projection' block, please check your input file");
                MessagePrinter::exitAsFem();
            }
            if(t_dofhandler.isValidDofName(matename)){
                MessagePrinter::printErrorTxt("scalar material name("+matename+") in your 'projection' block should be different with your dofs name, please check your input file");
                MessagePrinter::exitAsFem();
            }
            t_projsystem.addScalarMateName2List(matename);
        }
    }

    if(t_json.contains("vectormate")){
        string matename;
        for(int i=0;i<static_cast<int>(t_json.at("vectormate").size());i++){
            if(!t_json.at("vectormate").at(i).is_string()){
                MessagePrinter::printErrorTxt(to_string(i+1)+"-th mate name is invalid in 'vectormate' of your 'projection' block, please check your input file");
                MessagePrinter::exitAsFem();
            }
            matename=t_json.at("vectormate").at(i);
            if(matename.size()<1){
                MessagePrinter::printErrorTxt(to_string(i+1)+"-th mate name is invalid or empty in 'vectormate' of your 'projection' block, please check your input file");
                MessagePrinter::exitAsFem();
            }
            if(t_dofhandler.isValidDofName(matename)){
                MessagePrinter::printErrorTxt("vector material name("+matename+") in your 'projection' block should be different with your dofs name, please check your input file");
                MessagePrinter::exitAsFem();
            }
            t_projsystem.addVectorMateName2List(matename);
        }
    }

    if(t_json.contains("rank2mate")){
        string matename;
        for(int i=0;i<static_cast<int>(t_json.at("rank2mate").size());i++){
            if(!t_json.at("rank2mate").at(i).is_string()){
                MessagePrinter::printErrorTxt(to_string(i+1)+"-th mate name is invalid in 'rank2mate' of your 'projection' block, please check your input file");
                MessagePrinter::exitAsFem();
            }
            matename=t_json.at("rank2mate").at(i);
            if(matename.size()<1){
                MessagePrinter::printErrorTxt(to_string(i+1)+"-th mate name is invalid or empty in 'rank2mate' of your 'projection' block, please check your input file");
                MessagePrinter::exitAsFem();
            }
            if(t_dofhandler.isValidDofName(matename)){
                MessagePrinter::printErrorTxt("rank-2 tensor material name("+matename+") in your 'projection' block should be different with your dofs name, please check your input file");
                MessagePrinter::exitAsFem();
            }
            t_projsystem.addRank2MateName2List(matename);
        }
    }

    if(t_json.contains("rank4mate")){
        string matename;
        for(int i=0;i<static_cast<int>(t_json.at("rank4mate").size());i++){
            if(!t_json.at("rank4mate").at(i).is_string()){
                MessagePrinter::printErrorTxt(to_string(i+1)+"-th mate name is invalid in 'rank4mate' of your 'projection' block, please check your input file");
                MessagePrinter::exitAsFem();
            }
            matename=t_json.at("rank4mate").at(i);
            if(matename.size()<1){
                MessagePrinter::printErrorTxt(to_string(i+1)+"-th mate name is invalid or empty in 'rank4mate' of your 'projection' block, please check your input file");
                MessagePrinter::exitAsFem();
            }
            if(t_dofhandler.isValidDofName(matename)){
                MessagePrinter::printErrorTxt("rank-4 tensor material name("+matename+") in your 'projection' block should be different with your dofs name, please check your input file");
                MessagePrinter::exitAsFem();
            }
            t_projsystem.addRank4MateName2List(matename);
        }
    }
    return true;
    
}
