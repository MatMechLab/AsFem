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
//+++ Date   : 2022.06.27
//+++ Purpose: read the boundary condition block from json
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/InputSystem.h"

bool InputSystem::readBCsBlock(nlohmann::json &t_Json,const FECell &t_FECell,const DofHandler &t_DofHandler,BCSystem &t_BCSystem){
    // now the json should already read 'bcs'
    BCBlock bcBlock;
    int nblocks=0;

    bool HasType=false;
    bool HasValue=false;
    bool HasParams=false;
    bool HasBoundary=false;
    bool HasDof=false;

    for(auto it=t_Json.begin();it!=t_Json.end();it++){
        bcBlock.reset();
        nblocks+=1;

        HasType=false;
        HasValue=false;
        HasParams=false;
        HasBoundary=false;
        HasDof=false;

        bcBlock.m_BCBlockIndex=nblocks;
        bcBlock.m_BCBlockName=it.key();
        nlohmann::json bcjson=t_Json.at(bcBlock.m_BCBlockName);//it.value();// now ejson gose into 'bc-x'
        
        if(bcjson.contains("type")){
            if(!bcjson.at("type").is_string()){
                MessagePrinter::printErrorTxt("invalid type name in bc sub block-"+to_string(nblocks)+", please check your input file");
                MessagePrinter::exitAsFem();
            }
            bcBlock.m_BCTypeName=bcjson.at("type");
            bcBlock.m_BCType=BCType::NULLBC;
            if(bcBlock.m_BCTypeName=="dirichlet"){
                bcBlock.m_BCType=BCType::DIRICHLETBC;
            }
            else if(bcBlock.m_BCTypeName=="rotateddirichlet"){
                bcBlock.m_BCType=BCType::ROTATEDDIRICHLETBC;
            }
            else if(bcBlock.m_BCTypeName=="nodaldirichlet"){
                bcBlock.m_BCType=BCType::NODALDIRICHLETBC;
            }
            else if(bcBlock.m_BCTypeName=="user1dirichlet"){
                bcBlock.m_BCType=BCType::USER1DIRICHLETBC;
            }
            else if(bcBlock.m_BCTypeName=="user2dirichlet"){
                bcBlock.m_BCType=BCType::USER2DIRICHLETBC;
            }
            else if(bcBlock.m_BCTypeName=="user3dirichlet"){
                bcBlock.m_BCType=BCType::USER3DIRICHLETBC;
            }
            else if(bcBlock.m_BCTypeName=="user4dirichlet"){
                bcBlock.m_BCType=BCType::USER4DIRICHLETBC;
            }
            else if(bcBlock.m_BCTypeName=="user5dirichlet"){
                bcBlock.m_BCType=BCType::USER5DIRICHLETBC;
            }
            else if(bcBlock.m_BCTypeName=="poisson2dbenchmark"){
                bcBlock.m_BCType=BCType::POISSON2DBENCHMARKBC;
            }
            else if(bcBlock.m_BCTypeName=="neumann"){
                bcBlock.m_BCType=BCType::NEUMANNBC;
            }
            else if(bcBlock.m_BCTypeName=="nodalneumann"){
                bcBlock.m_BCType=BCType::NODALNEUMANNBC;
            }
            else if(bcBlock.m_BCTypeName=="pressure"){
                bcBlock.m_BCType=BCType::PRESSUREBC;
            }
            else if(bcBlock.m_BCTypeName=="traction"){
                bcBlock.m_BCType=BCType::TRACTIONBC;
            }
            // for user-defined-bc
            else if(bcBlock.m_BCTypeName=="user1" && bcBlock.m_BCTypeName!="user1dirichlet"){
                bcBlock.m_BCType=BCType::USER1BC;
            }
            else if(bcBlock.m_BCTypeName=="user2" && bcBlock.m_BCTypeName!="user2dirichlet"){
                bcBlock.m_BCType=BCType::USER2BC;
            }
            else if(bcBlock.m_BCTypeName=="user3" && bcBlock.m_BCTypeName!="user3dirichlet"){
                bcBlock.m_BCType=BCType::USER3BC;
            }
            else if(bcBlock.m_BCTypeName=="user4" && bcBlock.m_BCTypeName!="user4dirichlet"){
                bcBlock.m_BCType=BCType::USER4BC;
            }
            else if(bcBlock.m_BCTypeName=="user5" && bcBlock.m_BCTypeName!="user5dirichlet"){
                bcBlock.m_BCType=BCType::USER5BC;
            }
            else{
                MessagePrinter::printErrorTxt("type="+bcBlock.m_BCTypeName
                    +" is invalid in '"+bcBlock.m_BCBlockName+"' of your 'bcs' block, please check your input file");
                MessagePrinter::exitAsFem();
            }
            HasType=true;
        }
        else{
            MessagePrinter::printErrorTxt("type= can\'t be found in '"+bcBlock.m_BCBlockName+"' in your 'bcs' block,"
                                          " please check your input file");
            MessagePrinter::exitAsFem();
            HasType=false;
        }// end-of-type-reading

        if(bcjson.contains("dofs")){
            if(bcjson.at("dofs").size()<1){
                MessagePrinter::printErrorTxt("invalid dofs name or empty dofs name in '"
                                              +bcBlock.m_BCBlockName+"' of your 'bcs' subblock, please check your input file");
                MessagePrinter::exitAsFem();
            }
            else{
                string dofname;
                for(int i=0;i<static_cast<int>(bcjson.at("dofs").size());i++){
                    if(!bcjson.at("dofs").at(i).is_string()){
                        MessagePrinter::printErrorTxt("dof-"+to_string(i+1)+"'s name is invalid in 'dofs' of '"+bcBlock.m_BCBlockName+"' in your 'bcs' subblock, please check your input file");
                        MessagePrinter::exitAsFem();
                    }
                    dofname=bcjson.at("dofs").at(i);
                    if(dofname.size()<1){
                        MessagePrinter::printErrorTxt("dof-"+to_string(i+1)+"'s name is invalid or empty in 'dofs' of '"+bcBlock.m_BCBlockName+"' in your 'bcs' subblock, please check your input file");
                        MessagePrinter::exitAsFem();
                    }
                    if(!t_DofHandler.isValidDofName(dofname)){
                        MessagePrinter::printErrorTxt("dof-"+to_string(i+1)+"'s name is invalid in 'dofs' of '"+bcBlock.m_BCBlockName+"' in your 'bcs' subblock,"
                                                      " it must be one of the names in your 'dofs' block, please check your input file");
                        MessagePrinter::exitAsFem();
                    }
                    bcBlock.m_DofsName.push_back(dofname);
                    bcBlock.m_DofIDs.push_back(t_DofHandler.getDofIDViaName(dofname));
                }
                HasDof=true;
            }
        }
        else{
            HasDof=false;
            MessagePrinter::printErrorTxt("can\'t find dofs name in '"+bcBlock.m_BCBlockName+"' of your 'bcs' subblock, please check your input file");
            MessagePrinter::exitAsFem();
        }// end-of-dofs-reading

        if(bcjson.contains("bcvalue")){
            double bcvalue=0.0;
            if(bcjson.at("bcvalue").is_number_integer()){
                bcvalue=bcjson.at("bcvalue");
                bcBlock.m_BCValue=1.0*bcvalue;
                HasValue=true;
            }
            else if(bcjson.at("bcvalue").is_number_float()){
                bcvalue=bcjson.at("bcvalue");
                bcBlock.m_BCValue=1.0*bcvalue;
                HasValue=true;
            }
            else if(bcjson.at("bcvalue").is_number_unsigned()){
                bcvalue=bcjson.at("bcvalue");
                bcBlock.m_BCValue=1.0*bcvalue;
                HasValue=true;
            }
            else if(bcjson.at("bcvalue").is_string()){
                string str;
                vector<double> numbers;
                str=bcjson.at("bcvalue");
                if(str.find("t")!=string::npos||
                   str.find("*t")!=string::npos||
                   str.find("t*")!=string::npos){
                    if(StringUtils::isValidTimeDependentExpression(str)){
                        numbers=StringUtils::splitStrNum(str);
                        if(numbers.size()<1){
                            MessagePrinter::printErrorTxt("can\'t find numbers in 'bcvalue' of '"+bcBlock.m_BCBlockName
                                                         +"' in your 'bcs' subblock, please check your input file");
                            MessagePrinter::exitAsFem();
                        }
                        bcBlock.m_BCValue=numbers[0];
                        bcBlock.m_IsTimeDependent=true;
                        HasValue=true;
                    }
                    else{
                        MessagePrinter::printErrorTxt("you are trying to define an invalid time-dependent expression in '"+bcBlock.m_BCBlockName+
                                                      "' in your 'bcs' subblock, please check your input file");
                        MessagePrinter::exitAsFem();
                    }
                }
                else{
                    numbers=StringUtils::splitStrNum(str);
                    if(numbers.size()<1){
                        MessagePrinter::printErrorTxt("can\'t find numbers in 'bcvalue' of '"+bcBlock.m_BCBlockName
                                                         +"' in your 'bcs' subblock, please check your input file");
                        MessagePrinter::exitAsFem();
                    }
                    bcBlock.m_BCValue=numbers[0];
                    HasValue=true;
                }
            }
            else{
                HasValue=false;
                MessagePrinter::printErrorTxt("unsupported expression in 'bcvalue' of '"
                    +bcBlock.m_BCBlockName+"' in your 'bcs' subblock, "
                    "it should be either a float number or an expression('1.0' or '1.0*t'), "
                    "please check your input file");
                MessagePrinter::exitAsFem();
            }
        }//end-of-bcvalue-reading

        if(bcjson.contains("side")){
            if(bcjson.at("side").size()<1){
                MessagePrinter::printErrorTxt("invalid side name in '"+bcBlock.m_BCBlockName+"' of your 'bcs' subblock, please check your input file");
                MessagePrinter::exitAsFem();
            }
            string sidename;
            for(int i=0;i<static_cast<int>(bcjson.at("side").size());i++){
                if(!bcjson.at("side").at(i).is_string()){
                    MessagePrinter::printErrorTxt("side name is invalid in 'sides' of '"+bcBlock.m_BCBlockName+"' in your 'bcs' subblock, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                sidename=bcjson.at("side").at(i);
                if(sidename.size()<1){
                    MessagePrinter::printErrorTxt("side name is invalid or empty in 'side' of '"+bcBlock.m_BCBlockName+"' in your 'bcs' subblock, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                if(!t_FECell.isFECellBCElmtPhyNameValid(sidename)){
                   MessagePrinter::printErrorTxt("side name is invalid in 'side' of '"+bcBlock.m_BCBlockName+"' in your 'bcs' subblock, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                bcBlock.m_BoundaryNameList.push_back(sidename);
            }
            HasBoundary=true;
        }
        else{
            HasBoundary=false;
            MessagePrinter::printErrorTxt("can\'t find 'side' in '"+bcBlock.m_BCBlockName+"' of your 'bcs' subblock, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        } //end-of-'side'-reading

        if(!bcjson.contains("parameters")){
            bcBlock.m_JsonParams.clear();
        }
        else{
            HasParams=true;
            bcBlock.m_JsonParams=bcjson.at("parameters");
        }//end-of-parameters-reading

        if(HasParams){}

        if(!HasType || !HasValue || !HasBoundary || !HasDof){
            MessagePrinter::printErrorTxt("information of your 'bcs' subblock('"+bcBlock.m_BCBlockName
                                          +"') is not complete,please check your input file");
            MessagePrinter::exitAsFem();
        }
        else{
            t_BCSystem.addBCBlock2List(bcBlock);
        }
    }

    return true;
}