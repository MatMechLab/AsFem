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
//+++ Date   : 2022.06.27
//+++ Purpose: read the boundary condition block from json
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/InputSystem.h"

bool InputSystem::readBCsBlock(nlohmann::json &t_json,const Mesh &t_mesh,const DofHandler &t_dofhandler,BCSystem &t_bcsystem){
    // now the json should already read 'bcs'
    BCBlock bcBlock;
    int nblocks=0;

    bool HasType=false;
    bool HasValue=false;
    bool HasParams=false;
    bool HasBoundary=false;
    bool HasDof=false;

    for(auto it=t_json.begin();it!=t_json.end();it++){
        bcBlock.reset();
        nblocks+=1;

        HasType=false;
        HasValue=false;
        HasParams=false;
        HasBoundary=false;
        HasDof=false;

        bcBlock.m_bcBlockIndex=nblocks;
        bcBlock.m_bcBlockName=it.key();
        nlohmann::json bcjson=t_json.at(bcBlock.m_bcBlockName);//it.value();// now ejson gose into 'bc-x'
        
        if(bcjson.contains("type")){
            if(!bcjson.at("type").is_string()){
                MessagePrinter::printErrorTxt("invalid type name in bc sub block-"+to_string(nblocks)+", please check your input file");
                MessagePrinter::exitAsFem();
            }
            bcBlock.m_bcTypeName=bcjson.at("type");
            bcBlock.m_bcType=BCType::NULLBC;
            if(bcBlock.m_bcTypeName=="dirichlet"){
                bcBlock.m_bcType=BCType::DIRICHLETBC;
            }
            else if(bcBlock.m_bcTypeName=="rotateddirichlet"){
                bcBlock.m_bcType=BCType::ROTATEDDIRICHLETBC;
            }
            else if(bcBlock.m_bcTypeName=="nodaldirichlet"){
                bcBlock.m_bcType=BCType::NODALDIRICHLETBC;
            }
            else if(bcBlock.m_bcTypeName=="user1dirichlet"){
                bcBlock.m_bcType=BCType::USER1DIRICHLETBC;
            }
            else if(bcBlock.m_bcTypeName=="user2dirichlet"){
                bcBlock.m_bcType=BCType::USER2DIRICHLETBC;
            }
            else if(bcBlock.m_bcTypeName=="user3dirichlet"){
                bcBlock.m_bcType=BCType::USER3DIRICHLETBC;
            }
            else if(bcBlock.m_bcTypeName=="user4dirichlet"){
                bcBlock.m_bcType=BCType::USER4DIRICHLETBC;
            }
            else if(bcBlock.m_bcTypeName=="user5dirichlet"){
                bcBlock.m_bcType=BCType::USER5DIRICHLETBC;
            }
            else if(bcBlock.m_bcTypeName=="poisson2dbenchmark"){
                bcBlock.m_bcType=BCType::POISSON2DBENCHMARKBC;
            }
            else if(bcBlock.m_bcTypeName=="neumann"){
                bcBlock.m_bcType=BCType::NEUMANNBC;
            }
            else if(bcBlock.m_bcTypeName=="nodalneumann"){
                bcBlock.m_bcType=BCType::NODALNEUMANNBC;
            }
            else if(bcBlock.m_bcTypeName=="pressure"){
                bcBlock.m_bcType=BCType::PRESSUREBC;
            }
            else if(bcBlock.m_bcTypeName=="traction"){
                bcBlock.m_bcType=BCType::TRACTIONBC;
            }
            // for user-defined-bc
            else if(bcBlock.m_bcTypeName=="user1" && bcBlock.m_bcTypeName!="user1dirichlet"){
                bcBlock.m_bcType=BCType::USER1BC;
            }
            else if(bcBlock.m_bcTypeName=="user2" && bcBlock.m_bcTypeName!="user2dirichlet"){
                bcBlock.m_bcType=BCType::USER2BC;
            }
            else if(bcBlock.m_bcTypeName=="user3" && bcBlock.m_bcTypeName!="user3dirichlet"){
                bcBlock.m_bcType=BCType::USER3BC;
            }
            else if(bcBlock.m_bcTypeName=="user4" && bcBlock.m_bcTypeName!="user4dirichlet"){
                bcBlock.m_bcType=BCType::USER4BC;
            }
            else if(bcBlock.m_bcTypeName=="user5" && bcBlock.m_bcTypeName!="user5dirichlet"){
                bcBlock.m_bcType=BCType::USER5BC;
            }
            else if(bcBlock.m_bcTypeName=="user6" && bcBlock.m_bcTypeName!="user6dirichlet"){
                bcBlock.m_bcType=BCType::USER6BC;
            }
            else if(bcBlock.m_bcTypeName=="user7" && bcBlock.m_bcTypeName!="user7dirichlet"){
                bcBlock.m_bcType=BCType::USER7BC;
            }
            else if(bcBlock.m_bcTypeName=="user8" && bcBlock.m_bcTypeName!="user8dirichlet"){
                bcBlock.m_bcType=BCType::USER8BC;
            }
            else if(bcBlock.m_bcTypeName=="user9" && bcBlock.m_bcTypeName!="user9dirichlet"){
                bcBlock.m_bcType=BCType::USER9BC;
            }
            else if(bcBlock.m_bcTypeName=="user10" && bcBlock.m_bcTypeName!="user10dirichlet"){
                bcBlock.m_bcType=BCType::USER10BC;
            }
            else{
                MessagePrinter::printErrorTxt("type="+bcBlock.m_bcTypeName
                    +" is invalid in ["+bcBlock.m_bcBlockName+"] of your [bcs] block, please check your input file");
                MessagePrinter::exitAsFem();
            }
            HasType=true;
        }
        else{
            MessagePrinter::printErrorTxt("type= can\'t be found in ["+bcBlock.m_bcBlockName+"] in your [bcs] block,"
                                          " please check your input file");
            MessagePrinter::exitAsFem();
            HasType=false;
        }// end-of-type-reading

        if(bcjson.contains("dofs")){
            if(bcjson.at("dofs").size()<1){
                MessagePrinter::printErrorTxt("invalid dofs name or empty dofs name in ["
                                              +bcBlock.m_bcBlockName+"] of your [bcs] subblock, please check your input file");
                MessagePrinter::exitAsFem();
            }
            else{
                string dofname;
                for(int i=0;i<static_cast<int>(bcjson.at("dofs").size());i++){
                    if(!bcjson.at("dofs").at(i).is_string()){
                        MessagePrinter::printErrorTxt("dof-"+to_string(i+1)+"'s name is invalid in 'dofs' of ["+bcBlock.m_bcBlockName+"] in your [bcs] subblock, please check your input file");
                        MessagePrinter::exitAsFem();
                    }
                    dofname=bcjson.at("dofs").at(i);
                    if(dofname.size()<1){
                        MessagePrinter::printErrorTxt("dof-"+to_string(i+1)+"'s name is invalid or empty in 'dofs' of ["+bcBlock.m_bcBlockName+"] in your [bcs] subblock, please check your input file");
                        MessagePrinter::exitAsFem();
                    }
                    if(!t_dofhandler.isValidDofName(dofname)){
                        MessagePrinter::printErrorTxt("dof-"+to_string(i+1)+"'s name is invalid in 'dofs' of ["+bcBlock.m_bcBlockName+"] in your [bcs] subblock,"
                                                      " it must be one of the names in your 'dofs' block, please check your input file");
                        MessagePrinter::exitAsFem();
                    }
                    bcBlock.m_dofsName.push_back(dofname);
                    bcBlock.m_dofIDs.push_back(t_dofhandler.getDofIDViaName(dofname));
                }
                HasDof=true;
            }
        }
        else{
            HasDof=false;
            MessagePrinter::printErrorTxt("can\'t find dofs name in ["+bcBlock.m_bcBlockName+"] of your [bcs] subblock, please check your input file");
            MessagePrinter::exitAsFem();
        }// end-of-dofs-reading

        if(bcjson.contains("bcvalue")){
            double bcvalue=0.0;
            if(bcjson.at("bcvalue").is_number_integer()){
                bcvalue=bcjson.at("bcvalue");
                bcBlock.m_bcValue=1.0*bcvalue;
                HasValue=true;
            }
            else if(bcjson.at("bcvalue").is_number_float()){
                bcvalue=bcjson.at("bcvalue");
                bcBlock.m_bcValue=1.0*bcvalue;
                HasValue=true;
            }
            else if(bcjson.at("bcvalue").is_number_unsigned()){
                bcvalue=bcjson.at("bcvalue");
                bcBlock.m_bcValue=1.0*bcvalue;
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
                            MessagePrinter::printErrorTxt("can\'t find numbers in 'bcvalue' of ["+bcBlock.m_bcBlockName
                                                         +"] in your [bcs] subblock, please check your input file");
                            MessagePrinter::exitAsFem();
                        }
                        bcBlock.m_bcValue=numbers[0];
                        bcBlock.m_isTimeDependent=true;
                        HasValue=true;
                    }
                    else{
                        MessagePrinter::printErrorTxt("you are trying to define an invalid time-dependent expression in ["+bcBlock.m_bcBlockName+
                                                      "] in your [bcs] subblock, please check your input file");
                        MessagePrinter::exitAsFem();
                    }
                }
                else{
                    numbers=StringUtils::splitStrNum(str);
                    if(numbers.size()<1){
                        MessagePrinter::printErrorTxt("can\'t find numbers in 'bcvalue' of ["+bcBlock.m_bcBlockName
                                                         +"] in your [bcs] subblock, please check your input file");
                        MessagePrinter::exitAsFem();
                    }
                    bcBlock.m_bcValue=numbers[0];
                    HasValue=true;
                }
            }
            else{
                HasValue=false;
                MessagePrinter::printErrorTxt("unsupported expression in 'bcvalue' of ["
                    +bcBlock.m_bcBlockName+"] in your [bcs] subblock, "
                    "it should be either a float number or an expression('1.0' or '1.0*t'), "
                    "please check your input file");
                MessagePrinter::exitAsFem();
            }
        }//end-of-bcvalue-reading

        if(bcjson.contains("side")){
            if(bcjson.at("side").size()<1){
                MessagePrinter::printErrorTxt("invalid side name in ["+bcBlock.m_bcBlockName+"] of your [bcs] subblock, please check your input file");
                MessagePrinter::exitAsFem();
            }
            string sidename;
            for(int i=0;i<static_cast<int>(bcjson.at("side").size());i++){
                if(!bcjson.at("side").at(i).is_string()){
                    MessagePrinter::printErrorTxt("side name is invalid in 'sides' of ["+bcBlock.m_bcBlockName+"] in your [bcs] subblock, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                sidename=bcjson.at("side").at(i);
                if(sidename.size()<1){
                    MessagePrinter::printErrorTxt("side name is invalid or empty in 'side' of ["+bcBlock.m_bcBlockName+"] in your [bcs] subblock, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                if(!t_mesh.isBCElmtPhyNameValid(sidename)){
                   MessagePrinter::printErrorTxt("side name is invalid in 'side' of ["+bcBlock.m_bcBlockName+"] in your [bcs] subblock, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                bcBlock.m_boundaryNameList.push_back(sidename);
            }
            HasBoundary=true;
        }
        else{
            HasBoundary=false;
            MessagePrinter::printErrorTxt("can\'t find 'side' in ["+bcBlock.m_bcBlockName+"] of your [bcs] subblock, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        } //end-of-'side'-reading

        if(!bcjson.contains("parameters")){
            bcBlock.m_json_params.clear();
        }
        else{
            HasParams=true;
            bcBlock.m_json_params=bcjson.at("parameters");
        }//end-of-parameters-reading

        if(HasParams){}

        if(!HasType || !HasValue || !HasBoundary || !HasDof){
            MessagePrinter::printErrorTxt("information of your [bcs] subblock(["+bcBlock.m_bcBlockName
                                          +"]) is not complete,please check your input file");
            MessagePrinter::exitAsFem();
        }
        else{
            t_bcsystem.addBCBlock2List(bcBlock);
        }
    }

    return true;
}