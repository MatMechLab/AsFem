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

bool InputSystem::readICsBlock(nlohmann::json &t_json,const FECell &t_fecell,const DofHandler &t_dofhandler,ICSystem &t_icsystem){
    // now the json should already read 'bcs'
    ICBlock icBlock;
    int nblocks=0;

    bool HasType=false;
    bool HasParams=false;
    bool HasDomain=false;
    bool HasDof=false;

    for(auto it=t_json.begin();it!=t_json.end();it++){
        icBlock.reset();
        nblocks+=1;

        HasType=false;
        HasParams=false;
        HasDomain=false;
        HasDof=false;

        icBlock.m_ICBlockIndex=nblocks;
        icBlock.m_ICBlockName=it.key();
        nlohmann::json icjson=t_json.at(icBlock.m_ICBlockName);//it.value();// now ejson gose into 'ic-x'

        if(icjson.contains("type")){
            if(!icjson.at("type").is_string()){
                MessagePrinter::printErrorTxt("invalid type name in '"+icBlock.m_ICBlockName+"' in 'ics' subblock, please check your input file");
                MessagePrinter::exitAsFem();
            }
            icBlock.m_ICTypeName=icjson.at("type");
            icBlock.m_ICType=ICType::NULLIC;
            if(icBlock.m_ICTypeName=="const"){
                icBlock.m_ICType=ICType::CONSTIC;
            }
            else if(icBlock.m_ICTypeName=="random"){
                icBlock.m_ICType=ICType::RANDOMIC;
            }
            else if(icBlock.m_ICTypeName=="rectangle"){
                icBlock.m_ICType=ICType::RECTANGLEIC;
            }
            else if(icBlock.m_ICTypeName=="cubic"){
                icBlock.m_ICType=ICType::CUBICIC;
            }
            else if(icBlock.m_ICTypeName=="circle"){
                icBlock.m_ICType=ICType::CIRCLEIC;
            }
            else if(icBlock.m_ICTypeName=="smoothcircle"){
                icBlock.m_ICType=ICType::SMOOTHCIRCLEIC;
            }
            else if(icBlock.m_ICTypeName=="spherical"){
                icBlock.m_ICType=ICType::SPHERICALIC;
            }
            // for user-defined ics
            else if(icBlock.m_ICTypeName=="user1"){
                icBlock.m_ICType=ICType::USER1IC;
            }
            else if(icBlock.m_ICTypeName=="user2"){
                icBlock.m_ICType=ICType::USER2IC;
            }
            else if(icBlock.m_ICTypeName=="user3"){
                icBlock.m_ICType=ICType::USER3IC;
            }
            else if(icBlock.m_ICTypeName=="user4"){
                icBlock.m_ICType=ICType::USER4IC;
            }
            else if(icBlock.m_ICTypeName=="user5"){
                icBlock.m_ICType=ICType::USER5IC;
            }
            else if(icBlock.m_ICTypeName=="user6"){
                icBlock.m_ICType=ICType::USER6IC;
            }
            else if(icBlock.m_ICTypeName=="user7"){
                icBlock.m_ICType=ICType::USER7IC;
            }
            else if(icBlock.m_ICTypeName=="user8"){
                icBlock.m_ICType=ICType::USER8IC;
            }
            else if(icBlock.m_ICTypeName=="user9"){
                icBlock.m_ICType=ICType::USER9IC;
            }
            else if(icBlock.m_ICTypeName=="user10"){
                icBlock.m_ICType=ICType::USER10IC;
            }
            else{
                MessagePrinter::printErrorTxt("type="+icBlock.m_ICTypeName
                    +" is invalid in '"+icBlock.m_ICBlockName+"' of your 'ics' subblock, please check your input file");
                MessagePrinter::exitAsFem();
            }
            HasType=true;
        }//end-of-type-reading

        if(icjson.contains("dofs")){
            if(icjson.at("dofs").size()<1){
                MessagePrinter::printErrorTxt("invalid dofs name or empty dofs name in '"
                                              +icBlock.m_ICBlockName+"' of your 'ics' subblock, please check your input file");
                MessagePrinter::exitAsFem();
            }
            else{
                string dofname;
                for(int i=0;i<static_cast<int>(icjson.at("dofs").size());i++){
                    if(!icjson.at("dofs").at(i).is_string()){
                        MessagePrinter::printErrorTxt("dof-"+to_string(i+1)+"'s name is invalid in 'dofs' of '"+icBlock.m_ICBlockName+"' in your 'bcs' subblock, please check your input file");
                        MessagePrinter::exitAsFem();
                    }
                    dofname=icjson.at("dofs").at(i);
                    if(dofname.size()<1){
                        MessagePrinter::printErrorTxt("dof-"+to_string(i+1)+"'s name is invalid or empty in 'dofs' of '"+icBlock.m_ICBlockName+"' in your 'bcs' subblock, please check your input file");
                        MessagePrinter::exitAsFem();
                    }
                    if(!t_dofhandler.isValidDofName(dofname)){
                        MessagePrinter::printErrorTxt("dof-"+to_string(i+1)+"'s name is invalid in 'dofs' of '"+icBlock.m_ICBlockName+"' in your 'bcs' subblock,"
                                                      " it must be one of the names in your 'dofs' block, please check your input file");
                        MessagePrinter::exitAsFem();
                    }
                    icBlock.m_DofNames.push_back(dofname);
                    icBlock.m_DofIDs.push_back(t_dofhandler.getDofIDViaName(dofname));
                }
                HasDof=true;
            }
        }
        else{
            HasDof=false;
            MessagePrinter::printErrorTxt("can\'t find dofs name in '"+icBlock.m_ICBlockName+"' of your 'ics' subblock, please check your input file");
            MessagePrinter::exitAsFem();
        }// end-of-dofs-reading

        if(icjson.contains("domain")){
            if(icjson.at("domain").size()<1){
                MessagePrinter::printErrorTxt("invalid or empty domain name in '"+icBlock.m_ICBlockName+"' of your 'ics' subblock, please check your input file");
                MessagePrinter::exitAsFem();
            }
            string domain;
            for(int i=0;i<static_cast<int>(icjson.at("domain").size());i++){
                if(!icjson.at("domain").at(i).is_string()){
                    MessagePrinter::printErrorTxt("domain name is invalid in 'domain' of '"+icBlock.m_ICBlockName+"' in your 'ics' subblock, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                domain=icjson.at("domain").at(i);
                if(domain.size()<1){
                    MessagePrinter::printErrorTxt("domain name is invalid or empty in 'domain' of '"+icBlock.m_ICBlockName+"' in your 'ics' subblock, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                if(!t_fecell.isFECellBulkElmtPhyNameValid(domain)){
                   MessagePrinter::printErrorTxt("domain name is invalid in 'domain' of '"+icBlock.m_ICBlockName+"' in your 'ics' subblock, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                icBlock.m_DomainNameList.push_back(domain);
            }
            HasDomain=true;
        }
        else{
            HasDomain=false;
            MessagePrinter::printErrorTxt("can\'t find 'domain' in '"+icBlock.m_ICBlockName+"' of your 'ics' subblock, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        } //end-of-'domain'-reading

        if(!icjson.contains("parameters")){
            icBlock.m_Params.clear();
        }
        else{
            HasParams=true;
            icBlock.m_Params=icjson.at("parameters");
        }//end-of-parameters-reading

        if(icjson.contains("icvalue")){
            if(!icjson.at("icvalue").is_number()){
                MessagePrinter::printErrorTxt("invalid ic value in '"+icBlock.m_ICBlockName+"' in 'ics' subblock, please check your input file");
                MessagePrinter::exitAsFem();
            }
            icBlock.m_ICValue=icjson.at("icvalue");
        }
        else{
            icBlock.m_ICValue=0.0;
        }//end-of-parameters-reading

        if(!HasParams){
            MessagePrinter::printWarningTxt("no parameters found in '"+icBlock.m_ICBlockName+"', then no parameters will be used in this ic");
        }

        if(!HasType || !HasDomain || !HasDof){
            MessagePrinter::printErrorTxt("information of your 'ics' subblock('"+icBlock.m_ICBlockName
                                          +"') is not complete,please check your input file");
            MessagePrinter::exitAsFem();
        }
        else{
            t_icsystem.addICBlock2List(icBlock);
        }

    }

    return HasType;
}
