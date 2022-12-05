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

bool InputSystem::readICsBlock(nlohmann::json &t_json,const Mesh &t_mesh,const DofHandler &t_dofhandler,ICSystem &t_icsystem){
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

        icBlock.m_icBlockIndex=nblocks;
        icBlock.m_icBlockName=it.key();
        nlohmann::json icjson=t_json.at(icBlock.m_icBlockName);//it.value();// now ejson gose into 'ic-x'

        if(icjson.contains("type")){
            if(!icjson.at("type").is_string()){
                MessagePrinter::printErrorTxt("invalid type name in ["+icBlock.m_icBlockName+"] in [ics] subblock, please check your input file");
                MessagePrinter::exitAsFem();
            }
            icBlock.m_icTypeName=icjson.at("type");
            icBlock.m_icType=ICType::NULLIC;
            if(icBlock.m_icTypeName=="const"){
                icBlock.m_icType=ICType::CONSTIC;
            }
            else if(icBlock.m_icTypeName=="random"){
                icBlock.m_icType=ICType::RANDOMIC;
            }
            else if(icBlock.m_icTypeName=="rectangle"){
                icBlock.m_icType=ICType::RECTANGLEIC;
            }
            else if(icBlock.m_icTypeName=="cubic"){
                icBlock.m_icType=ICType::CUBICIC;
            }
            else if(icBlock.m_icTypeName=="circle"){
                icBlock.m_icType=ICType::CIRCLEIC;
            }
            else if(icBlock.m_icTypeName=="smoothcircle"){
                icBlock.m_icType=ICType::SMOOTHCIRCLEIC;
            }
            else if(icBlock.m_icTypeName=="spherical"){
                icBlock.m_icType=ICType::SPHERICALIC;
            }
            // for user-defined ics
            else if(icBlock.m_icTypeName=="user1"){
                icBlock.m_icType=ICType::USER1IC;
            }
            else if(icBlock.m_icTypeName=="user2"){
                icBlock.m_icType=ICType::USER2IC;
            }
            else if(icBlock.m_icTypeName=="user3"){
                icBlock.m_icType=ICType::USER3IC;
            }
            else if(icBlock.m_icTypeName=="user4"){
                icBlock.m_icType=ICType::USER4IC;
            }
            else if(icBlock.m_icTypeName=="user5"){
                icBlock.m_icType=ICType::USER5IC;
            }
            else if(icBlock.m_icTypeName=="user6"){
                icBlock.m_icType=ICType::USER6IC;
            }
            else if(icBlock.m_icTypeName=="user7"){
                icBlock.m_icType=ICType::USER7IC;
            }
            else if(icBlock.m_icTypeName=="user8"){
                icBlock.m_icType=ICType::USER8IC;
            }
            else if(icBlock.m_icTypeName=="user9"){
                icBlock.m_icType=ICType::USER9IC;
            }
            else if(icBlock.m_icTypeName=="user10"){
                icBlock.m_icType=ICType::USER10IC;
            }
            else{
                MessagePrinter::printErrorTxt("type="+icBlock.m_icTypeName
                    +" is invalid in ["+icBlock.m_icBlockName+"] of your [ics] subblock, please check your input file");
                MessagePrinter::exitAsFem();
            }
            HasType=true;
        }//end-of-type-reading

        if(icjson.contains("dofs")){
            if(icjson.at("dofs").size()<1){
                MessagePrinter::printErrorTxt("invalid dofs name or empty dofs name in ["
                                              +icBlock.m_icBlockName+"] of your [ics] subblock, please check your input file");
                MessagePrinter::exitAsFem();
            }
            else{
                string dofname;
                for(int i=0;i<static_cast<int>(icjson.at("dofs").size());i++){
                    if(!icjson.at("dofs").at(i).is_string()){
                        MessagePrinter::printErrorTxt("dof-"+to_string(i+1)+"'s name is invalid in 'dofs' of ["+icBlock.m_icBlockName+"] in your [bcs] subblock, please check your input file");
                        MessagePrinter::exitAsFem();
                    }
                    dofname=icjson.at("dofs").at(i);
                    if(dofname.size()<1){
                        MessagePrinter::printErrorTxt("dof-"+to_string(i+1)+"'s name is invalid or empty in 'dofs' of ["+icBlock.m_icBlockName+"] in your [bcs] subblock, please check your input file");
                        MessagePrinter::exitAsFem();
                    }
                    if(!t_dofhandler.isValidDofName(dofname)){
                        MessagePrinter::printErrorTxt("dof-"+to_string(i+1)+"'s name is invalid in 'dofs' of ["+icBlock.m_icBlockName+"] in your [bcs] subblock,"
                                                      " it must be one of the names in your 'dofs' block, please check your input file");
                        MessagePrinter::exitAsFem();
                    }
                    icBlock.m_dofsName.push_back(dofname);
                    icBlock.m_dofIDs.push_back(t_dofhandler.getDofIDViaName(dofname));
                }
                HasDof=true;
            }
        }
        else{
            HasDof=false;
            MessagePrinter::printErrorTxt("can\'t find dofs name in ["+icBlock.m_icBlockName+"] of your [ics] subblock, please check your input file");
            MessagePrinter::exitAsFem();
        }// end-of-dofs-reading

        if(icjson.contains("domain")){
            if(icjson.at("domain").size()<1){
                MessagePrinter::printErrorTxt("invalid or empty domain name in ["+icBlock.m_icBlockName+"] of your [ics] subblock, please check your input file");
                MessagePrinter::exitAsFem();
            }
            string domain;
            for(int i=0;i<static_cast<int>(icjson.at("domain").size());i++){
                if(!icjson.at("domain").at(i).is_string()){
                    MessagePrinter::printErrorTxt("domain name is invalid in 'domain' of ["+icBlock.m_icBlockName+"] in your [ics] subblock, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                domain=icjson.at("domain").at(i);
                if(domain.size()<1){
                    MessagePrinter::printErrorTxt("domain name is invalid or empty in 'domain' of ["+icBlock.m_icBlockName+"] in your [ics] subblock, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                if(!t_mesh.isBulkElmtPhyNameValid(domain)){
                   MessagePrinter::printErrorTxt("domain name is invalid in 'domain' of ["+icBlock.m_icBlockName+"] in your [ics] subblock, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                icBlock.m_domainNameList.push_back(domain);
            }
            HasDomain=true;
        }
        else{
            HasDomain=false;
            MessagePrinter::printErrorTxt("can\'t find 'domain' in ["+icBlock.m_icBlockName+"] of your [ics] subblock, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        } //end-of-'domain'-reading

        if(!icjson.contains("parameters")){
            icBlock.m_json_params.clear();
        }
        else{
            HasParams=true;
            icBlock.m_json_params=icjson.at("parameters");
        }//end-of-parameters-reading

        if(icjson.contains("icvalue")){
            if(!icjson.at("icvalue").is_number()){
                MessagePrinter::printErrorTxt("invalid ic value in ["+icBlock.m_icBlockName+"] in [ics] subblock, please check your input file");
                MessagePrinter::exitAsFem();
            }
            icBlock.m_icvalue=icjson.at("icvalue");
        }
        else{
            icBlock.m_icvalue=0.0;
        }//end-of-parameters-reading

        if(!HasParams){
            MessagePrinter::printWarningTxt("no parameters found in ["+icBlock.m_icBlockName+"], then no parameters will be used in this ic");
        }

        if(!HasType || !HasDomain || !HasDof){
            MessagePrinter::printErrorTxt("information of your [ics] subblock(["+icBlock.m_icBlockName
                                          +"]) is not complete,please check your input file");
            MessagePrinter::exitAsFem();
        }
        else{
            t_icsystem.addICBlock2List(icBlock);
        }

    }

    return HasType;
}
