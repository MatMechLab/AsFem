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
//+++ Date   : 2022.10.01
//+++ Purpose: read the postprocessor block 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/InputSystem.h"

bool InputSystem::readPostprocessBlock(nlohmann::json &t_json,
                                       const Mesh &t_mesh,
                                       const DofHandler &t_dofhandler,
                                       Postprocessor &t_postprocessor){
    // the json already contains "postprocess" block
    bool HasType=false;
    PostprocessorBlock ppsblock;
    for(auto it=t_json.begin();it!=t_json.end();it++){
        ppsblock.init();

        HasType=false;
        
        ppsblock.m_block_name=it.key();
        nlohmann::json ejson=t_json.at(ppsblock.m_block_name);
        if(ejson.contains("type")){
            if(!ejson.at("type").is_string()){
                MessagePrinter::printErrorTxt("invalid type name in postprocess block("+ppsblock.m_block_name+"), please check your input file");
                MessagePrinter::exitAsFem();
            }
            ppsblock.m_pps_typename=ejson.at("type");
            ppsblock.m_pps_type=PostprocessorType::NULLPPS;
            // for nodal type pps
            if(ppsblock.m_pps_typename=="nodalvalue"){
                ppsblock.m_pps_type=PostprocessorType::NODALVALUE;
            }
            else if(ppsblock.m_pps_typename=="nodalscalarmate"){
                ppsblock.m_pps_type=PostprocessorType::NODALSCALARMATERIALVALUE;
            }
            else if(ppsblock.m_pps_typename=="nodalvectormate"){
                ppsblock.m_pps_type=PostprocessorType::NODALVECTORMATERIALVALUE;
            }
            else if(ppsblock.m_pps_typename=="nodalrank2mate"){
                ppsblock.m_pps_type=PostprocessorType::NODALRANK2MATERIALVALUE;
            }
            else if(ppsblock.m_pps_typename=="nodalrank4mate"){
                ppsblock.m_pps_type=PostprocessorType::NODALRANK4MATERIALVALUE;
            }
            // for side integral type pps
            else if(ppsblock.m_pps_typename=="area"){
                ppsblock.m_pps_type=PostprocessorType::AREA;
            }
            else if(ppsblock.m_pps_typename=="sideaveragevalue"){
                ppsblock.m_pps_type=PostprocessorType::SIDEAVERAGEVALUE;
            }
            else if(ppsblock.m_pps_typename=="sideaveragescalarmate"){
                ppsblock.m_pps_type=PostprocessorType::SIDEAVERAGESCALARMATERIALVALUE;
            }
            else if(ppsblock.m_pps_typename=="sideaveragevectormate"){
                ppsblock.m_pps_type=PostprocessorType::SIDEAVERAGEVECTORMATERIALVALUE;
            }
            else if(ppsblock.m_pps_typename=="sideaveragerank2mate"){
                ppsblock.m_pps_type=PostprocessorType::SIDEAVERAGERANK2MATERIALVALUE;
            }
            else if(ppsblock.m_pps_typename=="sideaveragerank4mate"){
                ppsblock.m_pps_type=PostprocessorType::SIDEAVERAGERANK4MATERIALVALUE;
            }
            else if(ppsblock.m_pps_typename=="sideintegralvalue"){
                ppsblock.m_pps_type=PostprocessorType::SIDEINTEGRATEVALUE;
            }
            else if(ppsblock.m_pps_typename=="sideintegralscalarmate"){
                ppsblock.m_pps_type=PostprocessorType::SIDEINTEGRATESCALARMATERIALVALUE;
            }
            else if(ppsblock.m_pps_typename=="sideintegralvectormate"){
                ppsblock.m_pps_type=PostprocessorType::SIDEINTEGRATEVECTORMATERIALVALUE;
            }
            else if(ppsblock.m_pps_typename=="sideintegralrank2mate"){
                ppsblock.m_pps_type=PostprocessorType::SIDEINTEGRATERANK2MATERIALVALUE;
            }
            else if(ppsblock.m_pps_typename=="sideintegralrank4mate"){
                ppsblock.m_pps_type=PostprocessorType::SIDEINTEGRATERANK4MATERIALVALUE;
            }
            else if(ppsblock.m_pps_typename=="user1sideintegral"){
                ppsblock.m_pps_type=PostprocessorType::USER1SIDEINTEGRALPPS;
            }
            // for volume integral type pps
            else if(ppsblock.m_pps_typename=="volume"){
                ppsblock.m_pps_type=PostprocessorType::VOLUME;
            }
            else if(ppsblock.m_pps_typename=="volumeaveragevalue"){
                ppsblock.m_pps_type=PostprocessorType::VOLUMEAVERAGEVALUE;
            }
            else if(ppsblock.m_pps_typename=="volumeaveragescalarmate"){
                ppsblock.m_pps_type=PostprocessorType::VOLUMEAVERAGESCALARMATERIALVALUE;
            }
            else if(ppsblock.m_pps_typename=="volumeaveragevectormate"){
                ppsblock.m_pps_type=PostprocessorType::VOLUMEAVERAGEVECTORMATERIALVALUE;
            }
            else if(ppsblock.m_pps_typename=="volumeaveragerank2mate"){
                ppsblock.m_pps_type=PostprocessorType::VOLUMEAVERAGERANK2MATERIALVALUE;
            }
            else if(ppsblock.m_pps_typename=="volumeaveragerank4mate"){
                ppsblock.m_pps_type=PostprocessorType::VOLUMEAVERAGERANK4MATERIALVALUE;
            }
            else if(ppsblock.m_pps_typename=="volumeintegralvalue"){
                ppsblock.m_pps_type=PostprocessorType::VOLUMEINTEGRATEVALUE;
            }
            else if(ppsblock.m_pps_typename=="volumeintegralscalarmate"){
                ppsblock.m_pps_type=PostprocessorType::VOLUMEINTEGRATESCALARMATERIALVALUE;
            }
            else if(ppsblock.m_pps_typename=="volumeintegralvectormate"){
                ppsblock.m_pps_type=PostprocessorType::VOLUMEINTEGRATEVECTORMATERIALVALUE;
            }
            else if(ppsblock.m_pps_typename=="volumeintegralrank2mate"){
                ppsblock.m_pps_type=PostprocessorType::VOLUMEINTEGRATERANK2MATERIALVALUE;
            }
            else if(ppsblock.m_pps_typename=="volumeintegralrank4mate"){
                ppsblock.m_pps_type=PostprocessorType::VOLUMEINTEGRATERANK4MATERIALVALUE;
            }
            else if(ppsblock.m_pps_typename=="user1volumeintegral"){
                ppsblock.m_pps_type=PostprocessorType::USER1VOLUMEINTEGRALPPS;
            }
            else{
                MessagePrinter::printErrorTxt("unsupported element type in "+ppsblock.m_block_name+", please check your input file");
                MessagePrinter::exitAsFem();
            }//end-of-element-type-name
            HasType=true;
        }
        else{
            MessagePrinter::printErrorTxt("can\'t find type name in postprocess block("+ppsblock.m_block_name+"), please check your input file");
            MessagePrinter::exitAsFem();
        }//end-of-'type'-reading

        if(ejson.contains("dof")){
            if(ejson.at("dof").size()<1){
                MessagePrinter::printErrorTxt("invalid dofs name or empty dofs name in postprocess block("+ppsblock.m_block_name+"), please check your input file");
                MessagePrinter::exitAsFem();
            }
            else{
                string dofname;
                if(!ejson.at("dof").is_string()){
                    MessagePrinter::printErrorTxt("dof's name is invalid in 'dof' of the postprocess block("+ppsblock.m_block_name+"), please check your input file");
                    MessagePrinter::exitAsFem();
                }
                dofname=ejson.at("dof");
                if(dofname.size()<1){
                    MessagePrinter::printErrorTxt("dof's name is invalid or empty in 'dof' of the postprocess block("+ppsblock.m_block_name+"), please check your input file");
                    MessagePrinter::exitAsFem();
                }
                if(!t_dofhandler.isValidDofName(dofname)){
                    MessagePrinter::printErrorTxt("dof's name is invalid in 'dof' of the postprocess block("+ppsblock.m_block_name+"), it must be one of the names in your 'dofs' block, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                ppsblock.m_dofname=dofname;
                ppsblock.m_dofid=t_dofhandler.getDofIDViaName(dofname);
            }
        }

        if(ejson.contains("domain")){
            if(ejson.at("domain").size()<1){
                MessagePrinter::printErrorTxt("invalid or empty domain name in ["+ppsblock.m_block_name+"] of your [postprocess] subblock, please check your input file");
                MessagePrinter::exitAsFem();
            }
            string domain;
            for(int i=0;i<static_cast<int>(ejson.at("domain").size());i++){
                if(!ejson.at("domain").at(i).is_string()){
                    MessagePrinter::printErrorTxt("domain name is invalid in 'domain' of ["+ppsblock.m_block_name+"] in your [postprocess] subblock, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                domain=ejson.at("domain").at(i);
                if(domain.size()<1){
                    MessagePrinter::printErrorTxt("domain name is invalid or empty in 'domain' of ["+ppsblock.m_block_name+"] in your [postprocess] subblock, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                if(!t_mesh.isBulkElmtPhyNameValid(domain)){
                   MessagePrinter::printErrorTxt("domain name is invalid in 'domain' of ["+ppsblock.m_block_name+"] in your [postprocess] subblock, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                ppsblock.m_domainnamelist.push_back(domain);
            }
        }
        else{
            ppsblock.m_domainnamelist.clear();
            ppsblock.m_domainnamelist.push_back("alldomain");
        } //end-of-'domain'-reading

        if(ejson.contains("side")){
            if(ejson.at("side").size()<1){
                MessagePrinter::printErrorTxt("invalid or empty side name in ["+ppsblock.m_block_name+"] of your [postprocess] subblock, please check your input file");
                MessagePrinter::exitAsFem();
            }
            string sidename;
            for(int i=0;i<static_cast<int>(ejson.at("side").size());i++){
                if(!ejson.at("side").at(i).is_string()){
                    MessagePrinter::printErrorTxt("side name is invalid in 'side' of ["+ppsblock.m_block_name+"] in your [postprocess] subblock, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                sidename=ejson.at("side").at(i);
                if(sidename.size()<1){
                    MessagePrinter::printErrorTxt("side name is invalid or empty in 'side' of ["+ppsblock.m_block_name+"] in your [postprocess] subblock, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                if(!t_mesh.isBCElmtPhyNameValid(sidename)){
                   MessagePrinter::printErrorTxt("side name is invalid in 'side' of ["+ppsblock.m_block_name+"] in your [postprocess] subblock, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                ppsblock.m_sidenamelist.push_back(sidename);
            }
        }
        else{
            ppsblock.m_sidenamelist.clear();
        } //end-of-'domain'-reading

        if(ejson.contains("parameters")){
            // we will keep the parameter's json into m_json_params, this will be used later
            // for material property calculation (UMAT)
            ppsblock.m_parameters=ejson.at("parameters");
        }
        
        if(!HasType){
            MessagePrinter::printErrorTxt("information in ["+ppsblock.m_block_name
                                         +"] of your [postprocess] subblock is not complete, please check your input file");
            return false;
        }
        t_postprocessor.addPPSBlock2List(ppsblock);
    }// end-of-iter-of-'pps'-block

    return HasType;
}