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
//+++ Date   : 2022.05.11
//+++ Purpose: the element block for input file reading
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "nlohmann/json.hpp"

#include "ElmtSystem/ElmtType.h"
#include "MateSystem/MateType.h"
#include "Utils/MessagePrinter.h"

/**
 * This class defines the basic information for one single element block in input file
 */
class ElmtBlock{
public:
    /**
     * constructor
     */
    ElmtBlock(){
        m_elmtblock_index=0;

        m_elmt_blockname.clear();
        m_elmt_typename.clear();
        m_elmttype=ElmtType::NULLELMT;

        m_dof_names.clear();
        m_dof_ids.clear();

        m_domain_namelist.clear();

        m_mate_typename.clear();
        m_matetype=MateType::NULLMATE;

        m_json_params.clear();
    }
    /**
     * reset the content of current element block
     */
    void reset(){
        m_elmtblock_index=0;

        m_elmt_blockname.clear();
        m_elmt_typename.clear();
        m_elmttype=ElmtType::NULLELMT;

        m_dof_names.clear();
        m_dof_ids.clear();

        m_domain_namelist.clear();

        m_mate_typename.clear();
        m_matetype=MateType::NULLMATE;

        m_json_params.clear();
    }
    /**
     * print out the information of current element block
     */
    void printElmtBlockInfo()const{
        MessagePrinter::printDashLine();
        MessagePrinter::printNormalTxt(" bulk element block-"+to_string(m_elmtblock_index)+" info");
        MessagePrinter::printNormalTxt("  block name = "+m_elmt_blockname+", type name = "+m_elmt_typename);
        string str;

        str="";
        for(const auto &it:m_dof_names) str+=it+" ";
        MessagePrinter::printNormalTxt("  dofs name = "+str);

        str="";
        for(const auto &it:m_dof_ids) str+=to_string(it)+" ";
        MessagePrinter::printNormalTxt("  dofs id = "+str);

        MessagePrinter::printNormalTxt("  material type name = "+m_mate_typename);

        str="";
        for(const auto &it:m_domain_namelist) str+=it+" ";
        MessagePrinter::printNormalTxt("  domain = "+str);

        if(m_json_params.size()>0){
            MessagePrinter::printNormalTxt("  parameters are:");
            char buff[69];
            for(auto it=m_json_params.begin();it!=m_json_params.end();it++){
                if(it.value().is_boolean()){
                    if(it.value()==true){
                        snprintf(buff,69,"    %16s =    true",it.key().c_str());
                    }
                    else{
                        snprintf(buff,69,"    %16s =    false",it.key().c_str());
                    }
                }
                else if(it.value().is_string()){
                    string substr=it.value();
                    snprintf(buff,69,"    %16s = %16s",it.key().c_str(),substr.c_str());
                }
                else if(it.value().is_array()){
                    if(it.value().size()<3){
                        MessagePrinter::printErrorTxt("Invalid vector size in "+it.key()+" of your parameters in "+m_elmt_blockname+", please check your input file");
                        MessagePrinter::exitAsFem();
                    }
                    snprintf(buff,69,"    %15s = %13.5e %13.5e %13.5e",it.key().c_str(),static_cast<double>(it.value().at(0)),
                                                                                        static_cast<double>(it.value().at(1)),
                                                                                        static_cast<double>(it.value().at(2)));
                }
                else{
                    snprintf(buff,69,"    %16s = %14.5e",it.key().c_str(),static_cast<double>(it.value()));
                }
                str=buff;
                MessagePrinter::printNormalTxt(str);
            }
        }
    }
public:
    int m_elmtblock_index;/**< the index for the order of current element block in the input file, start from 1 */
    string m_elmt_blockname;/**< the string name of current element block(single block) */
    string m_elmt_typename;/**< the string name of current element block */
    ElmtType m_elmttype;/**< the type of current element */
    vector<string> m_dof_names;/**< the string name list of the dofs used in current element */
    vector<int> m_dof_ids;/**< the dof ids related to the dof names */
    vector<string> m_domain_namelist;/**< the physical name vector of the domain for current element */
    //*** for materials
    string m_mate_typename;/**< string name for material type of current element */
    MateType m_matetype;/**< the type of material used in current element */
    nlohmann::json m_json_params;/**< json class for material paramters of current element */

};