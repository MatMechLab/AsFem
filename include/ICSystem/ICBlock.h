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
//+++ Date   : 2020.07.11
//+++ Purpose: Define [ics] sub block for our input file
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "nlohmann/json.hpp"
#include "ICSystem/ICType.h"
#include "Utils/MessagePrinter.h"

using std::vector;
using std::string;

/**
 * This class defines the basic information for one single ic block in the input file
 */
class ICBlock{
public:
    /**
     * constructor
     */
    ICBlock(){
        m_icBlockIndex=0;
        m_icBlockName.clear();
        m_icTypeName.clear();
        m_icType=ICType::NULLIC;
        m_dofsName.clear();
        m_dofIDs.clear();
        m_json_params.clear();
        m_icvalue=0.0;
        m_domainNameList.clear();
    }
    /**
     * reset the bc block info
     */
    void reset(){
        m_icBlockIndex=0;
        m_icBlockName.clear();
        m_icTypeName.clear();
        m_icType=ICType::NULLIC;
        m_dofsName.clear();
        m_dofIDs.clear();
        m_json_params.clear();
        m_icvalue=0.0;
        m_domainNameList.clear();
    }

    /**
     * print out the information of current ic block
     */
    void printICBlockInfo()const{
        MessagePrinter::printNormalTxt(" ics block-"+to_string(m_icBlockIndex)+" info");
        MessagePrinter::printNormalTxt("  block name = "+m_icBlockName+", type name = "+m_icTypeName);
        string str;

        str="";
        for(const auto &it:m_dofsName) str+=it+" ";
        MessagePrinter::printNormalTxt("  dofs name = "+str);

        str="";
        for(const auto &it:m_dofIDs) str+=to_string(it)+" ";
        MessagePrinter::printNormalTxt("  dofs id = "+str);

        str="";
        for(const auto &it:m_domainNameList){
                str+=it+" ";
        }
        MessagePrinter::printNormalTxt("  domain = "+str);
        
        MessagePrinter::printNormalTxt("  ic value = "+to_string(m_icvalue));
        if(m_json_params.size()>0){
            char buff[69];
            MessagePrinter::printNormalTxt("  parameters are:");
            for(auto it=m_json_params.begin();it!=m_json_params.end();it++){
                snprintf(buff,69,"    %16s = %14.5e",it.key().c_str(),static_cast<double>(it.value()));
                str=buff;
                MessagePrinter::printNormalTxt(str);
            }
        }
        MessagePrinter::printDashLine();
    }

public:
    int            m_icBlockIndex;/**< index of current ic block */
    string         m_icBlockName;/**< for the name of current ic block */
    string         m_icTypeName;/**< the for ic type name of current ic block */
    ICType         m_icType;/**< ic type of current block */
    vector<string> m_dofsName;/**< list of dofs */
    vector<int>    m_dofIDs;/**< list of dofs id */
    double         m_icvalue;/**< the initial condition value */
    nlohmann::json m_json_params;/**< json class for material paramters of current ic block */
    vector<string> m_domainNameList;/**< domain name list */
    
};