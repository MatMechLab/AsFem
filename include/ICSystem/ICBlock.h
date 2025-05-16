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
        m_ICBlockIndex=0;
        m_ICBlockName.clear();
        m_ICTypeName.clear();
        m_ICType=ICType::NULLIC;
        m_DofNames.clear();
        m_DofIDs.clear();
        m_Params.clear();
        m_ICValue=0.0;
        m_DomainNameList.clear();
    }
    /**
     * reset the bc block info
     */
    void reset(){
        m_ICBlockIndex=0;
        m_ICBlockName.clear();
        m_ICTypeName.clear();
        m_ICType=ICType::NULLIC;
        m_DofNames.clear();
        m_DofIDs.clear();
        m_Params.clear();
        m_ICValue=0.0;
        m_DomainNameList.clear();
    }

    /**
     * print out the information of current ic block
     */
    void printICBlockInfo()const{
        MessagePrinter::printNormalTxt(" ics block-"+to_string(m_ICBlockIndex)+" info");
        MessagePrinter::printNormalTxt("  block name = "+m_ICBlockName+", type name = "+m_ICTypeName);
        string str;

        str="";
        for(const auto &it:m_DofNames) str+=it+" ";
        MessagePrinter::printNormalTxt("  dofs name = "+str);

        str="";
        for(const auto &it:m_DofIDs) str+=to_string(it)+" ";
        MessagePrinter::printNormalTxt("  dofs id = "+str);

        str="";
        for(const auto &it:m_DomainNameList){
                str+=it+" ";
        }
        MessagePrinter::printNormalTxt("  domain = "+str);
        
        MessagePrinter::printNormalTxt("  ic value = "+to_string(m_ICValue));
        if(m_Params.size()>0){
            char buff[69];
            MessagePrinter::printNormalTxt("  parameters are:");
            for(auto it=m_Params.begin();it!=m_Params.end();it++){
                snprintf(buff,69,"    %16s = %14.5e",it.key().c_str(),static_cast<double>(it.value()));
                str=buff;
                MessagePrinter::printNormalTxt(str);
            }
        }
        MessagePrinter::printDashLine();
    }

public:
    int            m_ICBlockIndex;/**< index of current ic block */
    string         m_ICBlockName;/**< for the name of current ic block */
    string         m_ICTypeName;/**< the for ic type name of current ic block */
    ICType         m_ICType;/**< ic type of current block */
    vector<string> m_DofNames;/**< list of dofs */
    vector<int>    m_DofIDs;/**< list of dofs id */
    double         m_ICValue;/**< the initial condition value */
    nlohmann::json m_Params;/**< json class for material paramters of current ic block */
    vector<string> m_DomainNameList;/**< domain name list */
    
};