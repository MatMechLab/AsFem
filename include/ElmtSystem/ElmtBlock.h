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
        m_ElmtBlockIndex=0;

        m_ElmtBlockName.clear();
        m_ElmtTypeName.clear();
        m_ElmtType=ElmtType::NULLELMT;

        m_DofNames.clear();
        m_DofIDs.clear();

        m_DomainNameList.clear();

        m_MateTypeName.clear();
        m_MateType=MateType::NULLMATE;

        m_JsonParams.clear();
    }
    /**
     * reset the content of current element block
     */
    void reset(){
        m_ElmtBlockIndex=0;

        m_ElmtBlockName.clear();
        m_ElmtTypeName.clear();
        m_ElmtType=ElmtType::NULLELMT;

        m_DofNames.clear();
        m_DofIDs.clear();

        m_DomainNameList.clear();

        m_MateTypeName.clear();
        m_MateType=MateType::NULLMATE;

        m_JsonParams.clear();
    }
    /**
     * print out the information of current element block
     */
    void printElmtBlockInfo()const{
        MessagePrinter::printDashLine();
        MessagePrinter::printNormalTxt(" bulk element block-"+to_string(m_ElmtBlockIndex)+" info");
        MessagePrinter::printNormalTxt("  block name = "+m_ElmtBlockName+", type name = "+m_ElmtTypeName);
        string str;

        str="";
        for(const auto &it:m_DofNames) str+=it+" ";
        MessagePrinter::printNormalTxt("  dofs name = "+str);

        str="";
        for(const auto &it:m_DofIDs) str+=to_string(it)+" ";
        MessagePrinter::printNormalTxt("  dofs id = "+str);

        MessagePrinter::printNormalTxt("  material type name = "+m_MateTypeName);

        str="";
        for(const auto &it:m_DomainNameList) str+=it+" ";
        MessagePrinter::printNormalTxt("  domain = "+str);

        if(m_JsonParams.size()>0){
            MessagePrinter::printNormalTxt("  parameters are:");
            char buff[69];
            for(auto it=m_JsonParams.begin();it!=m_JsonParams.end();it++){
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
                        MessagePrinter::printErrorTxt("Invalid vector size in "+it.key()+" of your parameters in "+m_ElmtBlockName+", please check your input file");
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
    int m_ElmtBlockIndex;/**< the index for the order of current element block in the input file, start from 1 */
    string m_ElmtBlockName;/**< the string name of current element block(single block) */
    string m_ElmtTypeName;/**< the string name of current element block */
    ElmtType m_ElmtType;/**< the type of current element */
    vector<string> m_DofNames;/**< the string name list of the dofs used in current element */
    vector<int> m_DofIDs;/**< the dof ids related to the dof names */
    vector<string> m_DomainNameList;/**< the physical name vector of the domain for current element */
    //*** for materials
    string m_MateTypeName;/**< string name for material type of current element */
    MateType m_MateType;/**< the type of material used in current element */
    nlohmann::json m_JsonParams;/**< json class for material paramters of current element */

};