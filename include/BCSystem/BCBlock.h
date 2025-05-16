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
//+++ Date   : 2020.07.10
//+++ Purpose: Define [bcs] sub block for our input file
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "nlohmann/json.hpp"
#include "BCSystem/BCType.h"
#include "Utils/MessagePrinter.h"

using namespace std;

/**
 * This class defines the basic information for one single bc block in the input file
 */
class BCBlock{
public:
    /**
     * constructor
     */
    BCBlock(){
        m_BCBlockIndex=0;
        m_BCBlockName.clear();
        m_BCTypeName.clear();
        m_BCType=BCType::NULLBC;
        m_DofsName.clear();
        m_DofIDs.clear();
        m_BCValue=0.0;
        m_JsonParams.clear();
        m_BoundaryNameList.clear();
        m_IsTimeDependent=false;
    }
    /**
     * reset the bc block info
     */
    void reset(){
        m_BCBlockIndex=0;
        m_BCBlockName.clear();
        m_BCTypeName.clear();
        m_BCType=BCType::NULLBC;
        m_DofsName.clear();
        m_DofIDs.clear();
        m_BCValue=0.0;
        m_JsonParams.clear();
        m_BoundaryNameList.clear();
        m_IsTimeDependent=false;
    }

    /**
     * print out the information of current element block
     */
    void printBCBlockInfo()const{
        MessagePrinter::printDashLine();
        MessagePrinter::printNormalTxt(" bcs block-"+to_string(m_BCBlockIndex)+" info");
        MessagePrinter::printNormalTxt("  block name = "+m_BCBlockName+", type name = "+m_BCTypeName);
        string str;

        str="";
        for(const auto &it:m_DofsName) str+=it+" ";
        MessagePrinter::printNormalTxt("  dofs name = "+str);

        str="";
        for(const auto &it:m_DofIDs) str+=to_string(it)+" ";
        MessagePrinter::printNormalTxt("  dofs id = "+str);

        str="";
        for(const auto &it:m_BoundaryNameList){
                str+=it+" ";
        }
        MessagePrinter::printNormalTxt("  boundary = "+str);

        char buff[69];
        if(m_IsTimeDependent){
            snprintf(buff,69,"  bc value = %14.5e (time dependent)",m_BCValue);
            str=buff;
        }
        else{
            snprintf(buff,69,"  bc value = %14.5e",m_BCValue);
            str=buff;
        }
        MessagePrinter::printNormalTxt(str);
        
        if(m_JsonParams.size()>0){
            MessagePrinter::printNormalTxt("  parameters are:");
            for(auto it=m_JsonParams.begin();it!=m_JsonParams.end();it++){
                string parname=it.key();
                if(m_JsonParams.at(parname).is_array()){
                    if(m_JsonParams.at(parname).size()<3){
                        MessagePrinter::printErrorTxt("Invalid vector size in "+it.key()+" of "+m_BCBlockName+", please check your input file");
                        MessagePrinter::exitAsFem();
                    }
                    snprintf(buff,69,"    %15s = %13.5e %13.5e %13.5e",it.key().c_str(),static_cast<double>(m_JsonParams.at(parname).at(0)),
                                                                                        static_cast<double>(m_JsonParams.at(parname).at(1)),
                                                                                        static_cast<double>(m_JsonParams.at(parname).at(2)));
                    str=buff;
                }
                else if(m_JsonParams.at(parname).is_boolean()){
                    if(m_JsonParams.at(parname)){
                        snprintf(buff,69,"    %15s =   true",it.key().c_str());
                    }
                    else{
                        snprintf(buff,69,"    %15s =   false",it.key().c_str());
                    }
                    str=buff;
                }
                else if(m_JsonParams.at(parname).is_string()){
                    string txt=it.value();
                    snprintf(buff,69,"    %15s =   %-20s",it.key().c_str(),txt.c_str());
                    str=buff;
                }
                else if(m_JsonParams.at(parname).is_number()||
                        m_JsonParams.at(parname).is_number_float()||
                        m_JsonParams.at(parname).is_number_integer()||
                        m_JsonParams.at(parname).is_number_unsigned()){
                    snprintf(buff,69,"    %15s = %13.5e",it.key().c_str(),static_cast<double>(it.value()));
                    str=buff;
                }
                else{
                    MessagePrinter::printErrorTxt("Unknown or unsupported options in parameters of "+m_BCBlockName+", please check your input file");
                    MessagePrinter::exitAsFem();
                }
                MessagePrinter::printNormalTxt(str);
            }
        }
    }

public:
    int            m_BCBlockIndex;/**< index of current bc block */
    string         m_BCBlockName;/**< for the name of current bc block */
    string         m_BCTypeName;/**< the for bc type name of current bc block */
    BCType         m_BCType;/**< bc type of current block */
    vector<string> m_DofsName;/**< list of dofs */
    vector<int>    m_DofIDs;/**< list of dofs id */
    double         m_BCValue;/**< bc value */
    nlohmann::json m_JsonParams;/**< json class for material paramters of current bc block */
    vector<string> m_BoundaryNameList;/**< it could be either an element set or a node set */
    bool           m_IsTimeDependent;/**< boolen for time dependent status */
    
};