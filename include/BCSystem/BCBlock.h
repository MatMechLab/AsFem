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
        m_bcBlockIndex=0;
        m_bcBlockName.clear();
        m_bcTypeName.clear();
        m_bcType=BCType::NULLBC;
        m_dofsName.clear();
        m_dofIDs.clear();
        m_bcValue=0.0;
        m_json_params.clear();
        m_boundaryNameList.clear();
        m_isTimeDependent=false;
    }
    /**
     * reset the bc block info
     */
    void reset(){
        m_bcBlockIndex=0;
        m_bcBlockName.clear();
        m_bcTypeName.clear();
        m_bcType=BCType::NULLBC;
        m_dofsName.clear();
        m_dofIDs.clear();
        m_bcValue=0.0;
        m_json_params.clear();
        m_boundaryNameList.clear();
        m_isTimeDependent=false;
    }

    /**
     * print out the information of current element block
     */
    void printBCBlockInfo()const{
        MessagePrinter::printDashLine();
        MessagePrinter::printNormalTxt(" bcs block-"+to_string(m_bcBlockIndex)+" info");
        MessagePrinter::printNormalTxt("  block name = "+m_bcBlockName+", type name = "+m_bcTypeName);
        string str;

        str="";
        for(const auto &it:m_dofsName) str+=it+" ";
        MessagePrinter::printNormalTxt("  dofs name = "+str);

        str="";
        for(const auto &it:m_dofIDs) str+=to_string(it)+" ";
        MessagePrinter::printNormalTxt("  dofs id = "+str);

        str="";
        for(const auto &it:m_boundaryNameList){
                str+=it+" ";
        }
        MessagePrinter::printNormalTxt("  boundary = "+str);

        char buff[69];
        if(m_isTimeDependent){
            snprintf(buff,69,"  bc value = %14.5e (time dependent)",m_bcValue);
            str=buff;
        }
        else{
            snprintf(buff,69,"  bc value = %14.5e",m_bcValue);
            str=buff;
        }
        MessagePrinter::printNormalTxt(str);
        
        if(m_json_params.size()>0){
            MessagePrinter::printNormalTxt("  parameters are:");
            for(auto it=m_json_params.begin();it!=m_json_params.end();it++){
                string parname=it.key();
                if(m_json_params.at(parname).is_array()){
                    if(m_json_params.at(parname).size()<3){
                        MessagePrinter::printErrorTxt("Invalid vector size in "+it.key()+" of "+m_bcBlockName+", please check your input file");
                        MessagePrinter::exitAsFem();
                    }
                    snprintf(buff,69,"    %15s = %13.5e %13.5e %13.5e",it.key().c_str(),static_cast<double>(m_json_params.at(parname).at(0)),
                                                                                        static_cast<double>(m_json_params.at(parname).at(1)),
                                                                                        static_cast<double>(m_json_params.at(parname).at(2)));
                    str=buff;
                }
                else if(m_json_params.at(parname).is_boolean()){
                    if(m_json_params.at(parname)){
                        snprintf(buff,69,"    %15s =   true",it.key().c_str());
                    }
                    else{
                        snprintf(buff,69,"    %15s =   false",it.key().c_str());
                    }
                    str=buff;
                }
                else if(m_json_params.at(parname).is_string()){
                    string txt=it.value();
                    snprintf(buff,69,"    %15s =   %-20s",it.key().c_str(),txt.c_str());
                    str=buff;
                }
                else if(m_json_params.at(parname).is_number()||
                        m_json_params.at(parname).is_number_float()||
                        m_json_params.at(parname).is_number_integer()||
                        m_json_params.at(parname).is_number_unsigned()){
                    snprintf(buff,69,"    %15s = %13.5e",it.key().c_str(),static_cast<double>(it.value()));
                    str=buff;
                }
                else{
                    MessagePrinter::printErrorTxt("Unknown or unsupported options in parameters of "+m_bcBlockName+", please check your input file");
                    MessagePrinter::exitAsFem();
                }
                MessagePrinter::printNormalTxt(str);
            }
        }
    }

public:
    int            m_bcBlockIndex;/**< index of current bc block */
    string         m_bcBlockName;/**< for the name of current bc block */
    string         m_bcTypeName;/**< the for bc type name of current bc block */
    BCType         m_bcType;/**< bc type of current block */
    vector<string> m_dofsName;/**< list of dofs */
    vector<int>    m_dofIDs;/**< list of dofs id */
    double         m_bcValue;/**< bc value */
    nlohmann::json m_json_params;/**< json class for material paramters of current bc block */
    vector<string> m_boundaryNameList;/**< it could be either an element set or a node set */
    bool           m_isTimeDependent;/**< boolen for time dependent status */
    
};