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
//+++ Date   : 2022.09.22
//+++ Purpose: defines postprocessor block for the input file in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "nlohmann/json.hpp"

#include "Utils/MessagePrinter.h"

using std::vector;
using std::string;

/**
 * This class defines the postprocessor block for the input file reading
 */
class PostprocessorBlock{
public:
    string m_block_name;/** the string name of the current pps block */
    PostprocessorType m_pps_type;/** the pps type of current pps block */
    string m_pps_typename;/** the type name of current pps block */
    vector<string> m_sidenamelist;/** the string name of current side */
    vector<string> m_domainnamelist;/** the string name of current domain */
    string m_dofname;/**< the string name of current dof */
    int m_dofid;/**< the dof id of current pps block */
    nlohmann::json m_parameters;/**< the parameters taken from input json file */

    /**
     * init all the internal variables
     */
    void init(){
        m_block_name.clear();
        m_pps_type=PostprocessorType::NULLPPS;
        m_pps_typename.clear();
        m_sidenamelist.clear();
        m_domainnamelist.clear();
        m_dofname.clear();
        m_dofid=-1;
        m_parameters.clear();
    }

    void printBlockInfo()const{
        string str;
        MessagePrinter::printDashLine();
        MessagePrinter::printNormalTxt(" postprocess-["+m_block_name+"] info");
        MessagePrinter::printNormalTxt("  type name = "+m_pps_typename);
        if(m_dofid>0){
            MessagePrinter::printNormalTxt("  dof name = "+m_dofname+"(id="+to_string(m_dofid)+")");
        }
        if(m_sidenamelist.size()){
            str="  boundary = ";
            for(const auto &it:m_sidenamelist){
                str+=it+" ";
            }
            MessagePrinter::printNormalTxt(str);
        }
        if(m_domainnamelist.size()){
            str="  domain = ";
            for(auto it:m_domainnamelist) str+=it+" ";
            MessagePrinter::printNormalTxt(str);
        }
        if(m_parameters.size()>0){
            string str;
            char buff[69];
            MessagePrinter::printNormalTxt("  parameters are:");
            for(auto it=m_parameters.begin();it!=m_parameters.end();it++){
                string parname=it.key();
                if(!m_parameters.at(parname).is_array()){
                    if(m_parameters.at(parname).is_number()||
                       m_parameters.at(parname).is_number_float()||
                       m_parameters.at(parname).is_number_integer()||
                       m_parameters.at(parname).is_number_unsigned()){
                        snprintf(buff,69,"    %15s = %13.5e",it.key().c_str(),static_cast<double>(it.value()));
                        str=buff;
                    }
                    else{
                        str=it.value();
                        snprintf(buff,69,"    %15s = %s",it.key().c_str(),str.c_str());
                        str=buff;
                    }
                    MessagePrinter::printNormalTxt(str);
                }
                else{
                    if(m_parameters.at(parname).size()<3){
                        MessagePrinter::printErrorTxt("Invalid vector size in "+it.key()+" of "+m_block_name+", please check your input file");
                        MessagePrinter::exitAsFem();
                    }
                    snprintf(buff,69,"    %15s = %13.5e %13.5e %13.5e",it.key().c_str(),static_cast<double>(m_parameters.at(parname).at(0)),
                                                                                        static_cast<double>(m_parameters.at(parname).at(1)),
                                                                                        static_cast<double>(m_parameters.at(parname).at(2)));
                    str=buff;
                    MessagePrinter::printNormalTxt(str);
                }
            }
        }
    }

};