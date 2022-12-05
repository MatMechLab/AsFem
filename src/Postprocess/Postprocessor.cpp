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
//+++ Purpose: implement postprocess system in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/Postprocessor.h"

Postprocessor::Postprocessor(){
    m_output_interval=1;
    m_pps_namelist.clear();
    m_pps_values.clear();

    m_inputfilename.clear();

    m_csv_filename.clear();

    m_pps_blocklist.clear();
    m_pps_blocksnum=0;

}
void Postprocessor::init(){
    m_nodes0.resize(27+1);
    m_nodes.resize(27+1);
}
void Postprocessor::releaseMemory(){
    m_nodes0.clear();
    m_nodes.clear();
}
//**********************************************
void Postprocessor::addPPSBlock2List(const PostprocessorBlock &t_block){
    if(m_pps_blocklist.size()<1){
        m_pps_blocklist.push_back(t_block);
        m_pps_blocksnum=1;
        m_pps_namelist.push_back(t_block.m_block_name);
        m_pps_values.push_back(0.0);
    }
    else{
        bool IsInList=false;
        for(const auto &it:m_pps_blocklist){
            if(it.m_block_name==t_block.m_block_name){
                IsInList=true;
                break;
            }
        }
        if(!IsInList){
            m_pps_blocklist.push_back(t_block);
            m_pps_blocksnum+=1;
            m_pps_namelist.push_back(t_block.m_block_name);
            m_pps_values.push_back(0.0);
        }
        else{
            MessagePrinter::printErrorTxt("duplicate pps block(name="+t_block.m_block_name+"), please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
}
//**********************************************
void Postprocessor::printInfo()const{
    MessagePrinter::printNormalTxt("Postprocess information summary");
    MessagePrinter::printNormalTxt("  number of postprocess block="+to_string(getPPSBlocksNum()));
    if(getPPSBlocksNum()){
        string str;
        str="  postprocess variables= ";
        for(auto it:m_pps_namelist){
            str+=it+" ";
        }
        MessagePrinter::printNormalTxt(str);
        for(const auto &it:m_pps_blocklist) it.printBlockInfo();
    }
    MessagePrinter::printStars();
}