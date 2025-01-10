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
//+++ Purpose: Define the initial condition system in AsFem
//+++          Here we can apply:
//+++               1) constant value ic
//+++               2) random value ic
//+++               3) other type or user defined ic
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ICSystem/ICSystem.h"

ICSystem::ICSystem(){
    m_icblock_list.clear();
    m_icblocks_num=0;
    m_localU.clean();
}

void ICSystem::init(const int &nodal_dofs){
    m_localU.resize(1+nodal_dofs,0.0);
}

void ICSystem::releaseMemory(){
    m_icblock_list.clear();
    m_icblocks_num=0;
    m_localU.clean();
}

void ICSystem::addICBlock2List(const ICBlock &icblock){
    if(m_icblock_list.size()<1){
        m_icblock_list.push_back(icblock);
        m_icblocks_num=1;
    }
    else{
        bool NotInList=true;
        for(const auto &it:m_icblock_list){
            if(it.m_ICBlockName==icblock.m_ICBlockName){
                NotInList=false;break;
            }
        }
        if(NotInList){
            m_icblock_list.push_back(icblock);
            m_icblocks_num+=1;
        }
        else{
            MessagePrinter::printErrorTxt("duplicate ["+icblock.m_ICBlockName+"] in your [ics] sub block,"
                                          " please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
}

void ICSystem::printICSystemInfo()const{
    MessagePrinter::printNormalTxt("Initial condition system information summary");
    for(const auto &it:m_icblock_list){
        it.printICBlockInfo();
    }
    MessagePrinter::printStars();
}