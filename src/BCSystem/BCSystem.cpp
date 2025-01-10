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
//+++ Purpose: Define the boundary condition system in AsFem
//+++          Here we can apply:
//+++               1) dirichlet bc, i.e. displacement, temperature ...
//+++               2) neuman bc, i.e. flux, force
//+++               3) robin bc as well as user-defined-bc (ubc)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/BCSystem.h"

BCSystem::BCSystem(){
    m_BCBlockList.clear();
    m_BCBlocksNum=0;
    m_DirichletPenalty=1.0e16;
}

void BCSystem::init(const int &dofs){
    m_LocalK.resize(dofs+1,dofs+1,0.0);
    m_LocalR.resize(dofs+1,0.0);
    m_Nodes0.resize(27+1);
    m_Nodes.resize(27+1);

    m_LocalElmtSoln.m_QpU.resize(dofs+1,0.0);
    m_LocalElmtSoln.m_QpUold.resize(dofs+1,0.0);
    m_LocalElmtSoln.m_QpUolder.resize(dofs+1,0.0);
    m_LocalElmtSoln.m_QpV.resize(dofs+1,0.0);

    m_LocalElmtSoln.m_QpGradU.resize(dofs+1,0.0);
    m_LocalElmtSoln.m_QpGradUold.resize(dofs+1,0.0);
    m_LocalElmtSoln.m_QpGradUolder.resize(dofs+1,0.0);

    m_LocalElmtSoln.m_QpGradV.resize(dofs+1,0.0);

    // m_Utemp.resize(468,0.0);

}

void BCSystem::releaseMemory(){

    m_BCBlockList.clear();
    m_BCBlocksNum=0;

    m_LocalK.clean();
    m_LocalR.clean();
    m_Nodes0.clear();
    m_Nodes.clear();

    m_LocalElmtSoln.m_QpU.clear();
    m_LocalElmtSoln.m_QpUold.clear();
    m_LocalElmtSoln.m_QpUolder.clear();
    m_LocalElmtSoln.m_QpV.clear();

    m_LocalElmtSoln.m_QpGradU.clear();
    m_LocalElmtSoln.m_QpGradUold.clear();
    m_LocalElmtSoln.m_QpGradUolder.clear();

    m_LocalElmtSoln.m_QpGradV.clear();
}

void BCSystem::addBCBlock2List(const BCBlock &t_BCBlock){
    if(m_BCBlockList.size()<1){
        m_BCBlockList.push_back(t_BCBlock);
        m_BCBlocksNum=1;
    }
    else{
        bool NotInList=true;
        for(const auto &it:m_BCBlockList){
            if(it.m_BCBlockName==t_BCBlock.m_BCBlockName){
                NotInList=false;break;
            }
        }
        if(NotInList){
            m_BCBlockList.push_back(t_BCBlock);
            m_BCBlocksNum+=1;
        }
        else{
            MessagePrinter::printErrorTxt("duplicate ["+t_BCBlock.m_BCBlockName+"] in your [bcs] sub block,"
                                          " please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
}

void BCSystem::printBCSystemInfo()const{
    MessagePrinter::printNormalTxt("Boundary condition system information summary");
    for(const auto &it:m_BCBlockList){
        it.printBCBlockInfo();
    }
    MessagePrinter::printStars();
}