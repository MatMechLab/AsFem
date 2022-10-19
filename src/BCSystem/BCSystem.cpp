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
//+++ Purpose: Define the boundary condition system in AsFem
//+++          Here we can apply:
//+++               1) dirichlet bc, i.e. displacement, temperature ...
//+++               2) neuman bc, i.e. flux, force
//+++               3) robin bc as well as user-defined-bc (ubc)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/BCSystem.h"

BCSystem::BCSystem(){
    m_bclock_list.clear();
    m_bcblocks_num=0;
    m_dirichlet_penalty=1.0e16;
}

void BCSystem::init(const int &dofs){
    m_localK.resize(dofs+1,dofs+1,0.0);
    m_localR.resize(dofs+1,0.0);
    m_nodes0.resize(27+1);
    m_nodes.resize(27+1);

    m_local_elmtsoln.m_gpU.resize(dofs+1,0.0);
    m_local_elmtsoln.m_gpUold.resize(dofs+1,0.0);
    m_local_elmtsoln.m_gpUolder.resize(dofs+1,0.0);
    m_local_elmtsoln.m_gpV.resize(dofs+1,0.0);

    m_local_elmtsoln.m_gpGradU.resize(dofs+1,0.0);
    m_local_elmtsoln.m_gpGradUold.resize(dofs+1,0.0);
    m_local_elmtsoln.m_gpGradUolder.resize(dofs+1,0.0);

    m_local_elmtsoln.m_gpGradV.resize(dofs+1,0.0);

    // m_Utemp.resize(468,0.0);

}

void BCSystem::releaseMemory(){

    m_bclock_list.clear();
    m_bcblocks_num=0;

    m_localK.clean();
    m_localR.clean();
    m_nodes0.clear();
    m_nodes.clear();

    m_local_elmtsoln.m_gpU.clear();
    m_local_elmtsoln.m_gpUold.clear();
    m_local_elmtsoln.m_gpUolder.clear();
    m_local_elmtsoln.m_gpV.clear();

    m_local_elmtsoln.m_gpGradU.clear();
    m_local_elmtsoln.m_gpGradUold.clear();
    m_local_elmtsoln.m_gpGradUolder.clear();

    m_local_elmtsoln.m_gpGradV.clear();
}

void BCSystem::addBCBlock2List(const BCBlock &t_bcblock){
    if(m_bclock_list.size()<1){
        m_bclock_list.push_back(t_bcblock);
        m_bcblocks_num=1;
    }
    else{
        bool NotInList=true;
        for(const auto &it:m_bclock_list){
            if(it.m_bcBlockName==t_bcblock.m_bcBlockName){
                NotInList=false;break;
            }
        }
        if(NotInList){
            m_bclock_list.push_back(t_bcblock);
            m_bcblocks_num+=1;
        }
        else{
            MessagePrinter::printErrorTxt("duplicate ["+t_bcblock.m_bcBlockName+"] in your [bcs] sub block,"
                                          " please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
}

void BCSystem::printBCSystemInfo()const{
    MessagePrinter::printNormalTxt("Boundary condition system information summary");
    for(const auto &it:m_bclock_list){
        it.printBCBlockInfo();
    }
    MessagePrinter::printStars();
}