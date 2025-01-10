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
//+++ Date   : 2022.05.09
//+++ Purpose: the bulkdof manager in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "DofHandler/BulkDofHandler.h"

void BulkDofHandler::init(){
    m_DofNameList.clear();
    m_DofIDList.clear();

    m_BulkElmtsNum=0;
    m_NodesNum=0;
    m_MaxDofsPerNode=0;
    m_TotalDofs=0;
    m_ActiveDofs=0;
    m_ElementalDofIDs_Global.clear();
    m_NodalDofIDs_Global.clear();
}

void BulkDofHandler::addDofName2List(const string &dofname){
    if(dofname.size()<1){
        MessagePrinter::printErrorTxt("invalid dof name, it is too short or an empty string, please check your input file");
        MessagePrinter::exitAsFem();
    }
    if(m_DofNameList.size()<1){
        m_DofNameList.push_back(dofname);
        m_DofIDList.clear();
        m_DofIDList.push_back(1);
        m_MaxDofsPerNode=1;
    }
    else{
        bool IsExist=false;
        for(const auto &it:m_DofNameList){
            if(it==dofname){
                IsExist=true;break;
            }
        }
        if(IsExist){
            MessagePrinter::printErrorTxt("can\'t add string to dof namelist, dofname="+dofname+" is already in your 'names' list");
            MessagePrinter::exitAsFem();
        }
        else{
            m_DofNameList.push_back(dofname);
            m_MaxDofsPerNode+=1;
            m_DofIDList.push_back(m_MaxDofsPerNode);
        }
    }
}