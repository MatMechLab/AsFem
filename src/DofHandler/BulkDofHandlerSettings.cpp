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
//+++ Date   : 2022.05.09
//+++ Purpose: the bulkdof manager in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "DofHandler/BulkDofHandler.h"

void BulkDofHandler::init(){
    m_dof_namelist.clear();
    m_dof_idlist.clear();

    m_bulkelmts=0;
    m_nodes=0;
    m_maxdofs_pernode=0;
    m_total_dofs=0;
    m_active_dofs=0;
    m_elmt_dofids.clear();
    m_nodal_dofids.clear();
}

void BulkDofHandler::addDofName2List(const string &dofname){
    if(dofname.size()<1){
        MessagePrinter::printErrorTxt("invalid dof name, it is too short or an empty string, please check your input file");
        MessagePrinter::exitAsFem();
    }
    if(m_dof_namelist.size()<1){
        m_dof_namelist.push_back(dofname);
        m_dof_idlist.clear();
        m_dof_idlist.push_back(1);
        m_maxdofs_pernode=1;
    }
    else{
        bool IsExist=false;
        for(const auto &it:m_dof_namelist){
            if(it==dofname){
                IsExist=true;break;
            }
        }
        if(IsExist){
            MessagePrinter::printErrorTxt("can\'t add string to dof namelist, dofname="+dofname+" is already in your 'names' list");
            MessagePrinter::exitAsFem();
        }
        else{
            m_dof_namelist.push_back(dofname);
            m_maxdofs_pernode+=1;
            m_dof_idlist.push_back(m_maxdofs_pernode);
        }
    }
}