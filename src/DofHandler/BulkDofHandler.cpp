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

BulkDofHandler::BulkDofHandler(){
    m_dof_namelist.clear();
    m_dof_idlist.clear();

    m_bulkelmts=0;
    m_nodes=0;
    m_maxdofs_pernode=0;
    m_maxdofs_perelmt=0;
    m_total_dofs=0;
    m_active_dofs=0;
    m_elmt_dofids.clear();
    m_nodal_dofids.clear();
}

void BulkDofHandler::releaseMemory(){
    m_dof_namelist.clear();
    m_dof_idlist.clear();

    m_bulkelmts=0;
    m_nodes=0;
    m_maxdofs_pernode=0;
    m_maxdofs_perelmt=0;
    m_total_dofs=0;
    m_active_dofs=0;
    m_elmt_dofids.clear();
    m_nodal_dofids.clear();
}
BulkDofHandler::~BulkDofHandler(){
    m_dof_namelist.clear();
    m_dof_idlist.clear();

    m_bulkelmts=0;
    m_nodes=0;
    m_maxdofs_pernode=0;
    m_maxdofs_perelmt=0;
    m_total_dofs=0;
    m_active_dofs=0;
    m_elmt_dofids.clear();
    m_nodal_dofids.clear();
}

void BulkDofHandler::printBulkDofsInfo()const{
    char buff[69];
    string str;
    MessagePrinter::printStars();
    MessagePrinter::printNormalTxt("bulk dofs information summary");
    snprintf(buff,69,"  total dofs=%9d, active dofs=%9d, max dofs per node=%2d",getTotalDofs(),
                                                                                getActiveDofs(),
                                                                                getMaxDofsPerNode());
    str=string(buff);
    MessagePrinter::printNormalTxt(str);
    MessagePrinter::printNormalTxt("  dofs id                                     dofs name");
    for(int i=1;i<=getMaxDofsPerNode();i++){
        snprintf(buff,69,"   %3d                              %18s",getIthDofID(i),getIthDofName(i).c_str());
        str=string(buff);
        MessagePrinter::printNormalTxt(str);
    }
    MessagePrinter::printStars();
}
//*****************************************
void BulkDofHandler::printBulkElementalDofsInfo(const bool &flag)const{
    if(flag){
        // if flag=true, we print out the dofs map into a txt file
        PetscMPIInt rank;
        MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
        if(rank==0){
            string str;
            std::ofstream out;
            out.open("dofs.map",std::ios::out);
            if (!out.is_open()){
                MessagePrinter::printErrorTxt("can\'t create dofs.map file, please ensure you have the write permission");
                MessagePrinter::exitAsFem();
            }
            out<<"*** total number of bulk elements="<<m_bulkelmts
               <<", total nodes="<<m_nodes<<endl;
            out<<"*** total dofs="<<m_total_dofs
               <<", active dofs="<<m_active_dofs
               <<", max dofs per node="<<m_maxdofs_pernode
               <<", max dofs per elmt="<<m_maxdofs_perelmt<<endl;
            for(int e=1;e<=m_bulkelmts;e++){
                str="*** "+to_string(e)+"-th element: ";
                for(const auto &dofid:m_elmt_dofids[e-1]) str+=to_string(dofid)+" ";
                out<<str<<endl;
            }
            out.close();
        }
    }
    else{
        cout<<"*** total number of bulk elements="<<m_bulkelmts
            <<", total nodes="<<m_nodes<<endl;
        cout<<"*** total dofs="<<m_total_dofs
            <<", active dofs="<<m_active_dofs
            <<", max dofs per node="<<m_maxdofs_pernode
            <<", max dofs per elmt="<<m_maxdofs_perelmt<<endl;
        string str;
        for(int e=1;e<=m_bulkelmts;e++){
            str="*** "+to_string(e)+"-th element: ";
            for(const auto &dofid:m_elmt_dofids[e-1]) str+=to_string(dofid)+" ";
            cout<<str<<endl;
        }
    }
}