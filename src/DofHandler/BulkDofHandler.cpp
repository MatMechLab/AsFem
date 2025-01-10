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

BulkDofHandler::BulkDofHandler(){
    m_DofNameList.clear();
    m_DofIDList.clear();

    m_BulkElmtsNum=0;
    m_NodesNum=0;
    m_MaxDofsPerNode=0;
    m_MaxDofsPerElmt=0;
    m_TotalDofs=0;
    m_ActiveDofs=0;
    m_ElementalDofIDs_Global.clear();
    m_NodalDofIDs_Global.clear();
}

void BulkDofHandler::releaseMemory(){
    m_DofNameList.clear();
    m_DofIDList.clear();

    m_BulkElmtsNum=0;
    m_NodesNum=0;
    m_MaxDofsPerNode=0;
    m_MaxDofsPerElmt=0;
    m_TotalDofs=0;
    m_ActiveDofs=0;
    m_ElementalDofIDs_Global.clear();
    m_NodalDofIDs_Global.clear();
}
BulkDofHandler::~BulkDofHandler(){
    m_DofNameList.clear();
    m_DofIDList.clear();

    m_BulkElmtsNum=0;
    m_NodesNum=0;
    m_MaxDofsPerNode=0;
    m_MaxDofsPerElmt=0;
    m_TotalDofs=0;
    m_ActiveDofs=0;
    m_ElementalDofIDs_Global.clear();
    m_NodalDofIDs_Global.clear();
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

    snprintf(buff,69,"  max nnz=%9d, max row nnz=%9d",getMaxNNZ(),getMaxRowNNZ());
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
            out<<"*** total number of bulk elements="<<m_BulkElmtsNum
               <<", total nodes="<<m_NodesNum<<endl;
            out<<"*** total dofs="<<m_TotalDofs
               <<", active dofs="<<m_ActiveDofs
               <<", max dofs per node="<<m_MaxDofsPerNode
               <<", max dofs per elmt="<<m_MaxDofsPerElmt<<endl;
            for(int e=1;e<=m_BulkElmtsNum;e++){
                str="*** "+to_string(e)+"-th element: ";
                for(const auto &dofid:m_ElementalDofIDs_Global[e-1]) str+=to_string(dofid)+" ";
                out<<str<<endl;
            }
            out.close();
        }
    }
    else{
        cout<<"*** total number of bulk elements="<<m_BulkElmtsNum
            <<", total nodes="<<m_NodesNum<<endl;
        cout<<"*** total dofs="<<m_TotalDofs
            <<", active dofs="<<m_ActiveDofs
            <<", max dofs per node="<<m_MaxDofsPerNode
            <<", max dofs per elmt="<<m_MaxDofsPerElmt<<endl;
        string str;
        for(int e=1;e<=m_BulkElmtsNum;e++){
            str="*** "+to_string(e)+"-th element: ";
            for(const auto &dofid:m_ElementalDofIDs_Global[e-1]) str+=to_string(dofid)+" ";
            cout<<str<<endl;
        }
    }
}