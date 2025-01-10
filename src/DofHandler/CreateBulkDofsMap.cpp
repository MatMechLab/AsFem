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
//+++ Date   : 2022.07.23
//+++ Purpose: generate the dofs map for bulk elements
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "DofHandler/BulkDofHandler.h"
#include "MPIUtils/MPIDataBus.h"

void BulkDofHandler::createBulkDofsMap(FECell &t_fecell,const ElmtSystem &t_elmtSystem){
    
    m_BulkElmtsNum=t_fecell.getFECellBulkElmtsNum();
    m_NodesNum=t_fecell.getFECellNodesNum();
    m_TotalDofs=m_NodesNum*m_MaxDofsPerNode;// maxdofs_pernode dependes on the dofs name 
                                            // one defined in the input file
    m_MaxDofsPerElmt=m_MaxDofsPerNode*t_fecell.getFECellNodesNumPerBulkElmt();


    int rank,size;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if(rank==0){
        // allocate memory for nodal and elemental dofs map
        m_NodalDofIDs_Global.resize(m_NodesNum,vector<int>(m_MaxDofsPerNode,0));
        m_ElementalDofIDs_Global.resize(m_BulkElmtsNum,vector<int>(m_MaxDofsPerElmt,0));
        m_ActiveDofs=0;
        int dofid,elmtid,nodeid;
        vector<int> vec;
        for(const auto &block:t_elmtSystem.getBulkElmtBlockList()){
            for(const auto &name:block.m_DomainNameList){
                vec=t_fecell.getFECellBulkElmtIDSetViaPhyName(name);
                for(int e=1;e<=static_cast<int>(vec.size());e++){
                    elmtid=vec[e-1];
                    for(int i=1;i<=t_fecell.getFECellNodesNumPerBulkElmt();i++){
                        nodeid=t_fecell.getFECellIthBulkElmtJthNodeID(elmtid,i);
                        for(int k=0;k<static_cast<int>(block.m_DofIDs.size());k++){
                            dofid=block.m_DofIDs[k];
                            if(m_NodalDofIDs_Global[nodeid-1][dofid-1]==0){
                                // if current nodal dof flag is zero, then we assign a non-zero value to it
                                // otherwise, it is duplicated, then skip it.
                                m_ActiveDofs+=1;
                                // m_NodalDofIDs[nodeid-1][dofid-1]=m_ActiveDofs;
                                m_NodalDofIDs_Global[nodeid-1][dofid-1]=(nodeid-1)*m_MaxDofsPerNode+dofid;
                            }
                        }
                    }
                }
            }
        }
        
        // check the dofs status for each node
        bool HasDofID=false;
        for(int i=0;i<m_NodesNum;i++){
            HasDofID=false;
            for(int j=0;j<m_MaxDofsPerNode;j++){
                if(m_NodalDofIDs_Global[i][j]){
                    HasDofID=true;break;
                }
            }
            if(!HasDofID){
                MessagePrinter::printErrorTxt("Node-"+to_string(i+1)+" hasen\'t been assigned by the dof, please check your code");
                MessagePrinter::exitAsFem();
            }
        }
        
        // now we can create the elemental dofs map
        vector<vector<int>> RowDofIDs;
        RowDofIDs.resize(m_ActiveDofs);
        for(int e=1;e<=m_BulkElmtsNum;e++){
            for(int j=1;j<=t_fecell.getFECellIthBulkElmtNodesNum(e);j++){
                nodeid=t_fecell.getFECellIthBulkElmtJthNodeID(e,j);
                for(int k=1;k<=m_MaxDofsPerNode;k++){
                    dofid=m_NodalDofIDs_Global[nodeid-1][k-1];
                    m_ElementalDofIDs_Global[e-1][(j-1)*m_MaxDofsPerNode+k-1]=dofid;
                }
            }
            // remove the zero elements, only keep the nonzero part
            m_ElementalDofIDs_Global[e-1].erase(
                remove(m_ElementalDofIDs_Global[e-1].begin(),m_ElementalDofIDs_Global[e-1].end(),0),
                m_ElementalDofIDs_Global[e-1].end()
                );
            m_ElementalDofIDs_Global[e-1].shrink_to_fit();
            if(static_cast<int>(m_ElementalDofIDs_Global[e-1].size())!=m_MaxDofsPerElmt){
                MessagePrinter::printErrorTxt("the length of element dof ids is not equal to MaxDofsPerElmt in CreateBulkDofsMap");
                MessagePrinter::exitAsFem();
            }
            for(int i=0;i<static_cast<int>(m_ElementalDofIDs_Global[e-1].size());i++){
                for(int j=0;j<static_cast<int>(m_ElementalDofIDs_Global[e-1].size());j++){
                    RowDofIDs[m_ElementalDofIDs_Global[e-1][i]-1].push_back(m_ElementalDofIDs_Global[e-1][j]);
                }
            }
        }
        // now we remove duplicate dofs in each row
        m_MaxNNZ=0;m_MaxRowNNZ=-1;
        for(int i=0;i<m_ActiveDofs;i++){
            RowDofIDs[i].erase(
                remove(RowDofIDs[i].begin(),RowDofIDs[i].end(),0),
                RowDofIDs[i].end()
                );
            RowDofIDs[i].shrink_to_fit();
            if(static_cast<int>(RowDofIDs[i].size())>m_MaxRowNNZ) m_MaxRowNNZ=static_cast<int>(RowDofIDs[i].size());
            m_MaxNNZ+=static_cast<int>(RowDofIDs[i].size());
        }
        RowDofIDs=vector<vector<int>>(0);// release memory


        /**
         * send out the dof map info to other ranks
         */
        m_ElementalDofIDs_Local.clear();
        map<int,vector<vector<int>>> Rank2ElementalDofIDsMap;
        int cpuid,tag;
        Rank2ElementalDofIDsMap.clear();
        for(int e=1;e<=t_fecell.getFECellBulkElmtsNum();e++){
            cpuid=t_fecell.getFECellIthBulkElmtRankID(e);
            if(cpuid==0){
                // for master rank
                m_ElementalDofIDs_Local.push_back(m_ElementalDofIDs_Global[e-1]);
            }
            else{
                Rank2ElementalDofIDsMap[cpuid].push_back(m_ElementalDofIDs_Global[e-1]);//stores the local element dof id of each local rank
            }
        }
        for(cpuid=1;cpuid<size;cpuid++){
            if(static_cast<int>(Rank2ElementalDofIDsMap[cpuid].size())!=t_fecell.getIthRanksBulkElmtsNum(cpuid)){
                MessagePrinter::printErrorTxt("Rank-"+to_string(cpuid)+" owns different local bulk elmts number from the one given in FECell,"+
                                              "please check your code in CreateBulkDofsMap.cpp");
                MessagePrinter::exitAsFem();
            }
            tag=cpuid*1000+0;
            MPIDataBus::sendIntegerToOthers(m_MaxNNZ,tag,cpuid);
            tag=cpuid*1000+1;
            MPIDataBus::sendIntegerToOthers(m_MaxRowNNZ,tag,cpuid);
            // send out nodal dof id and elemental dof id
            tag=cpuid*1000+2;
            MPIDataBus::sendAlignedVectorOfIntegerVectorToOthers(m_NodalDofIDs_Global,tag,cpuid);
            tag=cpuid*1000+3;
            MPIDataBus::sendAlignedVectorOfIntegerVectorToOthers(m_ElementalDofIDs_Global,tag,cpuid);
            tag=cpuid*1000+4;
            MPIDataBus::sendAlignedVectorOfIntegerVectorToOthers(Rank2ElementalDofIDsMap[cpuid],tag,cpuid);
        }
        Rank2ElementalDofIDsMap.clear();
    }// end-of-master-rank
    else{
        int tag;
        tag=rank*1000+0;
        MPIDataBus::receiveIntegerFromMaster(m_MaxNNZ,tag);
        tag=rank*1000+1;
        MPIDataBus::receiveIntegerFromMaster(m_MaxRowNNZ,tag);
        tag=rank*1000+2;
        MPIDataBus::receiveAlignedVectorOfIntegerVectorFromMaster(m_NodalDofIDs_Global,tag);
        tag=rank*1000+3;
        MPIDataBus::receiveAlignedVectorOfIntegerVectorFromMaster(m_ElementalDofIDs_Global,tag);
        tag=rank*1000+4;
        MPIDataBus::receiveAlignedVectorOfIntegerVectorFromMaster(m_ElementalDofIDs_Local,tag);
    }
    MPI_Bcast(&m_ActiveDofs,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // reset the dof id of each local element in the FECell class
    for (int e=1;e<=t_fecell.getLocalFECellBulkElmtsNum();e++) {
        if (t_fecell.getLocalBulkFECellVecRef()[e-1].ElmtDofIDs.size()==m_ElementalDofIDs_Local[e-1].size()){
            for (int i=0;i<static_cast<int>(t_fecell.getLocalBulkFECellVecRef()[e-1].ElmtDofIDs.size());i++) {
                t_fecell.getLocalBulkFECellVecRef()[e-1].ElmtDofIDs[i]=m_ElementalDofIDs_Local[e-1][i];
            }
        }
        else {
            t_fecell.getLocalBulkFECellVecRef()[e-1].ElmtDofIDs=m_ElementalDofIDs_Local[e-1];
        }
    }
}