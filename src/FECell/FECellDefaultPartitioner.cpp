//****************************************************************
//* This file is part of the AsFem framework
//* Advanced Simulation kit based on Finite Element Method (AsFem)
//* All rights reserved, Yang Bai/MM-Lab@CopyRight 2020-present
//* https://github.com/MatMechLab/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author  : Yang Bai
//+++ Date    : 2024.08.01
//+++ Function: the built-in fe cell partitioner
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FECell/FECellDefaultPartitioner.h"


void FECellDefaultPartitioner::partitionFECell(FECellData &t_CellData){

    int size,rank;
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    if(rank==0){
        /**
         * Now we start to distribute the gloabl mesh cell into different ranks
        */
        int cpuid;
        int iStart,iEnd,ranksize;
        vector<SingleMeshCell> LocalCellVector;
        vector<int> nodeids;
        int phyid,dim,phynum,datasize;
        string phyname;

        t_CellData.PhyID2MeshCellVectorMap_Local.clear();
        t_CellData.PhyName2MeshCellVectorMap_Local.clear();

        t_CellData.BulkCellPartionInfo_Global.resize(t_CellData.BulkElmtsNum,0);
        t_CellData.RanksElmtsNum_Global.resize(size,0);
        t_CellData.NodeIDs_Local.clear();

        for(cpuid=1;cpuid<size;cpuid++){
            /**
             * send out the physical group info to other ranks
             */
            MPIDataBus::sendIntegerToOthers(t_CellData.PhyGroupNum_Global,1000*cpuid+0,cpuid);
            for(phynum=0;phynum<t_CellData.PhyGroupNum_Global;phynum++){
                phyid=t_CellData.PhyIDVector_Global[phynum];
                phyname=t_CellData.PhyNameVector_Global[phynum];
                datasize=t_CellData.PhyGroupElmtsNumVector_Global[phynum];
                dim=t_CellData.PhyDimVector_Global[phynum];
                MPIDataBus::sendIntegerToOthers(phyid,1000*cpuid+10,cpuid);
                MPIDataBus::sendIntegerToOthers(dim,1000*cpuid+20,cpuid);
                MPIDataBus::sendIntegerToOthers(datasize,1000*cpuid+30,cpuid);
                MPIDataBus::sendStringToOthers(phyname,1000*cpuid+40,cpuid);
            }
        }

        for(cpuid=0;cpuid<size;cpuid++){
            /**
             * for elemental physical cell set
             */
            for(phynum=0;phynum<t_CellData.PhyGroupNum_Global;phynum++){
                phyid=t_CellData.PhyIDVector_Global[phynum];
                phyname=t_CellData.PhyNameVector_Global[phynum];

                datasize=static_cast<int>(t_CellData.PhyID2MeshCellVectorMap_Global[phyid].size());
                ranksize=datasize/size;
                iStart=cpuid*ranksize;
                iEnd=(cpuid+1)*ranksize;
                if(cpuid==size-1) iEnd=datasize;

                if (datasize<size) {
                    // for the case where the mesh number is less than the cpu numbers
                    ranksize=1;
                    iStart=cpuid*ranksize;
                    iEnd=(cpuid+1)*ranksize;
                    if (iStart>datasize-2) {
                        iStart=0;
                        iEnd=0;
                    }
                    if(cpuid==size-1) {
                        iStart=datasize-1;
                        iEnd=datasize;
                    }
                }
                if (datasize==0) {
                    iStart=0;iEnd=0;
                }

                LocalCellVector.clear();
                for(int e=iStart;e<iEnd;e++){
                    LocalCellVector.push_back(t_CellData.PhyID2MeshCellVectorMap_Global[phyid][e]);
                }

                if(cpuid==0){
                    t_CellData.PhyID2MeshCellVectorMap_Local[phyid]=LocalCellVector;
                    t_CellData.PhyName2MeshCellVectorMap_Local[phyname]=LocalCellVector;
                }
                else{
                    MPIDataBus::sendPhyID2MeshCellMapToOthers(phyid,LocalCellVector,1000*cpuid+60,cpuid);
                    MPIDataBus::sendPhyName2MeshCellMapToOthers(phyname,LocalCellVector,1000*cpuid+80,cpuid);
                }
            }// end-of-physical-group-loop
            /**
             * for all domain bulk cell
             */
            datasize=t_CellData.BulkElmtsNum;
            ranksize=datasize/size;
            iStart=cpuid*ranksize;
            iEnd=(cpuid+1)*ranksize;
            if(cpuid==size-1) iEnd=datasize;
            if (datasize<size) {
                // for the case where the mesh number is less than the cpu numbers
                ranksize=1;
                iStart=cpuid*ranksize;
                iEnd=(cpuid+1)*ranksize;
                if (iStart>datasize-2) {
                    iStart=0;
                    iEnd=0;
                }
                if(cpuid==size-1) {
                    iStart=datasize-1;
                    iEnd=datasize;
                }
            }
            if (datasize==0) {
                iStart=0;iEnd=0;
            }

            LocalCellVector.clear();
            for(int e=iStart;e<iEnd;e++){
                LocalCellVector.push_back(t_CellData.MeshCell_Total[e]);
                t_CellData.BulkCellPartionInfo_Global[e]=cpuid;
                t_CellData.RanksElmtsNum_Global[cpuid]+=1;
            }

            if(cpuid==0){
                t_CellData.MeshCell_Local=LocalCellVector;
            }
            else{
                MPIDataBus::sendMeshCellToOthers(LocalCellVector,1000*cpuid+100,cpuid);
            }

            /**
             * for nodal physical group sets
             */
            if(cpuid!=0) MPIDataBus::sendIntegerToOthers(t_CellData.NodalPhyGroupNum_Global,1000*cpuid+200,cpuid);
            for(phynum=0;phynum<t_CellData.NodalPhyGroupNum_Global;phynum++){
                phyid=t_CellData.NodalPhyIDVector_Global[phynum];
                datasize=t_CellData.NodalPhyGroupNodesNumVector_Global[phynum];
                phyname=t_CellData.NodalPhyNameVector_Global[phynum];

                if(cpuid!=0) MPIDataBus::sendIntegerToOthers(phyid,1000*cpuid+210,cpuid);
                if(cpuid!=0) MPIDataBus::sendIntegerToOthers(datasize,1000*cpuid+220,cpuid);
                if(cpuid!=0) MPIDataBus::sendStringToOthers(phyname,1000*cpuid+230,cpuid);

                ranksize=static_cast<int>(t_CellData.NodalPhyName2NodeIDVecMap_Global[phyname].size())/size;
                iStart=cpuid*ranksize;
                iEnd=(cpuid+1)*ranksize;
                if(cpuid==size-1) iEnd=static_cast<int>(t_CellData.NodalPhyName2NodeIDVecMap_Global[phyname].size());

                nodeids.clear();
                for(int e=iStart;e<iEnd;e++){
                    nodeids.push_back(t_CellData.NodalPhyName2NodeIDVecMap_Global[phyname][e]);
                }
                if(cpuid==0){
                    t_CellData.NodalPhyName2NodeIDVecMap_Local[phyname]=nodeids;
                }
                else{
                    MPIDataBus::sendPhyName2NodeIDVecMapToOthers(phyname,nodeids,1000*cpuid+240,cpuid);
                }
            }
        }// end-of-cpuid-loop

        /**
         * Now we send the BulkCellPartionInfo_Global and RanksElmtsNum_Global to other ranks,
         * these two vectors are shared among ranks
         */
        for(cpuid=1;cpuid<size;cpuid++){
            MPIDataBus::sentIntegerVectorToOthers(t_CellData.BulkCellPartionInfo_Global,1000*cpuid+260,cpuid);
            MPIDataBus::sentIntegerVectorToOthers(t_CellData.RanksElmtsNum_Global,1000*cpuid+280,cpuid);
        }
        /**
         * distribute and send the node ids
         */
        vector<int> NodeIDs;
        for(cpuid=0;cpuid<size;cpuid++) {
            datasize=t_CellData.NodesNum;
            ranksize=datasize/size;
            iStart=cpuid*ranksize;
            iEnd=(cpuid+1)*ranksize;
            if(cpuid==size-1) iEnd=datasize;
            NodeIDs.clear();
            for (int i=iStart;i<iEnd;i++) {
                NodeIDs.push_back(i+1);
            }
            if (cpuid==0) {
                t_CellData.NodeIDs_Local=NodeIDs;
            }
            else {
                MPIDataBus::sentIntegerVectorToOthers(NodeIDs,1000*cpuid+300,cpuid);
            }
        }
    }// end of rank-0
    else{
        /**
         * receive mesh cell from master rank
         */
        vector<int> nodeids;
        int phyid,dim,phynum,datasize;
        string phyname;

        t_CellData.PhyID2MeshCellVectorMap_Local.clear();
        t_CellData.PhyName2MeshCellVectorMap_Local.clear();

        MPIDataBus::receiveIntegerFromMaster(t_CellData.PhyGroupNum_Global,1000*rank+0);

        t_CellData.PhyIDVector_Global.clear();
        t_CellData.PhyDimVector_Global.clear();
        t_CellData.PhyGroupElmtsNumVector_Global.clear();
        t_CellData.PhyNameVector_Global.clear();
        t_CellData.PhyID2NameMap_Global.clear();
        t_CellData.PhyName2IDMap_Global.clear();
        for(phynum=0;phynum<t_CellData.PhyGroupNum_Global;phynum++){
            MPIDataBus::receiveIntegerFromMaster(phyid,1000*rank+10);
            MPIDataBus::receiveIntegerFromMaster(dim,1000*rank+20);
            MPIDataBus::receiveIntegerFromMaster(datasize,1000*rank+30);
            MPIDataBus::receiveStringFromMaster(phyname,1000*rank+40);

            t_CellData.PhyIDVector_Global.push_back(phyid);
            t_CellData.PhyDimVector_Global.push_back(dim);
            t_CellData.PhyGroupElmtsNumVector_Global.push_back(datasize);
            t_CellData.PhyNameVector_Global.push_back(phyname);

            t_CellData.PhyID2NameMap_Global[phyid]=phyname;
            t_CellData.PhyName2IDMap_Global[phyname]=phyid;

            MPIDataBus::receivePhyID2MeshCellMapFromMaster(t_CellData.PhyID2MeshCellVectorMap_Local,1000*rank+60);
            MPIDataBus::receivePhyName2MeshCellMapFromMaster(t_CellData.PhyName2MeshCellVectorMap_Local,1000*rank+80);
        }
        /**
         * for all domain
         */
        MPIDataBus::receiveMeshCellFromMaster(t_CellData.MeshCell_Local,1000*rank+100);

        t_CellData.NodeIDs_Local.clear();
        for (const auto &cell:t_CellData.MeshCell_Local) {
            for (const auto &id:cell.ElmtConn) {
                t_CellData.NodeIDs_Local.push_back(id);
            }
        }
        // remove the duplicate ids
        sort(t_CellData.NodeIDs_Local.begin(),t_CellData.NodeIDs_Local.end());
        t_CellData.NodeIDs_Local.erase(unique(t_CellData.NodeIDs_Local.begin(),t_CellData.NodeIDs_Local.end()),t_CellData.NodeIDs_Local.end());

        /**
         * for nodal physical group set
         */
        MPIDataBus::receiveIntegerFromMaster(t_CellData.NodalPhyGroupNum_Global,1000*rank+200);
        t_CellData.NodalPhyGroupNodesNumVector_Global.clear();
        t_CellData.NodalPhyIDVector_Global.clear();
        t_CellData.NodalPhyNameVector_Global.clear();
        t_CellData.NodalPhyID2NameMap_Global.clear();
        t_CellData.NodalPhyName2IDMap_Global.clear();
        t_CellData.NodalPhyName2NodeIDVecMap_Local.clear();
        for(phynum=0;phynum<t_CellData.NodalPhyGroupNum_Global;phynum++){
            MPIDataBus::receiveIntegerFromMaster(phyid,1000*rank+210);
            MPIDataBus::receiveIntegerFromMaster(datasize,1000*rank+220);
            MPIDataBus::receiveStringFromMaster(phyname,1000*rank+230);
            MPIDataBus::receivePhyName2NodeIDVecMapFromMaster(t_CellData.NodalPhyName2NodeIDVecMap_Local,1000*rank+240);

            t_CellData.NodalPhyIDVector_Global.push_back(phyid);
            t_CellData.NodalPhyGroupNodesNumVector_Global.push_back(datasize);
            t_CellData.NodalPhyNameVector_Global.push_back(phyname);

            t_CellData.NodalPhyID2NameMap_Global[phyid]=phyname;
            t_CellData.NodalPhyName2IDMap_Global[phyname]=phyid;
        }

        /**
         * Receive bulk element's partition info and ranks elemts number info from master rank
         */
        MPIDataBus::receiveIntegerVectorFromMaster(t_CellData.BulkCellPartionInfo_Global,1000*rank+260);
        MPIDataBus::receiveIntegerVectorFromMaster(t_CellData.RanksElmtsNum_Global,1000*rank+280);
        /**
         * for the local node ids
         */
        MPIDataBus::receiveIntegerVectorFromMaster(t_CellData.NodeIDs_Local,1000*rank+300);
    }// end-of-other-ranks
}