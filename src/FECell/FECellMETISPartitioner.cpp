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
//+++ Date    : 2025.01.22
//+++ Function: the fe cell partitioner based on METIS
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FECell/FECellMETISPartitioner.h"

#include "metis.h"

void FECellMETISPartitioner::partitionFECell(FECellData &t_CellData) {
    int size,rank;
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    if(rank == 0) {
        /**
         * using METIS to partition the bulk elements
         */
        vector<int> ElmtPtr;
        vector<int> ElmtNodeIDs;
        int ElmtsNum,NodesNum;

        ElmtsNum=t_CellData.BulkElmtsNum;
        NodesNum=t_CellData.NodesNum;

        ElmtPtr.resize(ElmtsNum+1,0);
        ElmtNodeIDs.clear();
        for (int e=1;e<=ElmtsNum;e++) {
            ElmtPtr[e]=ElmtPtr[e-1]+t_CellData.MeshCell_Total[e-1].NodesNumPerElmt;
            for (int i=1;i<=t_CellData.MeshCell_Total[e-1].NodesNumPerElmt;i++) {
                ElmtNodeIDs.push_back(t_CellData.MeshCell_Total[e-1].ElmtConn[i-1]);
            }
        }

        // int options[METIS_NOPTIONS];
        // METIS_SetDefaultOptions(options);
        // options[METIS_OPTION_NUMBERING]=1;// let the partition rank start from 0

        int CommonNodes,ObjVal;
        vector<int> NodesPartInfo;

        CommonNodes=1;
        if(t_CellData.MaxDim==3){
            CommonNodes=t_CellData.NodesNumPerSurfElmt;
        }
        else if(t_CellData.MaxDim==2){
            CommonNodes=t_CellData.NodesNumPerLineElmt;
        }

        if (size>ElmtsNum) {
            MessagePrinter::printErrorTxt("The mesh number(="+to_string(ElmtsNum)+") is less than the total cpus(="+to_string(size)+"),"
                                             "METIS partitioner will not work for this");
            MessagePrinter::exitAsFem();
        }

        t_CellData.RanksElmtsNum_Global.resize(size,0);
        t_CellData.BulkCellPartionInfo_Global.resize(ElmtsNum,0);
        NodesPartInfo.resize(NodesNum,0);

        if (size>1) {
            // int error=METIS_PartMeshDual(&ElmtsNum,
            //                              &NodesNum,
            //                              ElmtPtr.data(),
            //                              ElmtNodeIDs.data(),
            //                              NULL,
            //                              NULL,
            //                              &CommonNodes,
            //                              &size,
            //                              NULL,
            //                              options,
            //                              &ObjVal,
            //                              t_CellData.BulkCellPartionInfo_Global.data(),
            //                              NodesPartInfo.data());
            if(CommonNodes||ObjVal){}
            int error=1;
            /**
             * If you use the METIS pacakge from PETSc, the partition rank-id is start from 0,
             * however, if one call the METIS from the original METIS package, the id will start from 1 instead of 0. 
             * Be careful !!!
             */
            if(error==METIS_OK){
                MessagePrinter::printNormalTxt("METIS partition is done, everyting is ok");
            }
            else if(error==METIS_ERROR_INPUT){
                MessagePrinter::printNormalTxt("METIS partition fails, you have an input error");
                MessagePrinter::exitAsFem();
            }
            else if(error==METIS_ERROR_MEMORY){
                MessagePrinter::printNormalTxt("METIS partition fails, it could not allocate the required memory");
                MessagePrinter::exitAsFem();
            }
            else if(error==METIS_ERROR){
                MessagePrinter::printNormalTxt("METIS partition fails, you have other types of error");
                MessagePrinter::exitAsFem();
            }

            // If one call the METIS package from PETSc, the rank will automatically start from 0, instead of 1 !!!
            // then you don\'t need the following lines.
            int MaxRank,MinRank;
            MaxRank=-1;MinRank=1000000;
            for (int e=1;e<=ElmtsNum;e++) {
                if(t_CellData.BulkCellPartionInfo_Global[e-1]>MaxRank) MaxRank=t_CellData.BulkCellPartionInfo_Global[e-1];
                if(t_CellData.BulkCellPartionInfo_Global[e-1]<MinRank) MinRank=t_CellData.BulkCellPartionInfo_Global[e-1];
                // if (t_CellData.BulkCellPartionInfo_Global[e-1]<1||
                //     t_CellData.BulkCellPartionInfo_Global[e-1]>size) {
                //     // MessagePrinter::printErrorTxt("the rank id of "+to_string(e)+"-th bulk element is invalid(rank="+to_string(t_CellData.BulkCellPartionInfo_Global[e-1])+"), your METIS partitioner fails");
                //     // MessagePrinter::exitAsFem();
                // }
                // else {
                //     t_CellData.BulkCellPartionInfo_Global[e-1]-=1;// let the rank id start from 0, instead of 1 !!!
                // }
            }
            MessagePrinter::printNormalTxt("METIS partition: max rank="+to_string(MaxRank)+", min rank="+to_string(MinRank));
        }
        else {
            // if one has only 1 cpu, then no partition is required anymore
            t_CellData.BulkCellPartionInfo_Global.resize(ElmtsNum,0);
        }

        map<int,vector<SingleMeshCell>> Rank2LocalCellMap;
        int rankid;
        Rank2LocalCellMap.clear();
        for (int e=1;e<=ElmtsNum;e++) {
            rankid=t_CellData.BulkCellPartionInfo_Global[e-1];
            Rank2LocalCellMap[rankid].push_back(t_CellData.MeshCell_Total[e-1]);
            t_CellData.RanksElmtsNum_Global[rankid]+=1;
        }

        /**
         * Now, we can start to distribute the fecell
         */
        int cpuid;
        int iStart,iEnd,ranksize;
        vector<SingleMeshCell> LocalCellVector;
        vector<int> nodeids;
        int phyid,dim,phynum,datasize;
        string phyname;

        t_CellData.PhyID2MeshCellVectorMap_Local.clear();
        t_CellData.PhyName2MeshCellVectorMap_Local.clear();

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
             * For bulk fe cell
             */
            if (Rank2LocalCellMap[cpuid].size()<1) {
                MessagePrinter::printErrorTxt("the "+to_string(cpuid)+" owns 0 local FE cell, which is not allowed by METIS partitioner");
                MessagePrinter::exitAsFem();
            }
            if (cpuid==0) {
                t_CellData.MeshCell_Local=Rank2LocalCellMap[cpuid];
            }
            else {
                MPIDataBus::sendMeshCellToOthers(Rank2LocalCellMap[cpuid],1000*cpuid+100,cpuid);
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
    }// end-of-master-rank
    else {
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
    }
}