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
//+++ Date   : 2024.01.16
//+++ Purpose: Implement the message sender/receiver among mpi ranks
//+++          for mesh cell, and other data structure
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include "FECell/SingleMeshCell.h"

#include "mpi.h"

/**
 * This class implement the message send/receive API for mesh cell and other data structure among mpi ranks
*/
class MPIDataBus{
public:
    /**
     * Constructor
    */
    MPIDataBus() {}

    /**
     * Send mesh cell data from master rank to others, this should be called on master rank
     * @param meshcellvec the mesh cell vector to be sent
     * @param tag the id of the sending message
     * @param cpuid the rank id used to receive the message
    */
    static void sendMeshCellToOthers(const vector<SingleMeshCell> &meshcellvec,const int &tag,const int &cpuid);
    /**
     * Receive the mesh cell data from master rank, this should be called on each rank locally
     * @param the mesh cell vector used to store the data from master rank
     * @param tag the received message tag
     * @param cpuid the cpuid to receive the message
    */
    static void receiveMeshCellFromMaster(vector<SingleMeshCell> &meshcellvec,const int &tag);

    /**
     * Send the physical name and mesh cell to each rank
     * @param phyname the physical group name to be sent
     * @param meshcellvec the mesh cell vector to be sent
     * @param tag the id of the sending message
     * @param cpuid the rank id used to receive the message
    */
    static void sendPhyName2MeshCellMapToOthers(const string &phyname,const vector<SingleMeshCell> &meshcellvec,const int &tag,const int &cpuid);
    /**
     * Receive the physical group name and mesh cell for local map
     * @param localmap the local phyname-->meshcellvector map
     * @param tag the id of the sending message
    */
    static void receivePhyName2MeshCellMapFromMaster(map<string,vector<SingleMeshCell>> &localmap,const int &tag);

    /**
     * Send the physical id and mesh cell to each rank
     * @param phyid the physical group id to be sent
     * @param meshcellvec the mesh cell vector to be sent
     * @param tag the id of the sending message
     * @param cpuid the rank id used to receive the message
    */
    static void sendPhyID2MeshCellMapToOthers(const int &phyid,const vector<SingleMeshCell> &meshcellvec,const int &tag,const int &cpuid);
    /**
     * Receive the physical group name and mesh cell for local map
     * @param localmap the local id-->meshcellvector map
     * @param tag the id of the sending message
    */
    static void receivePhyID2MeshCellMapFromMaster(map<int,vector<SingleMeshCell>> &localmap,const int &tag);

    //*************************************************************
    //*** for nodal phy group transfer
    //*************************************************************
    /**
     * Send the nodal physical name and nodes id to each rank
     * @param phyname the physical group name to be sent
     * @param nodesid the mesh cell vector to be sent
     * @param tag the id of the sending message
     * @param cpuid the rank id used to receive the message
    */
    static void sendPhyName2NodeIDVecMapToOthers(const string &phyname,const vector<int> &nodesid,const int &tag,const int &cpuid);
    /**
     * Receive the physical group name and mesh cell for local map
     * @param localmap the local phyname-->nodeid map
     * @param tag the id of the sending message
    */
    static void receivePhyName2NodeIDVecMapFromMaster(map<string,vector<int>> &localmap,const int &tag);

};