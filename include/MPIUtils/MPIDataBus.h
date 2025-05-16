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
#include "FECell/MeshType.h"

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
     * Check if all the rank is release or some is occupied
     */
    static void printRankStatus();

    /**
     * Send integer from master rank to others, this should be called on master rank
     * @param Val the integer to be sent
     * @param Tag the id of the sending message
     * @param Cpuid the rank id used to receive the message
    */
    static void sendIntegerToOthers(const int &Val,const int &Tag,const int &Cpuid);
    /**
     * Receive the integer from master rank, this should be called on each rank locally
     * @param Val the integer value send from master rank
     * @param Tag the received message tag
    */
    static void receiveIntegerFromMaster(int &Val,const int &Tag);

    /**
     * Send string from master rank to others, this should be called on master rank
     * @param Txt the string to be sent
     * @param Tag the id of the sending message
     * @param Cpuid the rank id used to send the message
    */
    static void sendStringToOthers(const string &Txt,const int &Tag,const int &Cpuid);
    /**
     * Receive the string from master rank, this should be called on each rank locally
     * @param Txt the string used to store the data from master rank
     * @param Tag the received message tag
    */
    static void receiveStringFromMaster(string &Txt,const int &Tag);

    /**
     * Send the mesh type to other ranks
     * @param t_MeshType the meshtype to be sent
     * @param Tag the id of the sending message
     * @param Cpuid the rank id used to send the message
     */
     static void sendMeshTypeToOthers(const MeshType &t_MeshType,const int &Tag,const int &Cpuid);
     /**
      * Receive the string from master rank, this should be called on each rank locally
      * @param Txt the string used to store the data from master rank
      * @param Tag the received message tag*
      */
      static void receiveMeshTypeFromMater(MeshType &t_MeshType,const int &Tag);

    /**
     * Send the integer vector to other ranks
     * @param Vec the vector stores integers
     * @param Tag the id of the sending message
     * @param Cpuid the rank id used to receive the message
     */
    static void sentIntegerVectorToOthers(const vector<int> &Vec,const int &Tag,const int &Cpuid);
    /**
     * Receive the vector of integer vector from master
     * @param Vec the integer vector used to store the data from master rank
     * @param Tag the received message tag
     */
    static void receiveIntegerVectorFromMaster(vector<int> &Vec,const int &Tag);

    /**
     * Send the vector of integer vector to other ranks
     * @param Vec the vector stores integer vectors
     * @param Tag the id of the sending message
     * @param Cpuid the rank id used to receive the message
     */
    static void sentVectorOfIntegerVectorToOthers(const vector<vector<int>> &Vec,const int &Tag,const int &Cpuid);
    /**
     * Receive the vector of integer vector from master ranks
     * @param Vec the string used to store the data from master rank
     * @param Tag the received message tag
     */
    static void receiveVectorOfIntegerVectorFromMaster(vector<vector<int>> &Vec,const int &Tag);

    /**
     * Send the aligned vector (with the same width) of integer vector to other ranks
     * @param Vec the vector stores integer vectors
     * @param Tag the id of the sending message
     * @param Cpuid the rank id used to receive the message
     */
    static void sendAlignedVectorOfIntegerVectorToOthers(const vector<vector<int>> &Vec,const int &Tag,const int &Cpuid);
    /**
     * Receive the aligned vector (with the same width) of integer vector from master ranks
     * @param Vec the string used to store the data from master rank
     * @param Tag the received message tag
     */
    static void receiveAlignedVectorOfIntegerVectorFromMaster(vector<vector<int>> &Vec,const int &Tag);

    /**
     * Send the string vector to other ranks
     * @param Vec the vector stores strings
     * @param Tag the id of the sending message
     * @param Cpuid the rank id used to receive the message
     */
    static void sentStringVectorToOthers(const vector<string> &Vec,const int &Tag,const int &Cpuid);
    /**
     * Send the integer vector to other ranks
     * @param Vec the string used to store the data from master rank
     * @param Tag the received message tag
     */
    static void receiveStringVectorFromMaster(vector<string> &Vec,const int &Tag);

    /**
     * Send mesh cell data from master rank to others, this should be called on master rank
     * @param MeshCellVec the mesh cell vector to be sent
     * @param Tag the id of the sending message
     * @param Cpuid the rank id used to receive the message
    */
    static void sendMeshCellToOthers(const vector<SingleMeshCell> &MeshCellVec,const int &Tag,const int &Cpuid);
    /**
     * Receive the mesh cell data from master rank, this should be called on each rank locally
     * @param Meshcellvec the mesh cell vector used to store the data from master rank
     * @param Tag the received message tag
    */
    static void receiveMeshCellFromMaster(vector<SingleMeshCell> &Meshcellvec,const int &Tag);

    /**
     * Send the physical name and mesh cell to each rank
     * @param Phyname the physical group name to be sent
     * @param Meshcellvec the mesh cell vector to be sent
     * @param Tag the id of the sending message
     * @param Cpuid the rank id used to receive the message
    */
    static void sendPhyName2MeshCellMapToOthers(const string &Phyname,const vector<SingleMeshCell> &Meshcellvec,const int &Tag,const int &Cpuid);
    /**
     * Receive the physical group name and mesh cell for local map
     * @param Localmap the local phyname-->meshcellvector map
     * @param Tag the id of the sending message
    */
    static void receivePhyName2MeshCellMapFromMaster(map<string,vector<SingleMeshCell>> &Localmap,const int &Tag);

    /**
     * Send the physical id and mesh cell to each rank
     * @param Phyid the physical group id to be sent
     * @param Meshcellvec the mesh cell vector to be sent
     * @param Tag the id of the sending message
     * @param Cpuid the rank id used to receive the message
    */
    static void sendPhyID2MeshCellMapToOthers(const int &Phyid,const vector<SingleMeshCell> &Meshcellvec,const int &Tag,const int &Cpuid);
    /**
     * Receive the physical group name and mesh cell for local map
     * @param Localmap the local id-->meshcellvector map
     * @param Tag the id of the sending message
    */
    static void receivePhyID2MeshCellMapFromMaster(map<int,vector<SingleMeshCell>> &Localmap,const int &Tag);

    //*************************************************************
    //*** for nodal phy group transfer
    //*************************************************************
    /**
     * Send the nodal physical name and nodes id to each rank
     * @param Phyname the physical group name to be sent
     * @param Nodesid the mesh cell vector to be sent
     * @param Tag the id of the sending message
     * @param Cpuid the rank id used to receive the message
    */
    static void sendPhyName2NodeIDVecMapToOthers(const string &Phyname,const vector<int> &Nodesid,const int &Tag,const int &Cpuid);
    /**
     * Receive the physical group name and mesh cell for local map
     * @param Localmap the local phyname-->nodeid map
     * @param Tag the id of the sending message
    */
    static void receivePhyName2NodeIDVecMapFromMaster(map<string,vector<int>> &Localmap,const int &Tag);

};