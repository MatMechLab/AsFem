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
//+++ Date    : 2024.01.06
//+++ Function: finite element cell data structure
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once


#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <map>

#include "FECell/SingleMeshCell.h"


using std::vector;
using std::map;
using std::make_pair;
using std::string;

/**
 * This structure defines the basic data structure of the fe cell, i.e., mesh connectivity, dof maps, etc.
*/
struct FECellData{
    int MeshOrder;/**< the order of the fe cell mesh */
    int MaxDim;/**< the maximum dimension of the fe cell mesh */
    int MinDim;/**< the minimal dimension of the fe cell mesh */

    int NodesNum;/**< the total nodes number of the fe cell mesh */
    int NodesNumPerLineElmt;/**< nodes number per line element */
    int NodesNumPerSurfElmt;/**< nodes number per surf element */
    int NodesNumPerBulkElmt;/**< nodes number per bulk element */

    int ElmtsNum;/**< the elements number of the fe cell mesh */
    int BulkElmtsNum;/**< the bulk elements number of the fe cell mesh */
    int SurfElmtsNum;/**< the surface elements number of the fe cell mesh */
    int LineElmtsNum;/**< the line elements number of the fe cell mesh */
    string BulkMeshTypeName;/**< string name of the mesh type of bulk cell */

    int BulkElmtVTKCellType;/**< bulk element's vtk cell type */
    int SurfElmtVTKCellType;/**< surf element's vtk cell type */
    int LineElmtVTKCellType;/**< line element's vtk cell type */

    MeshType BulkElmtMeshType;/**< bulk element's mesh type */
    MeshType SurfElmtMeshType;/**< surf element's mesh type */
    MeshType LineElmtMeshType;/**< line element's mesh type */

    double Xmin;/**< xmin */
    double Xmax;/**< xmax */
    double Ymin;/**< ymin */
    double Ymax;/**< ymax */
    double Zmin;/**< zmin */
    double Zmax;/**< zmax */
    int Nx;/**< elements num along x-axis */
    int Ny;/**< elements num along y-axis */
    int Nz;/**< elements num along z-axis */

    vector<double> NodeCoords_Global;/**< the nodal coordinates of all the mesh point, which is only stored in master rank */

    int MaxDofsPerNode;/**< the maximum dofs number of each node */
    int TotalDofsNum;/**< the total dofs number of the fe cell mesh */
    int ActiveDofsNum;/**< the active dofs number of the fe cell mesh */

    //*** for elemental physical group
    int PhyGroupNum_Global;/**< physical group numbers among all ranks */
    vector<int> PhyDimVector_Global;/**< physical group dimension vector among all ranks */
    vector<int> PhyIDVector_Global;/**< physical group id vector among all ranks */
    vector<string> PhyNameVector_Global;/**< physical name vector among all ranks */
    vector<int> PhyGroupElmtsNumVector_Global;/**< physical group's elements number vector among all ranks */
    map<int,string> PhyID2NameMap_Global;/**< physical id to name map among all ranks */
    map<string,int> PhyName2IDMap_Global;/**< physical name to id map among all ranks */

    vector<SingleMeshCell> MeshCell_Total;/**< this vector stores the whole mesh info (only active on master rank) */
    vector<SingleMeshCell> MeshCell_Local;/**< local mesh cell vector of each rank */

    vector<int> NodeIDs_Local;/**< this vector stores the unique nodal index for each rank */

    map<string,vector<SingleMeshCell>> PhyName2MeshCellVectorMap_Global;/**< the global phy name to cell vector map */
    map<int,vector<SingleMeshCell>>    PhyID2MeshCellVectorMap_Global;/**< the global phy id to cell vector map */

    map<string,vector<SingleMeshCell>> PhyName2MeshCellVectorMap_Local;/**< the local phy name to cell vector map */
    map<int,vector<SingleMeshCell>>    PhyID2MeshCellVectorMap_Local;/**< the local phy id to cell vector map */

    map<string,vector<int>> PhyName2BulkFECellIDMap_Global;/**< the global phy anme to cell id vector map */
    map<int,vector<int>>    PhyID2BulkFECellIDMap_Global;/**< the global phy id to cell id vector map */
    
    //*** for nodal physical group information
    int NodalPhyGroupNum_Global;/**< nodal physical group numbers among all ranks */
    vector<int> NodalPhyIDVector_Global;/**< nodal physical group id vector among all ranks */
    vector<string> NodalPhyNameVector_Global;/**< nodal physical name vector among all ranks */
    vector<int> NodalPhyGroupNodesNumVector_Global;/**< nodal physical group's elements number vector among all ranks */
    map<int,string> NodalPhyID2NameMap_Global;/**< nodal physical id to name map among all ranks */
    map<string,int> NodalPhyName2IDMap_Global;/**< nodal physical name to id map among all ranks */

    map<string,vector<int>> NodalPhyName2NodeIDVecMap_Global;/**< the global nodal phy name to node id map*/
    map<string,vector<int>> NodalPhyName2NodeIDVecMap_Local;/**< the global nodal phy name to node id map*/

    //*** for partition message
    vector<int> BulkCellPartionInfo_Global;/**< the partition info of the bulk fe cell */
    vector<int> RanksElmtsNum_Global;/**< this vector stores the number of bulk elmts owned by each rank */

};