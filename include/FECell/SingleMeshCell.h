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
//+++ Date    : 2023.12.30
//+++ Function: defines the basic data structure of a single mesh cell
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <utility>

#include "FECell/MeshType.h"
#include "FECell/Nodes.h"

using std::vector;
using std::map;
using std::pair;
using std::string;

/**
 * The data structure of a single mesh cell, it stores the basic data info for a single mesh cell
 */
struct SingleMeshCell{
    int Dim;/**< the dimension number of current mesh cell */
    int NodesNumPerElmt;/**< nodes number per current mesh cell */
    vector<int> ElmtConn;/**< connectivity of current mesh cell */
    vector<int> ElmtDofIDs;/**< dof ids of current mesh cell */
    Nodes ElmtNodeCoords0;/**< node coordinates of current mesh cell (undeformed)*/
    Nodes ElmtNodeCoords ;/**< node coordinates of current mesh cell (deformed)*/
    int VTKCellType;/**< the vtk cell type of current mesh */ 
    MeshType CellMeshType;/**< the mesh type of current mesh cell */

    int            PhysicalGroupNums;/**< number of physical groups of current mesh cell, default one should be 1 */
    vector<string> PhysicalNameList;/**< pyhsical name list of current mesh cell, it can have more than one name */
    vector<int>    PhysicalIDList;/**< physical id list of current mesh cell, it can have more than one phy id */

};
