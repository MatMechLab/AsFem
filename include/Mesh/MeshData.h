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
//+++ Date   : 2022.04.20
//+++ Purpose: this class defines the basic data structure of FEM mesh
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <utility>

#include "Mesh/MeshType.h"

using std::vector;
using std::map;
using std::pair;
using std::string;

/**
 * the data structure for mesh, it stores the nodal coordinates,
 * elemental connectivity, boundary element information, physical group information, and so on.
 */
struct MeshData{
    // for the geometry information of your mesh
    int m_maxdim;/**< the maximum dimension of your mesh */
    int m_mindim;/**< the minimum dimension of your mesh */
    int m_nx;/**< number of element in x-direction */
    int m_ny;/**< number of element in y-direction */
    int m_nz;/**< number of element in z-direction */
    double m_xmin;/**< x-min of the regular domain */
    double m_xmax;/**< x-max of the regular domain */
    double m_ymin;/**< y-min of the regular domain */
    double m_ymax;/**< y-max of the regular domain */
    double m_zmin;/**< z-min of the regular domain */
    double m_zmax;/**< z-max of the regular domain */
    int m_order;/**< order of the mesh */
    vector<double> m_nodecoords0;/**< vector for the coordinates of nodes, undeformed ! */
    vector<double> m_nodecoords;/**< vector for the coordinates of nodes, deformed one! */
    vector<vector<int>> m_bulkelmt_connectivity;/**< stores the connectivity of bulk elements */
    vector<double>      m_bulkelmt_volume;/**< stores the volume of each bulk element */
    vector<vector<int>> m_pointelmt_connectivity;/**< stores the connectivity of line elements */
    vector<double>      m_pointelmt_volume;/**< stores the volume of each line element */
    vector<vector<int>> m_lineelmt_connectivity;/**< stores the connectivity of line elements */
    vector<double>      m_lineelmt_volume;/**< stores the volume of each line element */
    vector<vector<int>> m_surfaceelmt_connectivity;/**< stores the connectivity of line elements */
    vector<double>      m_surfaceelmt_volume;/**< stores the volume of each surface element */
    // for the number of different nodes
    int m_nodes;/**< the total number of nodes */
    int m_nodesperbulkelmt;/**< the nodes number per bulk element */
    int m_nodesperlineelmt;/**< the nodes number per line element */
    int m_nodespersurfaceelmt;/**< the nodes number per surface element */
    // for the number of different elements
    int m_elements;/**< the number of total elements */
    int m_bulkelmts;/**< the number of bulk elements */
    int m_pointelmts;/**< the number of point elements */
    int m_lineelmts;/**< the number of line elements */
    int m_surfaceelmts;/**< the number of surface elements */
    // for physical group information of elements
    int m_phygroups;/**< number of total physical groups */
    vector<int>                              m_phygroup_dimvec;/**< the vector for the dim of each phy group */
    vector<string>                           m_phygroup_phynamevec;/**< the vector stores the name of each phy group */
    vector<int>                              m_phygroup_phyidvec;/**< the vector of physical id for each phy group */
    vector<int>                              m_phygroup_elmtnumvec;/**< the vector for the element number of each phy group */
    vector<int>                              m_phygroup_nodesnumperelmtvec;/**< the vector for the nodes number of the element in each phy group */
    vector<pair<string,vector<vector<int>>>> m_phygroup_name2elmtconnvec;/**< vector for the name to element connectivity map */
    vector<pair<string,vector<int>>>         m_phygroup_name2bulkelmtidvec;/** vector for the name to bulk element id map */
    vector<pair<string,int>>                 m_phygroup_name2dimvec;/**< vector for the name to elmt dim map */
    vector<pair<string,int>>                 m_phygroup_name2phyidvec;/**< vector for the name to phyid map */
    vector<pair<int,string>>                 m_phygroup_phyid2namevec;/**< vector for the phyid to name map */

    int m_nodal_phygroups;/**< number of nodal physical groups */
    vector<pair<string,vector<int>>> m_nodephygroup_name2nodeidvec;/**< vector for the name to node ids map */
    vector<pair<string,int>>         m_nodephygroup_name2phyidvec;/**< vector for the name to phyid map */
    vector<pair<int,string>>         m_nodephygroup_phyid2namevec;/**< vector for the phyid to name map */
    vector<string>                   m_nodephygroup_phynamevec;/**< the vector stores the name of node phy group */
    vector<int>                      m_nodephygroup_phyidvec;/**< the vector of physical id for node phy group */

    MeshType m_bulkelmt_type;/**< the meshtype of bulk element */
    MeshType m_lineelmt_type;/**< the meshtype of line element */
    MeshType m_surfaceelmt_type;/**< the meshtype of surface element */

    int m_bulkelmt_vtktype;/**< the vtk cell type of bulk element */
    int m_lineelmt_vtktype;/**< the vtk cell type of line element */
    int m_surfaceelmt_vtktype;/**< the vtk cell type of surface element */
    string m_bulkelmt_typename;/**< name of the bulk element type */

};
