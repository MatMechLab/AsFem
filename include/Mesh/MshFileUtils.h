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
//+++ Date   : 2022.08.16
//+++ Purpose: This class offers the utils function for msh file
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "Utils/MessagePrinter.h"
#include "Mesh/MeshType.h"

/**
 * This class offers the utils function for msh file reading.
 */
class MshFileUtils{
public:
    /**
     * constructor
     */
    MshFileUtils();

    /**
     * get the current element's nodes number from the given element type(integer)
     * @param elmttype the (gmsh) element type
     */
    static int getElmtNodesNumFromElmtType(const int &elmttype);
    /**
     * get the current element's dimension from the given element type(integer)
     * @param elmttype the (gmsh) element type
     */
    static int getElmtDimFromElmtType(const int &elmttype);
    /**
     * get the current element's order from the given element type(integer)
     * @param elmttype the (gmsh) element type
     */
    static int getElmtOrderFromElmtType(const int &elmttype);
    /**
     * get the current element's vtk cell type from the given element type(integer)
     * @param elmttype the (gmsh) element type
     */
    static int getElmtVTKCellTypeFromElmtType(const int &elmttype);
    /**
     * get the current element's mesh type from the given element type(integer)
     * @param elmttype the (gmsh) element type
     */
    static MeshType getElmtMeshTypeFromElmtType(const int &elmttype);
    /**
     * get the current element's mesh type name from the given element type(integer)
     * @param elmttype the (gmsh) element type
     */
    static string getElmtMeshTypeNameFromElmtType(const int &elmttype);

    //********************************************
    //*** for lower dimension element
    //********************************************
    /**
     * get the current element's sub-element's nodes number from the given element type(integer)
     * @param elmttype the (gmsh) element type
     */
    static int getSubElmtNodesNumFromElmtType(const int &elmttype);
    /**
     * get the current element's sub-element's dimension from the given element type(integer)
     * @param elmttype the (gmsh) element type
     */
    static int getSubElmtDimFromElmtType(const int &elmttype);
    /**
     * get the current element's sub-element's order from the given element type(integer)
     * @param elmttype the (gmsh) element type
     */
    static int getSubElmtOrderFromElmtType(const int &elmttype);
    /**
     * get the current element's sub-element's mesh type from the given element type(integer)
     * @param elmttype the (gmsh) element type
     */
    static MeshType getSubElmtMeshTypeFromElmtType(const int &elmttype);
    //********************************************
    //*** for node reordering
    //********************************************
    /**
     * re-ordering gmsh nodes, make it consistent with AsFem
     * @param elmttype the (gmsh) element type
     * @param elmtconn the local element connectivity
     */
    static void reorderNodesIndex(const int &elmttype,vector<int> &elmtconn);
    /**
     * re-ordering gmsh nodes, make it consistent with AsFem
     * @param elmttype the (gmsh) element type
     * @param elmtconn the local element connectivity
     */
    static void reorderGmsh2NodesIndex(const int &elmttype,vector<int> &elmtconn);

};