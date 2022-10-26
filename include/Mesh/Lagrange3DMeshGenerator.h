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
//+++ Date   : 2022.05.04
//+++ Purpose: the lagrange 3d mesh generator
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "Mesh/MeshGeneratorBase.h"

/**
 * the mesh generator for 3d lagrange mesh
 */
class Lagrange3DMeshGenerator:public MeshGeneratorBase{
public:
    /**
     * constructor
     */
    Lagrange3DMeshGenerator();
    /**
     * deconstructor
     */
    ~Lagrange3DMeshGenerator();

    /**
     * function for the details of 3d lagrange mesh generation, if everything works fine, it will return true.
     * @param t_meshtype the type of mesh one want to use for the mesh generation
     * @param t_meshdata the mesh data structure, which should be updated within each mesh generator!
     */
    virtual bool generateMesh(const MeshType &t_meshtype,MeshData &t_meshdata) override;
private:
    bool m_mesh_generated=false;
    vector<vector<int>> leftconn,rightconn;
    vector<vector<int>> bottomconn,topconn;
    vector<vector<int>> backconn,frontconn;
    vector<int> leftnodes,rightnodes;
    vector<int> bottomnodes,topnodes;
    vector<int> backnodes,frontnodes;
};