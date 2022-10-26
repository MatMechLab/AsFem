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
//+++ Purpose: the base class for mesh generator
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <algorithm>

#include "Mesh/MeshData.h"
#include "Mesh/MeshType.h"

using std::sort;
using std::make_pair;

/**
 * this is the base class for mesh generator, all the generator should inherit this base class
 */
class MeshGeneratorBase{
public:
    /**
     * constructor
     */
    MeshGeneratorBase(){}

    /**
     * virtual function for the details of different mesh generation, the child class should
     * offer the implementations
     * @param t_meshtype the type of mesh one want to use for the mesh generation
     * @param t_meshdata the mesh data structure, which should be updated within each mesh generator!
     */
    virtual bool generateMesh(const MeshType &t_meshtype,MeshData &t_meshdata)=0;

};