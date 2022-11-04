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
//+++ Date   : 2022.05.06
//+++ Purpose: the mesh generator for AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "Mesh/Lagrange1DMeshGenerator.h"
#include "Mesh/Lagrange2DMeshGenerator.h"
#include "Mesh/Lagrange3DMeshGenerator.h"

#include "Utils/MessagePrinter.h"

/**
 * the mesh generator class for the buit-in mesh generation in AsFem
 */
class MeshGenerator:public Lagrange1DMeshGenerator,
                    public Lagrange2DMeshGenerator,
                    public Lagrange3DMeshGenerator{
public:
    /**
     * constructor
     */
    MeshGenerator();

    /**
     * built-in mesh generation, if everything works fine, it should return true
     * @param dim the dimension of mesh
     * @param meshtype the type of generated mesh
     * @param meshData the mesh data structure
     */
    bool createMesh(const int &dim,const MeshType &meshtype,MeshData &meshdata);

    /**
     * get the mesh generation status, true for success, false for failed case.
     */
    bool getMeshGeneratorStatus()const{return m_isMeshGenerated;}

private:
    bool m_isMeshGenerated;/**< boolean var for mesh generator status */
};