//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.06.26
//+++ Purpose: implement a general mesh io class, which can import
//+++          mesh file from gmsh, netgen, etc...
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "Mesh/Gmsh2IO.h"
#include "Mesh/Gmsh4IO.h"
#include "Mesh/AbaqusIO.h"


enum class MeshIOType{
    NULLTYPE,
    GMSH2,
    GMSH4,
    NETGEN,
    ABAQUS,
    VTU,
    VTK
};

class MeshIO:public Gmsh2IO,public Gmsh4IO, public AbaqusIO{
public:
    MeshIO();
    virtual bool ReadMeshFromFile(Mesh &mesh) override;
    virtual void SetMeshFileName(string filename) override;
    virtual string GetMeshFileName()const override;

private:
    MeshIOType _MeshIOType=MeshIOType::GMSH2;

    bool IsGmsh2MeshFile(string filename);
    bool IsGmsh4MeshFile(string filename);
    bool IsNetgenMeshFile(string filename);
    bool IsAbaqusMeshFile(string filename);

    ifstream _in;
};