//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.06.26
//+++ Purpose: define the base class for mesh import function
//+++          other meshio, i.e. gmsh and netgen should inherit
//+++          from this
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "Mesh/Mesh.h"

using namespace std;

class MeshIOBase{
public:
    // MeshIOBase();
    virtual bool ReadMeshFromFile(Mesh &mesh)=0;
    virtual void SetMeshFileName(string filename)=0;
    virtual string GetMeshFileName()const=0;

protected:
    string _MeshFileName="";
    bool _HasSetMeshFileName=false;
};