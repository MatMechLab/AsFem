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
//+++ Date   : 2020.06.28
//+++ Purpose: Define the final mesh class, so, all the lagrange
//+++          mesh, interface mesh, as well as the IGA mesh class
//+++          should be inherited here
//+++          For the newly implemented mesh class, it should also
//+++          be inherited here
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "Mesh/BulkMesh.h"

class Mesh:public BulkMesh{
public:
    Mesh();

    //************************************************************
    //*** for the basic getting
    //************************************************************
    int GetDim() const{return GetBulkMeshDim();}

private:
    bool _HasMeshCreated=false;

};