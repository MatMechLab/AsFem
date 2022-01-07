//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/M3Group/AsFem
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

#include "Mesh/LagrangeMesh.h"

class Mesh:public LagrangeMesh{
public:
    Mesh();

    //************************************************************
    //*** for the basic getting
    //************************************************************
    int GetDim() const{return GetBulkMeshDim();}
    bool CreateMesh(){return LagrangeMesh::CreateLagrangeMesh();}
    void SaveMesh(string filename="")const{LagrangeMesh::SaveLagrangeMesh(filename);}


    void PrintMeshInfo()const{PrintBulkMeshInfo();}
    void PrintMeshDetailInfo()const{PrintBulkMeshInfoDetails();}

private:
    bool _HasMeshCreated=false;

};