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
//+++ Date   : 2020.07.10
//+++ Purpose: Define the boundary condition system in AsFem
//+++          Here we can apply:
//+++               1) constant value ic
//+++               2) random value ic
//+++               3) other type or user defined ic
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <string>
#include <vector>


//******************************************
//*** for AsFem own header
//******************************************
#include "Utils/MessagePrinter.h"

#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"


#include "ICSystem/ICBlock.h"
#include "ICSystem/ICType.h"

using namespace std;

class ICSystem{
public:
    ICSystem();
    void ApplyIC(const Mesh &mesh,const DofHandler &dofHandler,Vec &U);

    //******************************************************
    //*** add ic block
    //******************************************************
    void AddICBlock2List(ICBlock &icblock);
    inline int GetICBlocksNum()const{return _nICBlocks;}
    inline vector<ICBlock> GetICBlockList()const {return _ICBlockList;}


    void PrintICSystemInfo()const;

private:
    int _nICBlocks;
    vector<ICBlock> _ICBlockList;

    PetscMPIInt _rank,_size;
    PetscRandom _rnd;

private:
    //****************************************************
    //*** Apply different initial conditions
    //****************************************************
    void RunICLibs(const ICType &ictype,const int &DofIndex,const vector<double> &Parameters,const vector<string> &DomainList,const Mesh &mesh,const DofHandler &dofHandler,Vec &U);
    // for built-in ic
    void ApplyConstantIC(const int &DofIndex,const vector<double> &Parameters,const vector<string> &DomainList,const Mesh &mesh,const DofHandler &dofHandler,Vec &U);
    void ApplyRandomIC(const int &DofIndex,const vector<double> &Parameters,const vector<string> &DomainList,const Mesh &mesh,const DofHandler &dofHandler,Vec &U);
    void ApplyRectangleIC(const int &DofIndex,const vector<double> &Parameters,const vector<string> &DomainList,const Mesh &mesh,const DofHandler &dofHandler,Vec &U);
    void ApplyCircleIC(const int &DofIndex,const vector<double> &Parameters,const vector<string> &DomainList,const Mesh &mesh,const DofHandler &dofHandler,Vec &U);
    void ApplySmoothCircleIC(const int &DofIndex,const vector<double> &Parameters,const vector<string> &DomainList,const Mesh &mesh,const DofHandler &dofHandler,Vec &U);
    void ApplyCubicIC(const int &DofIndex,const vector<double> &Parameters,const vector<string> &DomainList,const Mesh &mesh,const DofHandler &dofHandler,Vec &U);
    void ApplySphereIC(const int &DofIndex,const vector<double> &Parameters,const vector<string> &DomainList,const Mesh &mesh,const DofHandler &dofHandler,Vec &U);
    // for user-defined ic
    void User1IC(const int &DofIndex,const vector<double> &Parameters,const vector<string> &DomainList,const Mesh &mesh,const DofHandler &dofHandler,Vec &U);
};
