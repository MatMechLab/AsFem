//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_ICSYSTEM_H
#define ASFEM_ICSYSTEM_H

#include <iostream>
#include <iomanip>
#include <string>

#include "petsc.h"

//***********************************
//*** For AsFem's own header file
//***********************************
#include "MessagePrinter/MessagePrinter.h"

#include "ICType.h"
#include "ICBlock.h"

#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"

using namespace std;


class ICSystem{
public:
    ICSystem();
    void ApplyIC(Mesh &mesh,DofHandler &dofHandler,Vec &U);
    //*********************************************
    //*** add ICBlock from input file
    //*********************************************
    void AddICBlock(ICBlock &icblock);

    //*********************************************
    //*** some getting functions
    //*********************************************
    inline PetscInt GetICBlocksNum()const{return _nICBlocks;}
    ICBlock GetIthICBlock(const PetscInt &i)const{
        return _ICBlockList[i-1];
    }

    void PrintICSystemInfo()const;

private:
    //*******************************
    //*** Apply related ic
    //*******************************
    void ApplyConstIC(const vector<double> Params,const int &DofIndex,Mesh &mesh,DofHandler &dofHandler,Vec &U);
    void ApplyRandomIC(const vector<double> Params,const int &DofIndex,Mesh &mesh,DofHandler &dofHandler,Vec &U);
    void ApplyCircleIC(const vector<double> Params,const int &DofIndex,Mesh &mesh,DofHandler &dofHandler,Vec &U);

private:
    vector<ICBlock> _ICBlockList;
    PetscInt _nICBlocks;

    PetscMPIInt _rank,_size;
    PetscRandom _rnd;
};


#endif // ASFEM_ICSYSTEM_H