//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_EQUATIONSYSTEM_H
#define ASFEM_EQUATIONSYSTEM_H

#include <iostream>
#include <iomanip>
#include <string>

#include "petsc.h"

#include "DofHandler/DofHandler.h"

using namespace std;


class EquationSystem{
public:
    EquationSystem();

    void InitEquationSystem(const PetscInt &ndofs,const PetscInt &maxrownnz);
    void CreateSparsityPattern(DofHandler &dofHandler);

public:
    Mat _AMATRIX;
    Vec _RHS;

private:
    PetscInt _nDofs;
};



#endif // ASFEM_EQUATIONSYSTEM_H