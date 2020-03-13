//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_FECTRLINFO_H
#define ASFEM_FECTRLINFO_H

#include <iostream>

#include "petsc.h"

#include "TimeStepping/TimeSteppingType.h"

using namespace std;

class FeCtrlInfo{
public:
    double ctan[2];
    double dt=1.0e-6;
    double t=1.0;
    PetscInt CurrentStep=0;
    PetscInt FinalStep=0;
    bool IsDebug=true;
    bool IsDepDebug=false;
    bool IsProjection=false;

    TimeSteppingType timesteppingtype=TimeSteppingType::BackWardEuler;
    
};

#endif // ASFEM_FECTRLINFO_H