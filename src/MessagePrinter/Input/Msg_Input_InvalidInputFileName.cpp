//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************


#include "MessagePrinter/MessagePrinter.h"

void Msg_Input_InvalidInputFileName(){
    PetscPrintf(PETSC_COMM_WORLD,"*** Error: invalid input file name or format or ext name !!!          ***\n");
    PetscPrintf(PETSC_COMM_WORLD,"***        your input file must look like 'xxx.i' end with '.i'       ***\n");
    // PetscPrintf(PETSC_COMM_WORLD,"*************************************************************************\n");
}