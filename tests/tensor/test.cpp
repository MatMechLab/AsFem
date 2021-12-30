//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.12.30
//+++ Purpose: test cpp for the functionality of AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include <iostream>
#include "petsc.h"

#include "Utils/RankTwoTensor.h"
#include "Utils/RankFourTensor.h"

int main(int args,char *argv[]){
    PetscErrorCode ierr;
    ierr=PetscInitialize(&args,&argv,NULL,NULL);if (ierr) return ierr;


    RankFourTensor r4;
    r4.SetFromEandNu(120.0,0.3);

    r4.PrintVoigt();

    cout<<r4(4,4,4,4)<<endl;

    ierr=PetscFinalize();CHKERRQ(ierr);
    return ierr;
}
