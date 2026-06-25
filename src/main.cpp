//****************************************************************
//* This file is part of the AsFem framework
//* Advanced Simulation kit based on Finite Element Method (AsFem)
//* All rights reserved, Yang Bai/MM-Lab@CopyRight 2020-present
//* https://github.com/MatMechLab/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.12.28
//+++ Update : 2026.06.25
//+++ Purpose: the main program of the whole AsFem framework
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include <iostream>

#include "Application.h"
#include "FEProblem/FEProblem.h"

int main(int args,char *argv[]){
    Application myapp;

    if(myapp.init(args,argv)) return 1;

    const int Year=2026;
    const int Month=6;
    const int Day=25;
    const double Version=0.9;
    myapp.printAppInfo(Year,Month,Day,Version);


    FEProblem feProblem;
    feProblem.initFEProblem(args,argv);
    feProblem.run();
    feProblem.finalize();

    return myapp.finalize();
}
