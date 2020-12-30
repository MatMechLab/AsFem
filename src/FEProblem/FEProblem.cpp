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
//+++ Date   : 2020.06.29
//+++ Purpose: Define the FE problem analysis class in AsFem,
//+++          It is designed as the top level of the whole AsFem
//+++          framework
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FEProblem/FEProblem.h"

FEProblem::FEProblem(){
    _feJobType=FEJobType::STATIC;
}

void FEProblem::InitFEProblem(int args,char *argv[]){
    _feJobType=FEJobType::STATIC;
    _inputSystem.InitInputSystem(args,argv);
}


void FEProblem::Finalize(){
    _solutionSystem.ReleaseMem();
    _equationSystem.ReleaseMem();
    _nonlinearSolver.ReleaseMem();
    _timestepping.ReleaseMem();
}