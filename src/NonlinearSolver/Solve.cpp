//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "NonlinearSolver/NonlinearSolver.h"

bool NonlinearSolver::Solve(Mesh &mesh,DofHandler &dofHandler,
               ElmtSystem &elmtSystem,MateSystem &mateSystem,
               BCSystem &bcSystem,ICSystem &icSystem,
               Solution &solution,EquationSystem &equationSystem,
               FE &fe,FESystem &feSystem){
    return NewtonRaphson(mesh,dofHandler,
                         elmtSystem,mateSystem,
                         bcSystem,icSystem,
                         solution,equationSystem,fe,feSystem);
}