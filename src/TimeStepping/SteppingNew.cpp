//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "TimeStepping/TimeStepping.h"


void TimeStepping::SteppingNew(Mesh &mesh,DofHandler &dofHandler,
                BCSystem &bcSystem,ICSystem &icSystem,
                ElmtSystem &elmtSystem,MateSystem &mateSystem,
                EquationSystem &equationSystem,Solution &solution,
                FE &fe,FESystem &feSystem,
                OutputSystem &outputSystem,
                NonlinearSolver &nonlinearsolver,
                FeCtrlInfo &fectrl){
    if(_TimeSteppingType==TimeSteppingType::BackWardEuler){
        BackwardEuler(mesh,dofHandler,
                      bcSystem,icSystem,
                      elmtSystem,mateSystem,
                      equationSystem,solution,
                      fe,feSystem,
                      outputSystem,nonlinearsolver,fectrl);
    }
    else if(_TimeSteppingType==TimeSteppingType::CrankNicolson){
        CrankNicolson(mesh,dofHandler,
                    bcSystem,icSystem,
                    elmtSystem,mateSystem,
                    equationSystem,solution,
                    fe,feSystem,
                    outputSystem,nonlinearsolver,fectrl);
    }
}