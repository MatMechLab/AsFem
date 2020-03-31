//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "NonlinearSolver/NonlinearSolver.h"


bool NonlinearSolver::NewtonRaphson(Mesh &mesh,DofHandler &dofHandler,
    ElmtSystem &elmtSystem,MateSystem &mateSystem,
    BCSystem &bcSystem,ICSystem &icSystem,
    Solution &solution,EquationSystem &equationSystem,
    FE &fe,FESystem &feSystem){
    _Iters=0;_IsConvergent=false;

    double t=1.0,dt=1.0;
    double ctan[]={1.0,1.0};

    if(icSystem.GetICBlocksNum()){}

    icSystem.ApplyIC(mesh,dofHandler,solution._Uold);
    bcSystem.ApplyInitialBC(mesh,dofHandler,t,solution._Uold);
    
    VecCopy(solution._Uold,solution._Unew);

    // VecView(solution._Uold,PETSC_VIEWER_STDOUT_WORLD);
    

    while(_Iters<_MaxIters&& !_IsConvergent){
        //w = a ∗ x + y
        VecWAXPY(solution._V,-1.0,solution._Uold,solution._Unew);//V=-Uold+Unew
        
        VecScale(solution._V,dt);//V=V*1.0/dt

        
        
        feSystem.FormFE(6,t,dt,ctan,mesh,dofHandler,
                        fe,elmtSystem,
                        mateSystem,
                        solution._Unew,solution._V,
                        solution._Hist,solution._HistOld,
                        solution._Proj,
                        equationSystem._AMATRIX,equationSystem._RHS);
        
        // VecView(equationSystem._RHS,PETSC_VIEWER_STDOUT_WORLD);

        bcSystem.SetBCPenaltyFactor(feSystem.GetMaxAMatrixValue());
        
        bcSystem.ApplyBC(mesh,dofHandler,fe,t,ctan,
                         equationSystem._AMATRIX,equationSystem._RHS,solution._Unew);
        
        // MatView(equationSystem._AMATRIX,PETSC_VIEWER_STDOUT_WORLD);

        
        if(LinearSolve(equationSystem._AMATRIX,solution._dU,equationSystem._RHS)){
            VecNorm(equationSystem._RHS,NORM_2,&_Rnorm);
            VecNorm(solution._dU,NORM_2,&_dUnorm);;
            _Enorm=_Rnorm*_dUnorm;
            if(_Iters==0){
                _Rnorm0=_Rnorm;
                _Enorm0=_Enorm;
                _dUnorm0=_dUnorm;
            }
            //update solution
            //y = y + a ∗ x//VecAXPY(Vec y,PetscScalar a,Vec x);
            VecAXPY(solution._Unew,1.0,solution._dU);

            _Iters+=1;

            PrintIterationDetailsInfo();

            if(CheckConvergence()){
                _IsConvergent=true;
                break;
            }
        }
        else{
            _IsConvergent=false;
            PetscPrintf(PETSC_COMM_WORLD,"*** Error: solve Ax=F failed in Newton-Raphson iteration        !!!   ***\n");
            PetscPrintf(PETSC_COMM_WORLD,"***        please check your ksp settings or your model         !!!   ***\n");
            return false;
        }
    }
    PrintIterationInfo();
    return _IsConvergent;
}
