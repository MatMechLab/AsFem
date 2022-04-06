//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2022.04.06
//+++ Purpose: implement the hand-write newton-raphson method without
//+++          any line search from PETSc
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "NonlinearSolver/NonlinearSolver.h"

bool NonlinearSolver::NewtowRaphsonSolve(Mesh &mesh,DofHandler &dofHandler,
            ElmtSystem &elmtSystem,MateSystem &mateSystem,
            BCSystem &bcSystem,
            SolutionSystem &solutionSystem,EquationSystem &equationSystem,
            FE &fe,FESystem &feSystem,
            FEControlInfo &fectrlinfo){
    _monctx.iters=0;
    bool IsConverge=false;
    char buff[68];
    string str;
    while(_monctx.iters<_MaxIters && !IsConverge){
        // preset new boundary condition to unew
        bcSystem.ApplyPresetBC(mesh,dofHandler,FECalcType::ComputeResidual,
                                fectrlinfo.t+fectrlinfo.dt,
                                fectrlinfo.ctan,
                                solutionSystem._Unew,
                                equationSystem._AMATRIX,
                                equationSystem._RHS);
        // set the ctan array according to different time stepping method
        if(fectrlinfo._timesteppingtype==TimeSteppingType::STATIC){
            // for static analysis
            fectrlinfo.ctan[0]=1.0;
            fectrlinfo.ctan[1]=0.0;
            fectrlinfo.ctan[2]=0.0;
            VecCopy(solutionSystem._Unew,solutionSystem._Utemp);// in FormBulkFE, we always use Utemp for the calculation
        }
        else if(fectrlinfo._timesteppingtype==TimeSteppingType::BACKWARDEULER){
            fectrlinfo.ctan[0]=1.0;
            fectrlinfo.ctan[1]=1.0/fectrlinfo.dt;
            fectrlinfo.ctan[2]=0.0;
            // calculate the current velocity
            VecWAXPY(solutionSystem._V,-1.0,solutionSystem._U,solutionSystem._Unew);//V=-Uold+Unew
            VecScale(solutionSystem._V,fectrlinfo.ctan[1]);//V=V*1.0/dt
            VecCopy(solutionSystem._Unew,solutionSystem._Utemp);// Used in FormBulkFE
        }
        else if(fectrlinfo._timesteppingtype==TimeSteppingType::CRANCKNICLSON){
            fectrlinfo.ctan[0]=0.5;
            fectrlinfo.ctan[1]=1.0/fectrlinfo.dt;
            fectrlinfo.ctan[2]=0.0;
            // calculate the current velocity
            VecWAXPY(solutionSystem._V,-1.0,solutionSystem._U,solutionSystem._Unew);//V=-Uold+Unew
            VecScale(solutionSystem._V,fectrlinfo.ctan[1]);//V=V*1.0/dt

            // for Cranck-Nicolson Utemp=0.5*(Unew+U)
            VecWAXPY(solutionSystem._Utemp,1.0,solutionSystem._Unew,solutionSystem._U);// Utemp=Unew+Uold
            VecScale(solutionSystem._Utemp,0.5);// Utemp=0.5*(Unew+U)
        }
        else if(fectrlinfo._timesteppingtype==TimeSteppingType::BDF2){
            // the expression for BDF2 is:
            // [Un+2-(4/3)Un+1+(1/3)Un]/dt=(2/3)f(Un+1)

            if(fectrlinfo.CurrentStep<2){
                fectrlinfo.ctan[0]=1.0;
                fectrlinfo.ctan[1]=1.0/fectrlinfo.dt;
                fectrlinfo.ctan[2]=0.0;
                // calculate the current velocity
                VecWAXPY(solutionSystem._V,-1.0,solutionSystem._Unew,solutionSystem._U);//V=-Uold+Unew
                VecScale(solutionSystem._V,fectrlinfo.ctan[1]);//V=V*1.0/dt
            }
            else{
                fectrlinfo.ctan[0]=2.0/3.0;
                fectrlinfo.ctan[1]=1.0/fectrlinfo.dt;
                fectrlinfo.ctan[2]=0.0;
                // calculate the current velocity
                VecWAXPY(solutionSystem._V,-4.0/3.0,solutionSystem._U,solutionSystem._Unew);//V=Unew-(4/3)*U
                VecAXPY(solutionSystem._V,1.0/3.0,solutionSystem._Uold);// V=V+(1/3)Uold
                VecScale(solutionSystem._V,fectrlinfo.ctan[1]);//V=V*1.0/dt
            }
            // for Utemp used in FormBulkFE
            VecCopy(solutionSystem._Unew,solutionSystem._Utemp);
        }
        else{
            MessagePrinter::PrintErrorTxt("Unsupported time stepping method in nonlinear solver");
            MessagePrinter::AsFem_Exit();
        }

        feSystem.FormBulkFE(FECalcType::ComputeResidual,
                        fectrlinfo.t+fectrlinfo.dt,
                        fectrlinfo.dt,
                        fectrlinfo.ctan,
                        mesh,dofHandler,fe,
                        elmtSystem,mateSystem,
                        solutionSystem,
                        equationSystem._AMATRIX,equationSystem._RHS);
        
        feSystem.FormBulkFE(FECalcType::ComputeJacobian,
                        fectrlinfo.t+fectrlinfo.dt,
                        fectrlinfo.dt,
                        fectrlinfo.ctan,
                        mesh,dofHandler,fe,
                        elmtSystem,mateSystem,
                        solutionSystem,
                        equationSystem._AMATRIX,equationSystem._RHS);
    
        bcSystem.SetBCPenaltyFactor(feSystem.GetMaxAMatrixValue()*1.0e8);
        bcSystem.ApplyBC(mesh,dofHandler,fe,
                        FECalcType::ComputeResidual,
                        fectrlinfo.t+fectrlinfo.dt,
                        fectrlinfo.ctan,
                        solutionSystem._Unew,solutionSystem._V,
                        equationSystem._AMATRIX,equationSystem._RHS);
        bcSystem.ApplyBC(mesh,dofHandler,fe,
                        FECalcType::ComputeJacobian,
                        fectrlinfo.t+fectrlinfo.dt,
                        fectrlinfo.ctan,
                        solutionSystem._Unew,solutionSystem._V,
                        equationSystem._AMATRIX,equationSystem._RHS);

        MatScale(equationSystem._AMATRIX,-1.0);

        KSPSetOperators(_ksp,equationSystem._AMATRIX,equationSystem._AMATRIX);
        KSPSolve(_ksp,equationSystem._RHS,solutionSystem._Vold);

        // update solution and iteration info
        _monctx.iters+=1;
        // VecAXPY(Vec y,PetscScalar a,Vec x);y = y + a*x
        VecAYPX(solutionSystem._Unew,1.0,solutionSystem._Vold);

        VecNorm(equationSystem._RHS,NORM_2,&_monctx.rnorm);
        VecNorm(solutionSystem._Vold,NORM_2,&_monctx.dunorm);
        _monctx.enorm=_monctx.rnorm*_monctx.dunorm;
        if(_monctx.iters==1){
            _monctx.rnorm0=_monctx.rnorm;
            _monctx.dunorm0=_monctx.dunorm;
            _monctx.enorm0=_monctx.enorm;
        }
        if(fectrlinfo.IsDepDebug){
            snprintf(buff,68,"  AsFem NR solver: iters=%4d, |R|=%12.5e, |dU|=%12.5e",_monctx.iters,_monctx.rnorm,_monctx.dunorm);
            str=buff;
            MessagePrinter::PrintNormalTxt(str);
        }
        if(_monctx.rnorm<_RAbsTol || _monctx.rnorm<_RRelTol*_monctx.rnorm0){
            IsConverge=true;
            break;
        }
    }
    return IsConverge;
}