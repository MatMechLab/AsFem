//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.12.30
//+++ Purpose: here we call the TS+SNES solver from PETSc to solve
//+++          our nonlinear system equation R(x,xdot,t)->0
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "TimeStepping/TimeStepping.h"


bool TimeStepping::Solve(Mesh &mesh,DofHandler &dofHandler,
            ElmtSystem &elmtSystem,MateSystem &mateSystem,
            BCSystem &bcSystem,ICSystem &icSystem,
            SolutionSystem &solutionSystem,EquationSystem &equationSystem,
            FE &fe,FESystem &feSystem,
            OutputSystem &outputSystem,
            Postprocess &postprocessSystem,
            FEControlInfo &fectrlinfo,
            NonlinearSolver &nonlinearSolver){

    // initialize the feControl structure
    fectrlinfo.CurrentStep=1;
    fectrlinfo.dt=_Dt;
    fectrlinfo.t=0.0;

    char buff[68];
    string str;

    // apply the initial condition to the solution
    icSystem.ApplyIC(mesh,dofHandler,solutionSystem._U);
    VecCopy(solutionSystem._U,solutionSystem._Uold);
    VecCopy(solutionSystem._U,solutionSystem._Unew);
    VecCopy(solutionSystem._U,solutionSystem._Utemp);
    
    // initialize the history variables
    feSystem.FormBulkFE(FECalcType::InitMaterial,0.0,_Dt,fectrlinfo.ctan,mesh,dofHandler,fe,elmtSystem,mateSystem,solutionSystem,equationSystem._AMATRIX,equationSystem._RHS);
    solutionSystem.UpdateMaterials();

    // write result to the head of pvd file
    outputSystem.WritePVDFileHeader();
    if(fectrlinfo.IsProjection){
        feSystem.FormBulkFE(FECalcType::Projection,0.0,_Dt,fectrlinfo.ctan,mesh,dofHandler,fe,elmtSystem,mateSystem,solutionSystem,equationSystem._AMATRIX,equationSystem._RHS);
    }
    outputSystem.WriteResultToFile(0,mesh,dofHandler,solutionSystem);
    outputSystem.WriteResultToPVDFile(0.0,outputSystem.GetOutputFileName());
    MessagePrinter::PrintNormalTxt("Write result to "+outputSystem.GetOutputFileName());
    MessagePrinter::PrintDashLine();
   
    bool HasConvergeSolution; 
    for(double currenttime=0.0;currenttime<_FinalT;){
        HasConvergeSolution=false;
        while(fectrlinfo.dt>_DtMin){
            snprintf(buff,68,"TimeStepping: step=%8d,time=%12.5e,dt=%12.5e",fectrlinfo.CurrentStep,fectrlinfo.t+fectrlinfo.dt,fectrlinfo.dt);
            str=buff;
            MessagePrinter::PrintNormalTxt(str);
            if(nonlinearSolver.Solve(mesh,dofHandler,elmtSystem,mateSystem,bcSystem,solutionSystem,equationSystem,fe,feSystem,fectrlinfo)){
                // now the nonlinear solver converged 
                // the final solution is stored in Unew of solutionSystem
                // We first update the solution system
                VecCopy(solutionSystem._U,solutionSystem._Uold);
                VecCopy(solutionSystem._Unew,solutionSystem._U);
                VecCopy(solutionSystem._V,solutionSystem._Vold);

                // in FormBulkFE, we used Utemp
                VecCopy(solutionSystem._Unew,solutionSystem._Utemp);

                // then we update the history variables
                feSystem.FormBulkFE(FECalcType::UpdateMaterial,fectrlinfo.t,fectrlinfo.dt,fectrlinfo.ctan,mesh,dofHandler,fe,elmtSystem,mateSystem,solutionSystem,equationSystem._AMATRIX,equationSystem._RHS);
                // update the time and step
                currenttime+=fectrlinfo.dt;
                fectrlinfo.t+=fectrlinfo.dt;

                // now we do the projection
                if(fectrlinfo.IsProjection){
                    feSystem.FormBulkFE(FECalcType::Projection,fectrlinfo.t,fectrlinfo.dt,fectrlinfo.ctan,mesh,dofHandler,fe,elmtSystem,mateSystem,solutionSystem,equationSystem._AMATRIX,equationSystem._RHS);
                }
                // now, since the projection quantities are ready, we can do either the output or the postprocess
                if(fectrlinfo.CurrentStep%postprocessSystem.GetOutputIntervalNum()==0){
                    postprocessSystem.RunPostprocess(fectrlinfo.t,mesh,dofHandler,fe,solutionSystem);
                }
                // for result output
                if(fectrlinfo.CurrentStep%outputSystem.GetIntervalNum()==0){
                    outputSystem.WriteResultToFile(fectrlinfo.CurrentStep,mesh,dofHandler,solutionSystem);
                    outputSystem.WriteResultToPVDFile(fectrlinfo.t,outputSystem.GetOutputFileName());
                    MessagePrinter::PrintNormalTxt("Write result to "+outputSystem.GetOutputFileName(),MessageColor::BLUE);
                    MessagePrinter::PrintDashLine();
                }

                // update the step
                fectrlinfo.CurrentStep+=1;
                // update the materials
                solutionSystem.UpdateMaterials();

                // for adaptive time stepping
                if(IsAdaptive()){
                    if(nonlinearSolver.GetFinalInterations()<=_OptIters){
                        fectrlinfo.dt*=_GrowthFactor;
                        if(fectrlinfo.dt>_DtMax) fectrlinfo.dt=_DtMax;
                    }
                    else{
                        fectrlinfo.dt*=_CutBackFactor;
                        if(fectrlinfo.dt<_DtMin) fectrlinfo.dt=_DtMin;
                    }
                }

                HasConvergeSolution=true;
                // for converged case, we jump out the loop
                break;
            }// ===>end-of-if-nonlinearSolver-converged-case
            else{
                // nonlinearSolver diverged, then we will try to reduce delta t
                fectrlinfo.dt*=0.5;
                snprintf(buff,68,"TimeStepping failed: step=%8d,reduce dt to %14.5e",fectrlinfo.CurrentStep,fectrlinfo.dt);
                str=buff;
                MessagePrinter::PrintNormalTxt(str,MessageColor::RED);
                HasConvergeSolution=false;
            }//===>end-of-converge-diverge-if-condition
        }//===end-of-while-for-failed-try
        if(!HasConvergeSolution){
            MessagePrinter::PrintErrorTxt("The minimum delta t has been used for step="+to_string(fectrlinfo.CurrentStep)+", however, it still failed");
            MessagePrinter::PrintErrorTxt("Please check your code or your boundary/initial condition for the simulation");
            MessagePrinter::AsFem_Exit();
            break;
        }
    }
    outputSystem.WritePVDFileEnd();

    return true;
}
