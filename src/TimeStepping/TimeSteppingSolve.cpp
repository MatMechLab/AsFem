//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.12.29
//+++ Purpose: Do the time stepping for transient problem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "TimeStepping/TimeStepping.h"

bool TimeStepping::solve(Mesh &mesh,DofHandler &dofhandler,FE &fe,
                         ElmtSystem &elmtsystem,MateSystem &matesystem,
                         FESystem &fesystem,
                         BCSystem &bcsystem,
                         ICSystem &icsystem,
                         SolutionSystem &solutionsystem,
                         EquationSystem &equationsystem,
                         ProjectionSystem &projection,
                         FEControlInfo &fectrlinfo,
                         NonlinearSolver &nlsolver,
                         OutputSystem &output,
                         Postprocessor &postprocess){

    char buff[68];//77-12=65
    string str;
    int lastiters=1000;
    fectrlinfo.dt=m_data.m_dt0;
    fectrlinfo.CurrentStep=0;

    solutionsystem.m_u_current.setToZero();
    icsystem.applyInitialConditions(mesh,dofhandler,solutionsystem.m_u_current);
    solutionsystem.m_u_old.copyFrom(solutionsystem.m_u_current);
    solutionsystem.m_u_older.copyFrom(solutionsystem.m_u_current);

    // initialize the material
    fesystem.formBulkFE(FECalcType::INITMATERIAL,fectrlinfo.t,fectrlinfo.dt,fectrlinfo.ctan,
                        mesh,dofhandler,fe,
                        elmtsystem,matesystem,
                        solutionsystem,
                        equationsystem.m_amatrix,equationsystem.m_rhs);
    projection.executeProjection(mesh,dofhandler,elmtsystem,matesystem,fe,solutionsystem,fectrlinfo);
    MessagePrinter::printDashLine();
    MessagePrinter::printNormalTxt("Material properties have been initialized");
    MessagePrinter::printDashLine();
    
    output.savePVDHead();
    output.savePVDEnd();
    output.saveResults2File(0,mesh,dofhandler,solutionsystem,projection);
    output.savePVDResults(0.0);
    MessagePrinter::printDashLine(MessageColor::BLUE);
    MessagePrinter::printNormalTxt("Save results to "+output.getOutputFileName(),MessageColor::BLUE);
    MessagePrinter::printDashLine(MessageColor::BLUE);
    if(postprocess.hasPostprocess()){
        postprocess.prepareCSVFileHeader();
        postprocess.executePostprocess(mesh,dofhandler,fe,matesystem,projection,solutionsystem);
        postprocess.savePPSResults2CSVFile(0.0);
        MessagePrinter::printDashLine(MessageColor::BLUE);
        MessagePrinter::printNormalTxt("Save postprocess result to "+postprocess.getCSVFileName(),MessageColor::BLUE);
        MessagePrinter::printDashLine(MessageColor::BLUE);
    }
    MessagePrinter::printStars();

    for(fectrlinfo.t=0.0;fectrlinfo.t<=m_data.m_finaltime;){
        snprintf(buff,68,"Time=%13.5e, step=%8d, dt=%13.5e",fectrlinfo.t+fectrlinfo.dt,fectrlinfo.CurrentStep+1,fectrlinfo.dt);
        str=buff;
        MessagePrinter::printNormalTxt(str);
        if(nlsolver.solve(mesh,dofhandler,fe,
                          elmtsystem,matesystem,fesystem,
                          bcsystem,solutionsystem,equationsystem,
                          fectrlinfo)){
            // if the current nonlinear process success
            fectrlinfo.t+=fectrlinfo.dt;
            fectrlinfo.CurrentStep+=1;

            // update the material properties
            // update the solution
            solutionsystem.m_u_temp.copyFrom(solutionsystem.m_u_current);
            fesystem.formBulkFE(FECalcType::UPDATEMATERIAL,fectrlinfo.t,fectrlinfo.dt,fectrlinfo.ctan,
                                mesh,dofhandler,fe,
                                elmtsystem,matesystem,
                                solutionsystem,
                                equationsystem.m_amatrix,equationsystem.m_rhs);


            if(fectrlinfo.CurrentStep%output.getIntervalNum()==0){
                projection.executeProjection(mesh,dofhandler,elmtsystem,matesystem,fe,solutionsystem,fectrlinfo);
                output.saveResults2File(fectrlinfo.CurrentStep,mesh,dofhandler,solutionsystem,projection);
                output.savePVDResults(fectrlinfo.t);
                MessagePrinter::printDashLine(MessageColor::BLUE);
                MessagePrinter::printNormalTxt("Save results to "+output.getOutputFileName(),MessageColor::BLUE);
                MessagePrinter::printDashLine(MessageColor::BLUE);
            }
            if(postprocess.hasPostprocess()){
                if(fectrlinfo.CurrentStep%postprocess.getInterval()==0){
                    if(fectrlinfo.CurrentStep%output.getIntervalNum()!=0){
                        projection.executeProjection(mesh,dofhandler,elmtsystem,matesystem,fe,solutionsystem,fectrlinfo);
                    }
                    postprocess.executePostprocess(mesh,dofhandler,fe,matesystem,projection,solutionsystem);
                    postprocess.savePPSResults2CSVFile(fectrlinfo.t);
                    MessagePrinter::printDashLine(MessageColor::BLUE);
                    MessagePrinter::printNormalTxt("Save postprocess result to "+postprocess.getCSVFileName(),MessageColor::BLUE);
                    MessagePrinter::printDashLine(MessageColor::BLUE);
                }
            }
            MessagePrinter::printStars();
            if(isAdaptive()){
                if(nlsolver.getIterationNum()<=getOptimizeIters() && lastiters<=getOptimizeIters()){
                    fectrlinfo.dt*=getGrowthFactor();
                    if(fectrlinfo.dt>getMaxDt()) fectrlinfo.dt=getMaxDt();
                }
                else if(nlsolver.getIterationNum()<=getOptimizeIters() && lastiters>getOptimizeIters()){
                    // do not change current dt
                    fectrlinfo.dt*=1.0;
                }
                else if(nlsolver.getIterationNum()>getOptimizeIters()){
                    fectrlinfo.dt*=getCutbackFactor();
                    if(fectrlinfo.dt<getMinDt()) fectrlinfo.dt=getMinDt();
                }
            }
            // update the solution
            solutionsystem.m_u_older.copyFrom(solutionsystem.m_u_old);
            solutionsystem.m_u_old.copyFrom(solutionsystem.m_u_current);
            solutionsystem.m_v.setToZero();
            // store the previous step's iteration numbers
            lastiters=nlsolver.getIterationNum();
        }
        else{
            //  if the nonlinear solver failed we will try to reduce the dt
            fectrlinfo.dt*=getCutbackFactor();
            snprintf(buff,68," Transient solver failed, reduce dt to %13.5e",fectrlinfo.dt);
            str=buff;
            MessagePrinter::printWarningTxt(str);
            if(fectrlinfo.dt<getMinDt()){
                MessagePrinter::printErrorTxt("The minimum detal t is reached, however, your solver still fails. Please check either your code or your boundary conditions");
                return false;
            }
        }
    }
    if(fectrlinfo.CurrentStep%output.getIntervalNum()!=0){
        projection.executeProjection(mesh,dofhandler,elmtsystem,matesystem,fe,solutionsystem,fectrlinfo);
        output.saveResults2File(fectrlinfo.CurrentStep,mesh,dofhandler,solutionsystem,projection);
        output.savePVDResults(fectrlinfo.t);
        MessagePrinter::printDashLine(MessageColor::BLUE);
        MessagePrinter::printNormalTxt("Save results to "+output.getOutputFileName(),MessageColor::BLUE);
        MessagePrinter::printDashLine(MessageColor::BLUE);
        MessagePrinter::printStars();
    }
    
    return true;
}