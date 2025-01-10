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
//+++ Date   : 2020.12.29
//+++ Purpose: Do the time stepping for transient problem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "TimeStepping/TimeStepping.h"

bool TimeStepping::solve(FECell &fecell,DofHandler &dofhandler,FE &fe,
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
    fectrlinfo.Dt=m_Data.m_Dt0;
    fectrlinfo.CurrentStep=0;

    solutionsystem.m_Ucurrent.setToZero();
    icsystem.applyInitialConditions(fecell,dofhandler,solutionsystem.m_Ucurrent);
    solutionsystem.m_Uold.copyFrom(solutionsystem.m_Ucurrent);
    solutionsystem.m_Uolder.copyFrom(solutionsystem.m_Ucurrent);

    // initialize the material
    fesystem.formBulkFE(FECalcType::INITMATERIAL,fectrlinfo.T,fectrlinfo.Dt,fectrlinfo.Ctan,
                        fecell,dofhandler,fe,
                        elmtsystem,matesystem,
                        solutionsystem,
                        equationsystem.m_AMATRIX,equationsystem.m_RHS);
    projection.executeProjection(fecell,dofhandler,elmtsystem,matesystem,fe,solutionsystem,fectrlinfo);
    MessagePrinter::printDashLine();
    MessagePrinter::printNormalTxt("Material properties have been initialized");
    MessagePrinter::printDashLine();
    
    output.savePVDHead();
    output.savePVDEnd();
    output.saveResults2File(0,fecell,dofhandler,solutionsystem,projection);
    output.savePVDResults(0.0);
    MessagePrinter::printDashLine(MessageColor::BLUE);
    MessagePrinter::printNormalTxt("Save results to "+output.getOutputFileName(),MessageColor::BLUE);
    MessagePrinter::printDashLine(MessageColor::BLUE);
    if(postprocess.hasPostprocess()){
        postprocess.prepareCSVFileHeader();
        postprocess.executePostprocess(fecell,dofhandler,fe,matesystem,projection,solutionsystem);
        postprocess.savePPSResults2CSVFile(0.0);
        MessagePrinter::printDashLine(MessageColor::BLUE);
        MessagePrinter::printNormalTxt("Save postprocess result to "+postprocess.getCSVFileName(),MessageColor::BLUE);
        MessagePrinter::printDashLine(MessageColor::BLUE);
    }
    MessagePrinter::printStars();

    for(fectrlinfo.T=0.0;fectrlinfo.T<=m_Data.m_FinalTime;){
        snprintf(buff,68,"Time=%13.5e, step=%8d, dt=%13.5e",fectrlinfo.T+fectrlinfo.Dt,fectrlinfo.CurrentStep+1,fectrlinfo.Dt);
        str=buff;
        MessagePrinter::printNormalTxt(str);
        if(nlsolver.solve(fecell,dofhandler,fe,
                          elmtsystem,matesystem,fesystem,
                          bcsystem,solutionsystem,equationsystem,
                          fectrlinfo)){
            // if the current nonlinear process success
            fectrlinfo.T+=fectrlinfo.Dt;
            fectrlinfo.CurrentStep+=1;

            // update the material properties
            // update the solution
            solutionsystem.m_Utemp.copyFrom(solutionsystem.m_Ucurrent);
            fesystem.formBulkFE(FECalcType::UPDATEMATERIAL,fectrlinfo.T,fectrlinfo.Dt,fectrlinfo.Ctan,
                                fecell,dofhandler,fe,
                                elmtsystem,matesystem,
                                solutionsystem,
                                equationsystem.m_AMATRIX,equationsystem.m_RHS);


            if(fectrlinfo.CurrentStep%output.getIntervalNum()==0){
                projection.executeProjection(fecell,dofhandler,elmtsystem,matesystem,fe,solutionsystem,fectrlinfo);
                output.saveResults2File(fectrlinfo.CurrentStep,fecell,dofhandler,solutionsystem,projection);
                output.savePVDResults(fectrlinfo.T);
                MessagePrinter::printDashLine(MessageColor::BLUE);
                MessagePrinter::printNormalTxt("Save results to "+output.getOutputFileName(),MessageColor::BLUE);
                MessagePrinter::printDashLine(MessageColor::BLUE);
            }
            if(postprocess.hasPostprocess()){
                if(fectrlinfo.CurrentStep%postprocess.getInterval()==0){
                    if(fectrlinfo.CurrentStep%output.getIntervalNum()!=0){
                        projection.executeProjection(fecell,dofhandler,elmtsystem,matesystem,fe,solutionsystem,fectrlinfo);
                    }
                    postprocess.executePostprocess(fecell,dofhandler,fe,matesystem,projection,solutionsystem);
                    postprocess.savePPSResults2CSVFile(fectrlinfo.T);
                    MessagePrinter::printDashLine(MessageColor::BLUE);
                    MessagePrinter::printNormalTxt("Save postprocess result to "+postprocess.getCSVFileName(),MessageColor::BLUE);
                    MessagePrinter::printDashLine(MessageColor::BLUE);
                }
            }
            MessagePrinter::printStars();
            if(isAdaptive()){
                if(nlsolver.getIterationNum()<=getOptimizeIters() && lastiters<=getOptimizeIters()){
                    fectrlinfo.Dt*=getGrowthFactor();
                    if(fectrlinfo.Dt>getMaxDt()) fectrlinfo.Dt=getMaxDt();
                }
                else if(nlsolver.getIterationNum()<=getOptimizeIters() && lastiters>getOptimizeIters()){
                    // do not change current dt
                    fectrlinfo.Dt*=1.0;
                }
                else if(nlsolver.getIterationNum()>getOptimizeIters()){
                    fectrlinfo.Dt*=getCutbackFactor();
                    if(fectrlinfo.Dt<getMinDt()) fectrlinfo.Dt=getMinDt();
                }
            }
            // update the solution
            solutionsystem.m_Uolder.copyFrom(solutionsystem.m_Uold);
            solutionsystem.m_Uold.copyFrom(solutionsystem.m_Ucurrent);
            solutionsystem.m_V.setToZero();
            // store the previous step's iteration numbers
            lastiters=nlsolver.getIterationNum();
        }
        else{
            //  if the nonlinear solver failed we will try to reduce the dt
            fectrlinfo.Dt*=getCutbackFactor();
            snprintf(buff,68," Transient solver failed, reduce dt to %13.5e",fectrlinfo.Dt);
            str=buff;
            MessagePrinter::printWarningTxt(str);
            if(fectrlinfo.Dt<getMinDt()){
                MessagePrinter::printErrorTxt("The minimum detal t is reached, however, your solver still fails. Please check either your code or your boundary conditions");
                return false;
            }
        }
    }
    if(fectrlinfo.CurrentStep%output.getIntervalNum()!=0){
        projection.executeProjection(fecell,dofhandler,elmtsystem,matesystem,fe,solutionsystem,fectrlinfo);
        output.saveResults2File(fectrlinfo.CurrentStep,fecell,dofhandler,solutionsystem,projection);
        output.savePVDResults(fectrlinfo.T);
        MessagePrinter::printDashLine(MessageColor::BLUE);
        MessagePrinter::printNormalTxt("Save results to "+output.getOutputFileName(),MessageColor::BLUE);
        MessagePrinter::printDashLine(MessageColor::BLUE);
        MessagePrinter::printStars();
    }
    
    return true;
}