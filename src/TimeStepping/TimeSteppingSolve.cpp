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

bool TimeStepping::solve(FECell &t_FECell,
                         DofHandler &t_DofHandler,
                         FE &t_FE,
                         ElmtSystem &t_ElmtSystem,
                         MateSystem &t_MateSystem,
                         FESystem &t_FESystem,
                         BCSystem &t_BCSystem,
                         ICSystem &t_ICSystem,
                         SolutionSystem &t_SolnSystem,
                         EquationSystem &t_EqSystem,
                         ProjectionSystem &t_ProjSystem,
                         FEControlInfo &t_FECtrlInfo,
                         LinearSolver &t_LinearSolver,
                         NonlinearSolver &t_NLSolver,
                         OutputSystem &t_Output,
                         Postprocessor &t_PostProcess){

    char buff[68];//77-12=65
    string str;
    int lastiters=1000;
    t_FECtrlInfo.Dt=m_Data.m_Dt0;
    t_FECtrlInfo.CurrentStep=0;

    t_SolnSystem.m_Ucurrent.setToZero();
    t_ICSystem.applyInitialConditions(t_FECell,t_DofHandler,t_SolnSystem.m_Ucurrent);
    t_SolnSystem.m_Utemp.copyFrom(t_SolnSystem.m_Ucurrent);
    t_SolnSystem.m_Uold.copyFrom(t_SolnSystem.m_Ucurrent);
    t_SolnSystem.m_Uolder.copyFrom(t_SolnSystem.m_Ucurrent);

    // initialize the material
    t_FESystem.formBulkFE(FECalcType::INITMATERIAL,
                          t_FECtrlInfo.T,
                          t_FECtrlInfo.Dt,
                          t_FECtrlInfo.Ctan,
                          t_FECell,
                          t_DofHandler,
                          t_FE,
                          t_ElmtSystem,
                          t_MateSystem,
                          t_SolnSystem,
                          t_EqSystem.m_AMATRIX,
                          t_EqSystem.m_RHS);
    t_ProjSystem.executeProjection(t_FECell,t_DofHandler,t_ElmtSystem,t_MateSystem,t_FE,t_SolnSystem,t_FECtrlInfo);
    MessagePrinter::printDashLine();
    MessagePrinter::printNormalTxt("Material properties have been initialized");
    MessagePrinter::printDashLine();
    
    t_Output.savePVDHead();
    t_Output.savePVDEnd();
    t_Output.saveResults2File(0,t_FECell,t_DofHandler,t_SolnSystem,t_ProjSystem);
    t_Output.savePVDResults(0.0);
    MessagePrinter::printDashLine(MessageColor::BLUE);
    MessagePrinter::printNormalTxt("Save results to "+t_Output.getOutputFileName(),MessageColor::BLUE);
    MessagePrinter::printDashLine(MessageColor::BLUE);
    if(t_PostProcess.hasPostprocess()){
        t_PostProcess.prepareCSVFileHeader();
        t_PostProcess.executePostprocess(t_FECell,t_DofHandler,t_FE,t_MateSystem,t_ProjSystem,t_SolnSystem);
        t_PostProcess.savePPSResults2CSVFile(0.0);
        MessagePrinter::printDashLine(MessageColor::BLUE);
        MessagePrinter::printNormalTxt("Save postprocess result to "+t_PostProcess.getCSVFileName(),MessageColor::BLUE);
        MessagePrinter::printDashLine(MessageColor::BLUE);
    }
    MessagePrinter::printStars();

    for(t_FECtrlInfo.T=0.0;t_FECtrlInfo.T<=m_Data.m_FinalTime;){
        snprintf(buff,68,"Time=%13.5e, step=%8d, dt=%13.5e",t_FECtrlInfo.T+t_FECtrlInfo.Dt,t_FECtrlInfo.CurrentStep+1,t_FECtrlInfo.Dt);
        str=buff;
        MessagePrinter::printNormalTxt(str);
        if(t_NLSolver.solve(t_FECell,
                            t_DofHandler,
                          t_FE,
                          t_ElmtSystem,
                          t_MateSystem,
                          t_FESystem,
                          t_BCSystem,
                          t_SolnSystem,
                          t_EqSystem,
                          t_LinearSolver,
                          t_FECtrlInfo)){
            // if the current nonlinear process success
            t_FECtrlInfo.T+=t_FECtrlInfo.Dt;
            t_FECtrlInfo.CurrentStep+=1;

            // update the material properties
            // update the solution
            t_SolnSystem.m_Utemp.copyFrom(t_SolnSystem.m_Ucurrent);
            t_FESystem.formBulkFE(FECalcType::UPDATEMATERIAL,
                                  t_FECtrlInfo.T,
                                  t_FECtrlInfo.Dt,
                                  t_FECtrlInfo.Ctan,
                                  t_FECell,
                                  t_DofHandler,
                                  t_FE,
                                  t_ElmtSystem,
                                  t_MateSystem,
                                  t_SolnSystem,
                                  t_EqSystem.m_AMATRIX,
                                  t_EqSystem.m_RHS);


            if(t_FECtrlInfo.CurrentStep%t_Output.getIntervalNum()==0){
                t_ProjSystem.executeProjection(t_FECell,t_DofHandler,t_ElmtSystem,t_MateSystem,t_FE,t_SolnSystem,t_FECtrlInfo);
                t_Output.saveResults2File(t_FECtrlInfo.CurrentStep,t_FECell,t_DofHandler,t_SolnSystem,t_ProjSystem);
                t_Output.savePVDResults(t_FECtrlInfo.T);
                MessagePrinter::printDashLine(MessageColor::BLUE);
                MessagePrinter::printNormalTxt("Save results to "+t_Output.getOutputFileName(),MessageColor::BLUE);
                MessagePrinter::printDashLine(MessageColor::BLUE);
            }
            if(t_PostProcess.hasPostprocess()){
                if(t_FECtrlInfo.CurrentStep%t_PostProcess.getInterval()==0){
                    if(t_FECtrlInfo.CurrentStep%t_Output.getIntervalNum()!=0){
                        t_ProjSystem.executeProjection(t_FECell,t_DofHandler,t_ElmtSystem,t_MateSystem,t_FE,t_SolnSystem,t_FECtrlInfo);
                    }
                    t_PostProcess.executePostprocess(t_FECell,t_DofHandler,t_FE,t_MateSystem,t_ProjSystem,t_SolnSystem);
                    t_PostProcess.savePPSResults2CSVFile(t_FECtrlInfo.T);
                    MessagePrinter::printDashLine(MessageColor::BLUE);
                    MessagePrinter::printNormalTxt("Save postprocess result to "+t_PostProcess.getCSVFileName(),MessageColor::BLUE);
                    MessagePrinter::printDashLine(MessageColor::BLUE);
                }
            }
            MessagePrinter::printStars();
            if(isAdaptive()){
                if(t_NLSolver.getIterationNum()<=getOptimizeIters() && lastiters<=getOptimizeIters()){
                    t_FECtrlInfo.Dt*=getGrowthFactor();
                    if(t_FECtrlInfo.Dt>getMaxDt()) t_FECtrlInfo.Dt=getMaxDt();
                }
                else if(t_NLSolver.getIterationNum()<=getOptimizeIters() && lastiters>getOptimizeIters()){
                    // do not change current dt
                    t_FECtrlInfo.Dt*=1.0;
                }
                else if(t_NLSolver.getIterationNum()>getOptimizeIters()){
                    t_FECtrlInfo.Dt*=getCutbackFactor();
                    if(t_FECtrlInfo.Dt<getMinDt()) t_FECtrlInfo.Dt=getMinDt();
                }
            }
            // update the solution
            t_SolnSystem.m_Uolder.copyFrom(t_SolnSystem.m_Uold);
            t_SolnSystem.m_Uold.copyFrom(t_SolnSystem.m_Ucurrent);
            t_SolnSystem.m_V.setToZero();
            // store the previous step's iteration numbers
            lastiters=t_NLSolver.getIterationNum();
        }
        else{
            //  if the nonlinear solver failed we will try to reduce the dt
            t_SolnSystem.m_Ucurrent.copyFrom(t_SolnSystem.m_Uold);
            t_FECtrlInfo.Dt*=getCutbackFactor();
            snprintf(buff,68," Transient solver failed, reduce dt to %13.5e",t_FECtrlInfo.Dt);
            str=buff;
            MessagePrinter::printWarningTxt(str);
            if(t_FECtrlInfo.Dt<getMinDt()){
                MessagePrinter::printErrorTxt("The minimum detal t is reached, however, your solver still fails. Please check either your code or your boundary conditions");
                return false;
            }
        }
    }
    if(t_FECtrlInfo.CurrentStep%t_Output.getIntervalNum()!=0){
        t_ProjSystem.executeProjection(t_FECell,t_DofHandler,t_ElmtSystem,t_MateSystem,t_FE,t_SolnSystem,t_FECtrlInfo);
        t_Output.saveResults2File(t_FECtrlInfo.CurrentStep,t_FECell,t_DofHandler,t_SolnSystem,t_ProjSystem);
        t_Output.savePVDResults(t_FECtrlInfo.T);
        MessagePrinter::printDashLine(MessageColor::BLUE);
        MessagePrinter::printNormalTxt("Save results to "+t_Output.getOutputFileName(),MessageColor::BLUE);
        MessagePrinter::printDashLine(MessageColor::BLUE);
        MessagePrinter::printStars();
    }
    
    return true;
}