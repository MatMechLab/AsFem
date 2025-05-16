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
//+++ Date   : 2025.01.13
//+++ Purpose: read the linear solver block
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/InputSystem.h"

bool InputSystem::readLinearSolverBlock(nlohmann::json &t_json,LinearSolver &t_solver){
    // the json already contains "nlsolver" block
    bool HasType;
    if(t_json.contains("type")){
        if(!t_json.at("type").is_string()){
            MessagePrinter::printErrorTxt("the type name of your linear solver block is not a valid string");
            return false;
        }
        string solvertypename=t_json.at("type");
        if(solvertypename=="gmres"||
           solvertypename=="fgmres"||
           solvertypename=="lgmres"||
           solvertypename=="dgmres"||
           solvertypename=="pgmres"||
           solvertypename=="cg"||
           solvertypename=="chebyshev"||
           solvertypename=="richardson"||
           solvertypename=="bcgs"||
           solvertypename=="bicg"||
           solvertypename=="mumps"||
           solvertypename=="superlu"){
            t_solver.setKSPSolverTypeName(solvertypename);
        }
        else{
            MessagePrinter::printErrorTxt("type="+solvertypename+" is invalid in [linearsolver] block, please check your input file");
            MessagePrinter::exitAsFem();
        }
        HasType=true;
    }
    else{
        MessagePrinter::printWarningTxt("can\'t find 'type'= in your linear solver block, then the default solver (gmres) will be used");
        t_solver.setKSPSolverTypeName("gmres");
        HasType=true;
    }

    if(t_json.contains("preconditioner")){
        if(!t_json.at("preconditioner").is_string()){
            MessagePrinter::printErrorTxt("the preconditioner name of your linear solver block is not a valid string");
            return false;
        }
        string pcname=t_json.at("preconditioner");
        if(pcname=="jacobi"||
           pcname=="bjacobi"||
           pcname=="sor"||
           pcname=="eisenstat"||
           pcname=="icc"||
           pcname=="ilu"||
           pcname=="asm"||
           pcname=="gasm"||
           pcname=="gamg"||
           pcname=="bddc"||
           pcname=="ksp"||
           pcname=="composite"||
           pcname=="lu"||
           pcname=="cholesky"||
           pcname=="none"||
           pcname=="shell"){
            t_solver.setKSPPCTypeName(pcname);
        }
        else{
            MessagePrinter::printErrorTxt("preconditioner="+pcname+" is invalid in [linearsolver] block, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    else{
        MessagePrinter::printWarningTxt("can\'t find 'preconditioner' in your linear solver block, then the default linear solver (bjacobi) will be used");
        t_solver.setKSPPCTypeName("bjacobi");
    }

    if(t_json.contains("maxiters")){
        if(!t_json.at("maxiters").is_number_integer()){
            MessagePrinter::printErrorTxt("the maxiters in your linear solver block is not a valid integer,"
                                          "please check your input file");
            return false;
        }
        t_solver.setMaxIterations(t_json.at("maxiters"));
    }
    else{
        t_solver.setMaxIterations(10000);
    }

    if(t_json.contains("restarts")){
        if(!t_json.at("restarts").is_number_integer()){
            MessagePrinter::printErrorTxt("the restarts in your linear solver block is not a valid integer,"
                                          "please check your input file");
            return false;
        }
        t_solver.setGMRESRestartNumber(t_json.at("restarts"));
    }
    else{
        t_solver.setGMRESRestartNumber(2000);
    }
    //****************************************
    if(t_json.contains("tolerance")){
        if(!t_json.at("tolerance").is_number_float()){
            MessagePrinter::printErrorTxt("the tolerance in your linear solver block is not a valid float,"
                                          "please check your input file");
            return false;
        }
        t_solver.setSolverTolerance(t_json.at("tolerance"));
    }
    else{
        t_solver.setSolverTolerance(1.0e-16);
    }


    return HasType;
}