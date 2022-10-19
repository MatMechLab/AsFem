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
//+++ Date   : 2022.08.07
//+++ Purpose: read the nonlinear solver block 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/InputSystem.h"

bool InputSystem::readNLSolverBlock(nlohmann::json &t_json,NonlinearSolver &t_nlsolver){
    // the json already contains "nlsolver" block
    bool HasType;
    if(t_json.contains("type")){
        if(!t_json.at("type").is_string()){
            MessagePrinter::printErrorTxt("the type name of your nlsolver block is not a valid string");
            return false;
        }
        string solvertypename=t_json.at("type");
        if(solvertypename=="newton"){
            t_nlsolver.m_nlsolverblock.m_nlsolvertypename="newton with line search";
            t_nlsolver.m_nlsolverblock.m_nlsolvertype=NonlinearSolverType::NEWTONLS;
        }
        else if(solvertypename=="newton-ls"){
            t_nlsolver.m_nlsolverblock.m_nlsolvertypename="newton with line search";
            t_nlsolver.m_nlsolverblock.m_nlsolvertype=NonlinearSolverType::NEWTONLS;
        }
        else if(solvertypename=="newton-tr"){
            t_nlsolver.m_nlsolverblock.m_nlsolvertypename="newton with trust region";
            t_nlsolver.m_nlsolverblock.m_nlsolvertype=NonlinearSolverType::NEWTONTR;
        }
        else if(solvertypename=="newton-secant"){
            t_nlsolver.m_nlsolverblock.m_nlsolvertypename="newton secant";
            t_nlsolver.m_nlsolverblock.m_nlsolvertype=NonlinearSolverType::NEWTONSECANT;
        }
        else if(solvertypename=="newton-cg"){
            t_nlsolver.m_nlsolverblock.m_nlsolvertypename="newton cg";
            t_nlsolver.m_nlsolverblock.m_nlsolvertype=NonlinearSolverType::NEWTONCG;
        }
        else if(solvertypename=="newton-gmres"){
            t_nlsolver.m_nlsolverblock.m_nlsolvertypename="newton gmres";
            t_nlsolver.m_nlsolverblock.m_nlsolvertype=NonlinearSolverType::NEWTONGMRES;
        }
        else if(solvertypename=="richardson"){
            t_nlsolver.m_nlsolverblock.m_nlsolvertypename="richardson";
            t_nlsolver.m_nlsolverblock.m_nlsolvertype=NonlinearSolverType::RICHARDSON;
        }
        else if(solvertypename.find("bfgs")!=string::npos){
            t_nlsolver.m_nlsolverblock.m_nlsolvertypename="BFGS";
            t_nlsolver.m_nlsolverblock.m_nlsolvertype=NonlinearSolverType::BFGS;
        }
        else{
            MessagePrinter::printErrorTxt("type="+solvertypename+" is invalid in [nlsolver] block, please check your input file");
            MessagePrinter::exitAsFem();
        }
        HasType=true;
    }
    else{
        MessagePrinter::printWarningTxt("can\'t find 'type'= in your nlsolver block, then the default solver (newton) will be used");
        HasType=true;
    }

    if(t_json.contains("solver")){
        if(!t_json.at("solver").is_string()){
            MessagePrinter::printErrorTxt("the solver name of your nlsolver block is not a valid string");
            return false;
        }
        string solvername=t_json.at("solver");
        if(solvername=="gmres" || solvername=="default"){
            t_nlsolver.m_nlsolverblock.m_linearsolvername="gmres";
        }
        else if(solvername=="fgmres"){
            t_nlsolver.m_nlsolverblock.m_linearsolvername="fgmres";
        }
        else if(solvername=="cg"){
            t_nlsolver.m_nlsolverblock.m_linearsolvername="cg";
        }
        else if(solvername=="bicg"){
            t_nlsolver.m_nlsolverblock.m_linearsolvername="bicg";
        }
        else if(solvername=="richardson"){
            t_nlsolver.m_nlsolverblock.m_linearsolvername="richardson";
        }
        else if(solvername=="mumps"){
            t_nlsolver.m_nlsolverblock.m_linearsolvername="mumps";
        }
        else if(solvername=="superlu"){
            t_nlsolver.m_nlsolverblock.m_linearsolvername="superlu";
        }
    }
    else{
        MessagePrinter::printWarningTxt("can\'t find 'solver' in your nlsolver block, then the default linear solver (gmres) will be used");
        t_nlsolver.m_nlsolverblock.m_linearsolvername="gmres";
    }

    if(t_json.contains("maxiters")){
        if(!t_json.at("maxiters").is_number_integer()){
            MessagePrinter::printErrorTxt("the maxiters in your nlsolver block is not a valid integer,"
                                          "please check your input file");
            return false;
        }
        t_nlsolver.m_nlsolverblock.m_maxiters=t_json.at("maxiters");
    }
    else{
        t_nlsolver.m_nlsolverblock.m_maxiters=25;
    }
    //****************************************
    if(t_json.contains("abs-tolerance")){
        if(!t_json.at("abs-tolerance").is_number_float()){
            MessagePrinter::printErrorTxt("the abs-tolerance in your nlsolver block is not a valid float,"
                                          "please check your input file");
            return false;
        }
        t_nlsolver.m_nlsolverblock.m_abstol_r=t_json.at("abs-tolerance");
    }
    else{
        t_nlsolver.m_nlsolverblock.m_abstol_r=7.0e-7;
    }
    //**********************************************
    if(t_json.contains("rel-tolerance")){
        if(!t_json.at("rel-tolerance").is_number_float()){
            MessagePrinter::printErrorTxt("the rel-tolerance in your nlsolver block is not a valid float,"
                                          "please check your input file");
            return false;
        }
        t_nlsolver.m_nlsolverblock.m_reltol_r=t_json.at("rel-tolerance");
    }
    else{
        t_nlsolver.m_nlsolverblock.m_reltol_r=1.0e-9;
    }
    //**********************************************
    if(t_json.contains("rel-tolerance")){
        if(!t_json.at("rel-tolerance").is_number_float()){
            MessagePrinter::printErrorTxt("the rel-tolerance in your nlsolver block is not a valid float,"
                                          "please check your input file");
            return false;
        }
        t_nlsolver.m_nlsolverblock.m_reltol_r=t_json.at("rel-tolerance");
    }
    else{
        t_nlsolver.m_nlsolverblock.m_reltol_r=1.0e-9;
    }
    //**********************************************
    if(t_json.contains("s-tolerance")){
        if(!t_json.at("s-tolerance").is_number_float()){
            MessagePrinter::printErrorTxt("the s-tolerance in your nlsolver block is not a valid float,"
                                          "please check your input file");
            return false;
        }
        t_nlsolver.m_nlsolverblock.m_s_tol=t_json.at("s-tolerance");
    }
    else{
        t_nlsolver.m_nlsolverblock.m_s_tol=0.0;
    }
    //**********************************************
    if(t_json.contains("preconditioner")){
        if(!t_json.at("preconditioner").is_string()){
            MessagePrinter::printErrorTxt("the preconditioner in your nlsolver block is not a valid string,"
                                          "please check your input file");
            return false;
        }
        t_nlsolver.m_nlsolverblock.m_pctypename=t_json.at("preconditioner");
    }
    else{
        t_nlsolver.m_nlsolverblock.m_pctypename="lu";
    }




    return HasType;
}