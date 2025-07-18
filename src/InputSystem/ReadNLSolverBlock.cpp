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
        if(solvertypename=="asfem"){
            t_nlsolver.m_NlSolverBlock.m_NlSolverTypeName="newton-raphson from asfem";
            t_nlsolver.m_NlSolverBlock.m_NlSolverType=NonlinearSolverType::ASFEMNR;
        }
        else if(solvertypename=="newton"){
            t_nlsolver.m_NlSolverBlock.m_NlSolverTypeName="newton with line search";
            t_nlsolver.m_NlSolverBlock.m_NlSolverType=NonlinearSolverType::NEWTONLS;
        }
        else if(solvertypename=="newton-ls"){
            t_nlsolver.m_NlSolverBlock.m_NlSolverTypeName="newton with line search";
            t_nlsolver.m_NlSolverBlock.m_NlSolverType=NonlinearSolverType::NEWTONLS;
        }
        else if(solvertypename=="newton-al"){
            t_nlsolver.m_NlSolverBlock.m_NlSolverTypeName="newton with arc length method";
            t_nlsolver.m_NlSolverBlock.m_NlSolverType=NonlinearSolverType::NEWTONAL;
        }
        else if(solvertypename=="newton-tr"){
            t_nlsolver.m_NlSolverBlock.m_NlSolverTypeName="newton with trust region";
            t_nlsolver.m_NlSolverBlock.m_NlSolverType=NonlinearSolverType::NEWTONTR;
        }
        else if(solvertypename=="newton-secant"){
            t_nlsolver.m_NlSolverBlock.m_NlSolverTypeName="newton secant";
            t_nlsolver.m_NlSolverBlock.m_NlSolverType=NonlinearSolverType::NEWTONSECANT;
        }
        else if(solvertypename=="newton-cg"){
            t_nlsolver.m_NlSolverBlock.m_NlSolverTypeName="newton cg";
            t_nlsolver.m_NlSolverBlock.m_NlSolverType=NonlinearSolverType::NEWTONCG;
        }
        else if(solvertypename=="newton-gmres"){
            t_nlsolver.m_NlSolverBlock.m_NlSolverTypeName="newton gmres";
            t_nlsolver.m_NlSolverBlock.m_NlSolverType=NonlinearSolverType::NEWTONGMRES;
        }
        else if(solvertypename=="richardson"){
            t_nlsolver.m_NlSolverBlock.m_NlSolverTypeName="richardson";
            t_nlsolver.m_NlSolverBlock.m_NlSolverType=NonlinearSolverType::RICHARDSON;
        }
        else if(solvertypename=="ncg"){
            t_nlsolver.m_NlSolverBlock.m_NlSolverTypeName="ncg";
            t_nlsolver.m_NlSolverBlock.m_NlSolverType=NonlinearSolverType::NCG;
        }
        else if(solvertypename.find("bfgs")!=string::npos){
            t_nlsolver.m_NlSolverBlock.m_NlSolverTypeName="BFGS";
            t_nlsolver.m_NlSolverBlock.m_NlSolverType=NonlinearSolverType::BFGS;
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

    if(t_json.contains("maxiters")){
        if(!t_json.at("maxiters").is_number_integer()){
            MessagePrinter::printErrorTxt("the maxiters in your nlsolver block is not a valid integer,"
                                          "please check your input file");
            return false;
        }
        t_nlsolver.m_NlSolverBlock.m_MaxIters=t_json.at("maxiters");
    }
    else{
        t_nlsolver.m_NlSolverBlock.m_MaxIters=25;
    }
    //****************************************
    if(t_json.contains("abs-tolerance")){
        if(!t_json.at("abs-tolerance").is_number_float()){
            MessagePrinter::printErrorTxt("the abs-tolerance in your nlsolver block is not a valid float,"
                                          "please check your input file");
            return false;
        }
        t_nlsolver.m_NlSolverBlock.m_AbsTolR=t_json.at("abs-tolerance");
    }
    else{
        t_nlsolver.m_NlSolverBlock.m_AbsTolR=7.0e-7;
    }
    //**********************************************
    if(t_json.contains("rel-tolerance")){
        if(!t_json.at("rel-tolerance").is_number_float()){
            MessagePrinter::printErrorTxt("the rel-tolerance in your nlsolver block is not a valid float,"
                                          "please check your input file");
            return false;
        }
        t_nlsolver.m_NlSolverBlock.m_RelTolR=t_json.at("rel-tolerance");
    }
    else{
        t_nlsolver.m_NlSolverBlock.m_RelTolR=1.0e-9;
    }
    //**********************************************
    if(t_json.contains("s-tolerance")){
        if(!t_json.at("s-tolerance").is_number_float()){
            MessagePrinter::printErrorTxt("the s-tolerance in your nlsolver block is not a valid float,"
                                          "please check your input file");
            return false;
        }
        t_nlsolver.m_NlSolverBlock.m_STol=t_json.at("s-tolerance");
    }
    else{
        t_nlsolver.m_NlSolverBlock.m_STol=0.0;
    }

    return HasType;
}