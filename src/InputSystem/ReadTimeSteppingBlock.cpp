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

bool InputSystem::readTimeSteppingBlock(nlohmann::json &t_json,TimeStepping &t_timestepping){
    // the json already contains "nlsolver" block
    bool HasType;
    bool Adaptive=false;
    if(t_json.contains("type")){
        if(!t_json.at("type").is_string()){
            MessagePrinter::printErrorTxt("the type name of your timestepping block is not a valid string");
            return false;
        }
        string solvertypename=t_json.at("type");
        if(solvertypename=="backward-euler" || 
           solvertypename=="be" ||
           solvertypename=="BE"){
            t_timestepping.setTimeSteppingMethod(TimeSteppingType::BACKWARDEULER);
        }
        else if(solvertypename=="cranck-nicolson" ||
                solvertypename=="cn" ||
                solvertypename=="CN"){
            t_timestepping.setTimeSteppingMethod(TimeSteppingType::CRANCKNICOLSON);
        }
        else if(solvertypename=="static"||
                solvertypename=="STATIC"){
            t_timestepping.setTimeSteppingMethod(TimeSteppingType::STATIC);
        }
        else if(solvertypename=="bdf2"||
                solvertypename=="BDF2"){
            t_timestepping.setTimeSteppingMethod(TimeSteppingType::BDF2);
        }
        else{
            MessagePrinter::printErrorTxt("type="+solvertypename+" is invalid in [timestepping] block, please check your input file");
            MessagePrinter::exitAsFem();
        }
        HasType=true;
    }
    else{
        MessagePrinter::printWarningTxt("can\'t find 'type' in your timestepping block, then the default method will be used");
        t_timestepping.setTimeSteppingMethod(TimeSteppingType::BACKWARDEULER);
        HasType=true;
    }

    if(t_json.contains("adaptive")){
        if(!t_json.at("adaptive").is_boolean()){
            MessagePrinter::printErrorTxt("the adaptive option of your timestepping block is not a valid boolean");
            return false;
        }
        Adaptive=t_json.at("adaptive");
        t_timestepping.setAdaptiveFlag(Adaptive);
    }
    else{
        Adaptive=false;
        t_timestepping.setAdaptiveFlag(false);
    }

    if(t_json.contains("optimize-iters")){
        if(!t_json.at("optimize-iters").is_number()){
            MessagePrinter::printErrorTxt("the optimize iterations of your timestepping block is not a valid number");
            return false;
        }
        int iters=t_json.at("optimize-iters");
        if(iters<1){
            MessagePrinter::printErrorTxt("the optimize iterations of your timestepping block is not a valid number");
            return false;
        }
        t_timestepping.setOptimizeIters(iters);
    }
    else{
        t_timestepping.setOptimizeIters(3);
        if(Adaptive){
            MessagePrinter::printWarningTxt("can\'t find optimize iterations in your timestepping block, then the default value will be used");
        }
    }

    if(t_json.contains("dt0")){
        if(!t_json.at("dt0").is_number()){
            MessagePrinter::printErrorTxt("the dt0 in your timestepping block is not a valid float,"
                                          "please check your input file");
            return false;
        }
        double dt0=t_json.at("dt0");
        if(dt0<1.0e-13){
            MessagePrinter::printErrorTxt("the dt0 in your timestepping block is too small or negative,"
                                          "please check your input file");
            return false;
        }
        t_timestepping.setInitialDt(dt0);
    }
    else{
        MessagePrinter::printErrorTxt("can\'t find 'dt0' in your timestepping block, please check your input file");
        return false;
    }

    if(t_json.contains("dtmax")){
        if(!t_json.at("dtmax").is_number()){
            MessagePrinter::printErrorTxt("the dtmax in your timestepping block is not a valid float,"
                                          "please check your input file");
            return false;
        }
        double dtmax=t_json.at("dtmax");
        if(dtmax<1.0e-13){
            MessagePrinter::printErrorTxt("the dtmax in your timestepping block is too small or negative,"
                                          "please check your input file");
            return false;
        }
        t_timestepping.setMaxDt(dtmax);
    }
    else{
        if(Adaptive){
            MessagePrinter::printWarningTxt("can\'t find 'dtmax' in your timestepping block, then the default dtmax value will be used");
        }
        t_timestepping.setMaxDt(1.0e-2);
    }

    if(t_json.contains("dtmin")){
        if(!t_json.at("dtmin").is_number()){
            MessagePrinter::printErrorTxt("the dtmin in your timestepping block is not a valid float,"
                                          "please check your input file");
            return false;
        }
        double dtmin=t_json.at("dtmin");
        if(dtmin<1.0e-13){
            MessagePrinter::printErrorTxt("the dtmin in your timestepping block is too small or negative,"
                                          "please check your input file");
            return false;
        }
        t_timestepping.setMinDt(dtmin);
    }
    else{
        if(Adaptive){
            MessagePrinter::printWarningTxt("can\'t find 'dtmin' in your timestepping block, then the default dtmin value will be used");
        }
        t_timestepping.setMinDt(1.0e-13);
    }


    if(t_json.contains("end-time")){
        if(!t_json.at("end-time").is_number()){
            MessagePrinter::printErrorTxt("the end-time in your timestepping block is not a valid float,"
                                          "please check your input file");
            return false;
        }
        double endtime=t_json.at("end-time");
        if(endtime<1.0e-13){
            MessagePrinter::printErrorTxt("the end-time in your timestepping block is too small or negative,"
                                          "please check your input file");
            return false;
        }
        t_timestepping.setFinalTime(endtime);
    }
    else{
        MessagePrinter::printErrorTxt("can\'t find 'end-time' in your timestepping block, please check your input file");
        return false;
    }

    if(t_json.contains("growth-factor")){
        if(!t_json.at("growth-factor").is_number()){
            MessagePrinter::printErrorTxt("the growth-factor in your timestepping block is not a valid float,"
                                          "please check your input file");
            return false;
        }
        double growthfactor=t_json.at("growth-factor");
        if(growthfactor<1.0e-13){
            MessagePrinter::printErrorTxt("the growth-factor in your timestepping block is too small or negative,"
                                          "please check your input file");
            return false;
        }
        t_timestepping.setGrowthFactor(growthfactor);
    }
    else{
        if(Adaptive){
            MessagePrinter::printErrorTxt("can\'t find 'growth-factor' in your timestepping block, please check your input file");
        }
        return false;
    }

    if(t_json.contains("cutback-factor")){
        if(!t_json.at("cutback-factor").is_number()){
            MessagePrinter::printErrorTxt("the cutback-factor in your timestepping block is not a valid float,"
                                          "please check your input file");
            return false;
        }
        double cutbackfactor=t_json.at("cutback-factor");
        if(cutbackfactor<1.0e-13){
            MessagePrinter::printErrorTxt("the cutback-factor in your timestepping block is too small or negative,"
                                          "please check your input file");
            return false;
        }
        t_timestepping.setCutbackFactor(cutbackfactor);
    }
    else{
        if(Adaptive){
            MessagePrinter::printErrorTxt("can\'t find 'cutback-factor' in your timestepping block, please check your input file");
        }
        return false;
    }

    return HasType;
}