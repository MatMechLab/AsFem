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
//+++ Date   : 2020.12.29
//+++ Purpose: Define time stepping system for transient analysis
//+++          in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "TimeStepping/TimeStepping.h"

TimeStepping::TimeStepping(){
    _Dt=1.0e-5;
    _FinalT=1.0e-3;
    _Adaptive=false;
    _TotalSteps=-1;
    _TimeSteppingType=TimeSteppingType::BACKWARDEULER;
    _TimeSteppingTypeName="backward-euler";
    _GrowthFactor=1.1;
    _CutBackFactor=0.85;
    _OptIters=3;
    _DtMax=1.0e2;
    _DtMin=1.0e-12;
    _IterHist[0]=0;
    _IterHist[1]=1;
}

//****************************************************
void TimeStepping::SetOpitonsFromTimeSteppingBlock(TimeSteppingBlock &timeSteppingBlock){
    _Dt=timeSteppingBlock._Dt;
    _DtMin=timeSteppingBlock._DtMin;
    _DtMax=timeSteppingBlock._DtMax;
    _FinalT=timeSteppingBlock._FinalT;
    _TotalSteps=static_cast<long int>(_FinalT/_Dt);
    _TimeSteppingType=timeSteppingBlock._TimeSteppingType;
    _TimeSteppingTypeName=timeSteppingBlock._TimeSteppingTypeName;
    _Adaptive=timeSteppingBlock._Adaptive;
    _GrowthFactor=timeSteppingBlock._GrowthFactor;
    _CutBackFactor=timeSteppingBlock._CutBackFactor;
    _OptIters=timeSteppingBlock._OptIters;
}
//*******************************************************
void TimeStepping::PrintTimeSteppingInfo()const{
    MessagePrinter::PrintNormalTxt("Time stepping system information summary:");
    MessagePrinter::PrintNormalTxt("  stepping method ="+_TimeSteppingTypeName);
    if(_Adaptive){
        MessagePrinter::PrintNormalTxt("  adaptive is enabled, optimal iters ="+to_string(_OptIters));
        MessagePrinter::PrintNormalTxt("  adaptive growth factor="+to_string(_GrowthFactor)+", cut factor="+to_string(_CutBackFactor));
    }
    else{
        MessagePrinter::PrintNormalTxt("  adaptive is disabled");
    }
    char buff[20];
    snprintf(buff,20,"%14.5e",_Dt);
    string str;
    str="  init delta T="+string(buff);
    snprintf(buff,20,"%14.5e",_FinalT);
    str+=", final T="+string(buff);
    MessagePrinter::PrintNormalTxt(str);

    snprintf(buff,20,"%14.5e",_DtMax);
    str="  max delta T="+string(buff);
    snprintf(buff,20,"%14.5e",_DtMin);
    str+=", min delta T="+string(buff);
    MessagePrinter::PrintNormalTxt(str);
    MessagePrinter::PrintDashLine();
}
//****************************************
