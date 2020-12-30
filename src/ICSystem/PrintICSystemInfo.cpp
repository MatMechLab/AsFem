//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.12.30
//+++ Purpose: print out some basic info for ICSystem class
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "ICSystem/ICSystem.h"

void ICSystem::PrintICSystemInfo()const{
    MessagePrinter::PrintNormalTxt("Initial conditions information summary:");
    const int len=68;
    char buff[len];
    string str;
    for(auto it:_ICBlockList){
        snprintf(buff,len," +Initial condition block information:");
        str=buff;
        MessagePrinter::PrintNormalTxt(str);
        //*
        snprintf(buff,len,"   ic block name = [%40s]",it._ICBlockName.c_str());
        str=buff;
        MessagePrinter::PrintNormalTxt(str);
        //*
        snprintf(buff,len,"   ic type name  = %15s",it._ICTypeName.c_str());
        str=buff;
        MessagePrinter::PrintNormalTxt(str);
        //*
        snprintf(buff,len,"   dof name            = %15s, dof index=%2d",it._DofName.c_str(),it._DofID);
        str=buff;
        MessagePrinter::PrintNormalTxt(str);
        //*
        str="   ic parameters       =";
        for(auto value:it._Parameters) str+=to_string(value)+" ";
        MessagePrinter::PrintNormalTxt(str);
        //*
        str="   domain name       =";
        for(auto name:it._DomainNameList){
            str+=name;
        }
        MessagePrinter::PrintNormalTxt(str);
    }
    MessagePrinter::PrintDashLine();
}