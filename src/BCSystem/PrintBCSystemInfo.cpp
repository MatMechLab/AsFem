//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.12.26
//+++ Purpose: print out some basic info for BCSystem class
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/BCSystem.h"

void BCSystem::PrintBCSystemInfo()const{
    MessagePrinter::PrintNormalTxt("Boundary conditions information summary:");
    const int len=68;
    char buff[len];
    string str;
    for(auto it:_BCBlockList){
        snprintf(buff,len," +Boundary block information:");
        str=buff;
        MessagePrinter::PrintNormalTxt(str);
        //*
        snprintf(buff,len,"   boundary block name = [%40s]",it._BCBlockName.c_str());
        str=buff;
        MessagePrinter::PrintNormalTxt(str);
        //*
        if(it._IsTimeDependent){
            snprintf(buff,len,"   boundary type name  = %15s (Time dependent)",it._BCTypeName.c_str());
            str=buff;
            MessagePrinter::PrintNormalTxt(str);
        }
        else{
            snprintf(buff,len,"   boundary type name  = %15s",it._BCTypeName.c_str());
            str=buff;
            MessagePrinter::PrintNormalTxt(str);
        }
        //*
        snprintf(buff,len,"   dof name            = %15s, dof index=%2d",it._DofName.c_str(),it._DofID);
        str=buff;
        MessagePrinter::PrintNormalTxt(str);
        //*
        snprintf(buff,len,"   boundary value      = %14.6e",it._BCValue);
        str=buff;
        MessagePrinter::PrintNormalTxt(str);
        //*
        str="   boundary name       =";
        for(auto bcname:it._BoundaryNameList){
            str+=bcname+" ";
        }
        MessagePrinter::PrintNormalTxt(str);
    }
    MessagePrinter::PrintDashLine();
}