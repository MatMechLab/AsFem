//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.06.29
//+++ Purpose: Implement the welcome screen for the initial running
//+++          of AsFem.
//+++          This printer offers the summary information, i.e.
//+++            version of AsFem, PETSc
//+++            news report... and so on.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include <iostream>
#include <string>
#include <cstdio>

#include "petsc.h"
#include "Utils/MessagePrinter.h"

using namespace std;
void Welcome(const PetscInt &year,const PetscInt &month,const PetscInt &day,const PetscReal &version){
    PetscInt Major,Minor,SubMinor;
    PetscGetVersionNumber(&Major,&Minor,&SubMinor,NULL);
    char buff[50];
    string str;
    
    MessagePrinter::PrintStars(MessageColor::CYAN);
    MessagePrinter::PrintWelcomeTxt("Welcome to use AsFem");
    MessagePrinter::PrintWelcomeTxt("A Simple Finite Element Method Program");

    snprintf(buff,50,"Version: %-10.2f  Release @ %4d-%02d-%02d",version,year,month,day);
    str=buff;
    MessagePrinter::PrintWelcomeTxt(str);

    snprintf(buff,50,"PETSc version: %2d.%2d.%-2d",Major,Minor,SubMinor);
    str=buff;
    MessagePrinter::PrintWelcomeTxt(str);
    
    MessagePrinter::PrintWelcomeTxt("License: GPL-3.0");
    MessagePrinter::PrintWelcomeTxt("Author: Yang Bai");
    MessagePrinter::PrintWelcomeTxt("Contact: walkandthinker@gmail.com");

    MessagePrinter::PrintWelcomeTxt("QQ Group: 879908352");
    MessagePrinter::PrintWelcomeTxt("Website: https://github.com/yangbai90/AsFem");
    MessagePrinter::PrintWelcomeTxt("Feel free to use and discuss  .:.");
    MessagePrinter::PrintStars(MessageColor::CYAN);
}
