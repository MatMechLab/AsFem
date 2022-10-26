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
//+++ Date   : 2020.06.29
//+++ Purpose: Implement the welcome screen for the initial running
//+++          of AsFem.
//+++          This printer offers the summary information, i.e.
//+++            version of AsFem, PETSc
//+++            news report... and so on.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include <iostream>
#include <string>

#include "petsc.h"
#include "Utils/MessagePrinter.h"

using std::string;

/**
 * print the welcome message in your terminal
 * @param year the integer number of release year
 * @param month the integer number of release month
 * @param day the integer number of release day
 * @param version the double number of release version (2-digital)
 */
void welcome(const int &year,const int &month,const int &day,const double &version){
    PetscInt Major,Minor,SubMinor;
    PetscGetVersionNumber(&Major,&Minor,&SubMinor,NULL);
    char buff[50];
    string str;
    
    MessagePrinter::printStars(MessageColor::BLUE);
    MessagePrinter::printWelcomeTxt("Welcome to use AsFem                                      AAA");
    MessagePrinter::printWelcomeTxt("A Simple Finite Element Method Program                   // \\\\");

    snprintf(buff,50,"Version: %-10.2f  Release @ %4d-%02d-%02d",version,year,month,day);
    str=buff;
    MessagePrinter::printWelcomeTxt(str+"               //   \\\\");

    snprintf(buff,50,"PETSc version: %2d.%2d.%-2d",Major,Minor,SubMinor);
    str=buff;
    MessagePrinter::printWelcomeTxt(str+"                                //     \\\\");
    
    MessagePrinter::printWelcomeTxt("License: GPL-3.0                                      //       \\\\");
    MessagePrinter::printWelcomeTxt("Author: Yang Bai @ M3-Group                          //_________\\\\");
    MessagePrinter::printWelcomeTxt("Contact: yangbai90@outlook.com                      //-----------\\\\");

    MessagePrinter::printWelcomeTxt("QQ Group: 879908352                                //             \\\\");
    MessagePrinter::printWelcomeTxt("Website: https://github.com/M3Group/AsFem         //               \\\\");
    MessagePrinter::printWelcomeTxt("Feel free to use and discuss  .:.                **                 **");
    MessagePrinter::printStars(MessageColor::BLUE);

}
