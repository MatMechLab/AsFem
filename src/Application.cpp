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
//+++ Date   : 2024.01.31
//+++ Purpose: the application level management for AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Application.h"

Application::Application(){}

PetscErrorCode Application::init(int args,char *argv[]){
    PetscCall(PetscInitialize(&args,&argv,NULL,NULL));
    /**
     * setup the necessary petsc options
    */
    PetscCall(PetscOptionsCreate(&m_Options));
    PetscCall(PetscOptionsInsertString(m_Options,"-options_left no"));
    PetscCall(PetscOptionsPush(m_Options));
    PetscCall(PetscOptionsLeft(m_Options));
    return 0;
}
//*****************************************
PetscErrorCode Application::finalize(){
    // uncomment following two lines will raise the PETSc warning !
    // PetscCall(PetscOptionsPop());
    // PetscCall(PetscOptionsDestroy(&m_Options));
    PetscCall(PetscFinalize());
    return 0;
}
//***********************************************
void Application::printAppInfo(const int &year,const int &month,const int &day,const double &version)const{
    PetscInt Major,Minor,SubMinor;
    PetscGetVersionNumber(&Major,&Minor,&SubMinor,NULL);
    char buff[50];
    string str;
    
    MessagePrinter::printStars(MessageColor::BLUE);
    MessagePrinter::printWelcomeTxt("Welcome to use AsFem                                      AAA");
    
    MessagePrinter::printSingleTxt("*** ",MessageColor::BLUE);
    MessagePrinter::printSingleTxt("A",MessageColor::BLUE);
    MessagePrinter::printSingleTxt("dvanced ",MessageColor::WHITE);
    //
    MessagePrinter::printSingleTxt("S",MessageColor::BLUE);
    MessagePrinter::printSingleTxt("imulation ",MessageColor::WHITE);
    //
    MessagePrinter::printSingleTxt("kit based on ",MessageColor::WHITE);
    //
    MessagePrinter::printSingleTxt("F",MessageColor::BLUE);
    MessagePrinter::printSingleTxt("inite ",MessageColor::WHITE);
    //
    MessagePrinter::printSingleTxt("E",MessageColor::BLUE);
    MessagePrinter::printSingleTxt("lement ",MessageColor::WHITE);
    //
    MessagePrinter::printSingleTxt("M",MessageColor::BLUE);
    MessagePrinter::printSingleTxt("ethod ",MessageColor::WHITE);
    MessagePrinter::printSingleTxt("  // \\\\",MessageColor::BLUE);
    MessagePrinter::printSingleTxt("                ***",MessageColor::BLUE);
    MessagePrinter::printNewLine();

    snprintf(buff,50,"Version: %-10.2f  Release @ %4d-%02d-%02d",version,year,month,day);
    str=buff;
    MessagePrinter::printWelcomeTxt(str+"               //   \\\\");

    snprintf(buff,50,"PETSc version: %2d.%2d.%-2d",static_cast<int>(Major),static_cast<int>(Minor),static_cast<int>(SubMinor));
    str=buff;
    MessagePrinter::printWelcomeTxt(str+"                                //     \\\\");
    
    MessagePrinter::printWelcomeTxt("License: GPL-3.0                                      //       \\\\");
    MessagePrinter::printWelcomeTxt("Author: Yang Bai @ MM-Lab                            //_________\\\\");
    MessagePrinter::printWelcomeTxt("Contact: yangbai90@outlook.com                      //-----------\\\\");

    MessagePrinter::printWelcomeTxt("QQ Group: 879908352                                //             \\\\");
    MessagePrinter::printWelcomeTxt("Website: https://github.com/MatMechLab/AsFem      //               \\\\");
    MessagePrinter::printWelcomeTxt("Feel free to use and discuss  .:.                **                 **");
    MessagePrinter::printStars(MessageColor::BLUE);
}