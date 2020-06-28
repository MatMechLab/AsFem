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
//+++ Date   : 2020.06.28
//+++ Purpose: Implement a general message printer for AsFem
//+++          This class is general to provide almost all the
//+++          message print in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Utils/MessagePrinter.h"

MessagePrinter::MessagePrinter(){
}

void MessagePrinter::PrintTxt(string str){
    string _Head="*** ";
    string _End=" !!! ***";
    if(str.length()<_nWords-_Head.length()-_End.length()){
        PetscPrintf(PETSC_COMM_WORLD,"%s",_Head.c_str());
        PetscPrintf(PETSC_COMM_WORLD,"%s",str.c_str());

        int i1=static_cast<int>(_Head.size());
        int i2=static_cast<int>(_End.size());
        int i3=static_cast<int>(str.size());

        for(int i=0;i<_nWords-i1-i2-i3;i++){
            PetscPrintf(PETSC_COMM_WORLD," ");
        }
        PetscPrintf(PETSC_COMM_WORLD,"%s\n",_End.c_str());
    }
    else{
        string substr1,substr2;
        substr1=str.substr(0,_nWords-_Head.length()-_End.length());
        substr2=str.substr(_nWords-_Head.length()-_End.length());
        PetscPrintf(PETSC_COMM_WORLD,"%s",_Head.c_str());
        PetscPrintf(PETSC_COMM_WORLD,"%s",substr1.c_str());
        PetscPrintf(PETSC_COMM_WORLD,"%s\n",_End.c_str());

        PetscPrintf(PETSC_COMM_WORLD,"%s",_Head.c_str());
        PetscPrintf(PETSC_COMM_WORLD,"%s",substr2.c_str());
        int i1=static_cast<int>(_Head.size());
        int i2=static_cast<int>(_End.size());
        int i3=static_cast<int>(substr2.size());
        for(int i=0;i<_nWords-i1-i2-i3;i++){
            PetscPrintf(PETSC_COMM_WORLD," ");
        }
        PetscPrintf(PETSC_COMM_WORLD,"%s\n",_End.c_str());
    }
}

void MessagePrinter::PrintWelcomeTxt(string str){
    string _Head="*** ";
    string _End =" ***";
    
    PetscPrintf(PETSC_COMM_WORLD,"%s",_Head.c_str());
    PetscPrintf(PETSC_COMM_WORLD,"%s",str.c_str());

    int i1=static_cast<int>(_Head.size());
    int i2=static_cast<int>(_End.size());
    int i3=static_cast<int>(str.size());

    for(int i=0;i<_nWords-i1-i2-i3;i++){
        PetscPrintf(PETSC_COMM_WORLD," ");
    }
    PetscPrintf(PETSC_COMM_WORLD,"%s\n",_End.c_str());
    
}