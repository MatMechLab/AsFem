//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************


#ifndef ASFEM_OUTPUTBLOCK_H
#define ASFEM_OUTPUTBLOCK_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include "petsc.h"

#include "OutputFileType.h"

using namespace std;

//******************************************************************
//*** you output block in your input file should look like:
//*** [output]
//***   type=vtu(vtk,csv,txt...)
//***   folder=folder_name [default value should be empty]
//*** [end]

class OutputBlock{
public:
    OutputFileType _OutputFileType=OutputFileType::VTU;
    string _FolderName="";// default value is empty
    string _InputFileName="";
    string _OutputFilePrefix="";
    string _LogFileName="";
    bool _IsPrint=false;

    void Reset(){
        _OutputFileType=OutputFileType::VTU;
        _FolderName.clear();
        _InputFileName.clear();
        _OutputFilePrefix.clear();
        _LogFileName.clear();
        _IsPrint=false;
    }

    void PrintOutputBlock()const{
        PetscPrintf(PETSC_COMM_WORLD,"***-------------------------------------------------------------------***\n");
        
        PetscPrintf(PETSC_COMM_WORLD,"*** Output block information:                                         ***\n");
        
        if(_OutputFileType==OutputFileType::VTU){
            PetscPrintf(PETSC_COMM_WORLD,"***   output file type = VTU                                          ***\n");
        }
        else if(_OutputFileType==OutputFileType::VTK){
            PetscPrintf(PETSC_COMM_WORLD,"***   output file type = VTK                                          ***\n");
        }
        else if(_OutputFileType==OutputFileType::CSV){
            PetscPrintf(PETSC_COMM_WORLD,"***   output file type = CSV                                          ***\n");
        }
        if(_FolderName.size()<1){
            PetscPrintf(PETSC_COMM_WORLD,"***   output folder    = emtpy                                        ***\n");
        }
        else{
            PetscPrintf(PETSC_COMM_WORLD,"***   output folder    = %-42s   ***\n",_FolderName.c_str());
        }

        PetscPrintf(PETSC_COMM_WORLD,"***   output fileprefix= %-42s   ***\n",_OutputFilePrefix.c_str());
    }
};

#endif