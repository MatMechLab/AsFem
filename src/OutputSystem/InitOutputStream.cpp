//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "OutputSystem/OutputSystem.h"

void OutputSystem::InitOutputStream(){
    _InputFileName=_OutputBlock._InputFileName;

    int i=_InputFileName.find_last_of('.');
    if(_OutputBlock._FolderName.size()<1){
        _OutputFilePrefix=_InputFileName.substr(0,i);
        _OutputBlock._OutputFilePrefix=_OutputFilePrefix;
    }
    else{
        // _OutputFilePrefix=_OutputBlock._FolderName+"/"+_InputFileName.substr(0,i);
        _OutputFilePrefix=_InputFileName.substr(0,i);
        /*
        if(std::filesystem::create_directory(_OutputBlock._FolderName)){
            PetscPrintf(PETSC_COMM_WORLD,"***   create folder(%25s) sucessful        !!!   ***\n",_OutputBlock._FolderName.c_str());
        }
        else{
            if(!std::filesystem::exists(_OutputBlock._FolderName)){
                PetscPrintf(PETSC_COMM_WORLD,"***   create folder failed, make sure you have write permission !!!   ***\n");
                Msg_AsFem_Exit();
            }
        }
        */
        _OutputBlock._OutputFilePrefix=_OutputFilePrefix;
    }
    _LogFileName=_OutputFilePrefix+".log";
    _OutputBlock._LogFileName=_LogFileName;
    ofstream _LogFile;
    if(_IsLogOn){
        _LogFile.open(_LogFileName,ios::out);
        if(!_LogFile.is_open()){
            PetscPrintf(PETSC_COMM_WORLD,"*** Error: can\'t create a new log file(=%36s) !!! ***\n",_LogFileName.c_str());
            PetscPrintf(PETSC_COMM_WORLD,"***        please make sure you have write permission           !!!   ***\n");
            Msg_AsFem_Exit();
        }
    }
}