//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "InputSystem/InputSystem.h"

InputSystem::InputSystem(){
    _InputFileName.clear();
    _MeshFileName.clear();
    _HasInputFileName=false;
    _IsBuiltInMesh=true;
}
InputSystem::InputSystem(int args,char *argv[]){
    if(args==1){
        // ./asfem
        _InputFileName.clear();
        _MeshFileName.clear();
        _HasInputFileName=false;
        _IsBuiltInMesh=true;
    }
    else if(args==3){
        // ./asfem -i inputfilename.i
        if(string("-i").find(argv[1])!=string::npos){
            _InputFileName=argv[2];
            _HasInputFileName=true;
            if(_InputFileName.compare(_InputFileName.size()-2,2,".i")!=0){
                Msg_Input_InvalidInputFileName();
                Msg_AsFem_Exit();
            }
        }
        else{
            _HasInputFileName=false;
            Msg_Input_InvalidArgs();
            Msg_AsFem_Exit();
        }
    }
    else{
        _HasInputFileName=false;
        Msg_Input_InvalidArgs();
        Msg_AsFem_Exit();
    }
}
//********************************************************
void InputSystem::InitInputSystem(int args,char *argv[]){
    if(args==1){
        // ./asfem
        _InputFileName.clear();
        _MeshFileName.clear();
        _HasInputFileName=false;
        _IsBuiltInMesh=true;
    }
    else if(args==3){
        // ./asfem -i inputfilename.i
        if(string("-i").find(argv[1])!=string::npos){
            _InputFileName=argv[2];
            _HasInputFileName=true;
            if(_InputFileName.compare(_InputFileName.size()-2,2,".i")!=0){
                Msg_Input_InvalidInputFileName();
                Msg_AsFem_Exit();
            }
        }
        else{
            _HasInputFileName=false;
            Msg_Input_InvalidArgs();
            Msg_AsFem_Exit();
        }
    }
    else{
        _HasInputFileName=false;
        Msg_Input_InvalidArgs();
        Msg_AsFem_Exit();
    }
}