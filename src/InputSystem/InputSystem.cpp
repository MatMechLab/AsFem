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
//+++ Date   : 2020.06.30
//+++ Purpose: Implement the input file reader for AsFem
//+++          So, in this class, the input file for AsFem is only
//+++          accessible via this class
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/InputSystem.h"

InputSystem::InputSystem(){
    _InputFileName.clear();
    _MeshFileName.clear();
    _HasInputFileName=false;
    _IsBuiltInMesh=true;
    _IsReadOnly=false;
}
//**********************************
InputSystem::InputSystem(int args,char *argv[]){
    _InputFileName.clear();
    _MeshFileName.clear();
    _HasInputFileName=false;
    _IsBuiltInMesh=true;
    _IsReadOnly=false;

    if(args==1){
        // ./asfem
        _InputFileName.clear();
        _MeshFileName.clear();
        _HasInputFileName=false;
        _IsBuiltInMesh=true;
        _IsReadOnly=false;
    }
    else if(args==3){
        // ./asfem -i inputfilename.i
        if(string("-i").find(argv[1])!=string::npos){
            _InputFileName=argv[2];
            _HasInputFileName=true;
            if(_InputFileName.compare(_InputFileName.size()-2,2,".i")!=0){
                MessagePrinter::PrintErrorTxt("Invalid input file name! The input file should have the extension .i");
                MessagePrinter::AsFem_Exit();
            }
        }
        else{
            _HasInputFileName=false;
            MessagePrinter::PrintErrorTxt("Invalid input args. The second args should be '-i', after it should be the input file name");
            MessagePrinter::AsFem_Exit();
        }
    }
    else{
        _HasInputFileName=false;
        if(string("-i").find(argv[1])!=string::npos){
            _InputFileName=argv[2];
            _HasInputFileName=true;
            if(_InputFileName.compare(_InputFileName.size()-2,2,".i")!=0){
                MessagePrinter::PrintErrorTxt("Invalid input file name! The input file should have the extension .i");
                MessagePrinter::AsFem_Exit();
            }
        }
        else{
            _HasInputFileName=false;
            MessagePrinter::PrintErrorTxt("Invalid input args. The second args should be '-i', after it should be the input file name");
            MessagePrinter::AsFem_Exit();
        }

        for(int i=4-1;i<args;i++){
            cout<<"argv["<<i<<"]="<<argv[i]<<endl;
            if(string("--read-only").find(argv[i])!=string::npos){
                _IsReadOnly=true;
            }
        }
    }
}
//***************************************************
void InputSystem::InitInputSystem(int args,char *argv[]){
    if(args==1){
        // ./asfem
        _InputFileName.clear();
        _MeshFileName.clear();
        _HasInputFileName=false;
        _IsBuiltInMesh=true;
        _IsReadOnly=true;
    }
    else if(args==3){
        _InputFileName.clear();
        _MeshFileName.clear();
        _HasInputFileName=false;
        _IsBuiltInMesh=true;
        _IsReadOnly=false;
        // ./asfem -i inputfilename.i
        if(string("-i").find(argv[1])!=string::npos){
            _InputFileName=argv[2];
            _HasInputFileName=true;
            if(_InputFileName.compare(_InputFileName.size()-2,2,".i")!=0){
                MessagePrinter::PrintErrorTxt("Invalid input file name! The input file should have the extension .i");
                MessagePrinter::AsFem_Exit();
            }
        }
        else{
            _HasInputFileName=false;
            MessagePrinter::PrintErrorTxt("Invalid input args. The second args should be '-i', after it should be the input file name");
            MessagePrinter::AsFem_Exit();
        }
    }
    else{
        _InputFileName.clear();
        _MeshFileName.clear();
        _HasInputFileName=false;
        _IsBuiltInMesh=true;
        _IsReadOnly=false;
        _HasInputFileName=false;
        if(string("-i").find(argv[1])!=string::npos){
            _InputFileName=argv[2];
            _HasInputFileName=true;
            if(_InputFileName.compare(_InputFileName.size()-2,2,".i")!=0){
                MessagePrinter::PrintErrorTxt("Invalid input file name! The input file should have the extension .i");
                MessagePrinter::AsFem_Exit();
            }
        }
        else{
            _HasInputFileName=false;
            MessagePrinter::PrintErrorTxt("Invalid input args. The second args should be '-i', after it should be the input file name");
            MessagePrinter::AsFem_Exit();
        }
        for(int i=4-1;i<args;i++){
            if(string("--read-only").find(argv[i])!=string::npos){
                _IsReadOnly=true;
            }
        }
    }
}