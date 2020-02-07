//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/walkandthinker/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "OutputSystem/OutputSystem.h"

void OutputSystem::InitOutputStream(){
    int i=_InputFileName.find_last_of('.');
    _OutputFilePrefix=_InputFileName.substr(0,i);
    _LogFileName=_OutputFilePrefix+".log";
    ofstream _LogFile;
    if(_IsLogOn){
        _LogFile.open(_LogFileName,ios::out);
        if(!_LogFile.is_open()){
            cout<<"*** Error: can\'t create a new log file(="<<_LogFileName<<")!!!"<<endl;
            cout<<"***        please make sure you have write permission!***"<<endl;
            Msg_AsFem_Exit();
        }
    }
}