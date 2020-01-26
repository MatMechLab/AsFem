#include "OutputSystem/OutputSystem.h"

void OutputSystem::InitOutputStream(){
    int i=_InputFileName.find_last_of('.');
    _OutputFilePrefix=_InputFileName.substr(0,i);
    _LogFileName=_OutputFilePrefix+".log";
    if(_IsLogOn){
        _LogFile.open(_LogFileName,ios::out);
        if(!_LogFile.is_open()){
            cout<<"*** Error: can\'t create a new log file(="<<_LogFileName<<")!!!"<<endl;
            cout<<"***        please make sure you have write permission!***"<<endl;
            Msg_AsFem_Exit();
        }
    }
}