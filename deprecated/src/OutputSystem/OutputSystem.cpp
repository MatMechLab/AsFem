#include "OutputSystem/OutputSystem.h"

OutputSystem::OutputSystem(){
    _IsInit=false;
    _IsLogOn=false;
    _InputFileName.clear();
    _OutputFilePrefix.clear();
    _LogFileName.clear();
    _VTUFile.clear();_LogFile.clear();
}