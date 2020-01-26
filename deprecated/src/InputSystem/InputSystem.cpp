#include "InputSystem/InputSystem.h"

InputSystem::InputSystem()
{
    _InputFileName="";
    _HasInputFileName=false;
}
InputSystem::InputSystem(int args,char *argv[])
{
    if(args==1)
    {
        // ./asfem
        _InputFileName="";
        _HasInputFileName=false;
    }
    else if(args==3)
    {
        // ./asfem -i inputfilename.i
        if(string("-i").find(argv[1])!=string::npos)
        {
            _InputFileName=argv[2];
            _HasInputFileName=true;
        }
        else
        {
            _HasInputFileName=false;
            Msg_Input_InvalidArgs();
            Msg_AsFem_Exit();
        }
    }
    else
    {
        _HasInputFileName=false;
        Msg_Input_InvalidArgs();
        Msg_AsFem_Exit();
    }
}