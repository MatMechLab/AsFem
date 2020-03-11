//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "InputSystem/InputSystem.h"

bool InputSystem::ReadTimeSteppingBlock(ifstream &in,string str,int &linenum,TimeSteppingBlock &timesteppingblock){
    bool HasEndTime=false;
    bool HasType=false;
    bool HasDt=false;
    vector<string> namelist;
    vector<double> numbers;
    // now str already contains [dofs]
    getline(in,str);linenum+=1;
    str=StrToLower(str);
    namelist.clear();

    timesteppingblock.Reset();

    HasType=false;
    HasEndTime=false;
    
    while(str.find("[end]")==string::npos&&
          str.find("[END]")==string::npos){
        if(IsCommentLine(str)||str.length()<1){
            getline(in,str);linenum+=1;
            str=StrToLower(str);
            continue;
        }
        str=RemoveStrSpace(str);
        if(str.compare(0,5,"type=")==0||str.compare(0,5,"type=")==0){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            if(substr.find("backwardeuler")!=string::npos||
               substr.find("BACKWARDEULER")!=string::npos){
                   timesteppingblock._TimeSteppingMethodName="backward-euler";
                   timesteppingblock._TimeSteppingMethod=TimeSteppingType::BackWardEuler;
                   HasType=true;
            }
            else if(substr.find("be")!=string::npos||
                    substr.find("BE")!=string::npos){
                timesteppingblock._TimeSteppingMethodName="backward-euler";
                timesteppingblock._TimeSteppingMethod=TimeSteppingType::BackWardEuler;
                HasType=true;
            }
            else if(substr.find("cranknicolson")!=string::npos||
                    substr.find("CRANKNICOLSON")!=string::npos){
                timesteppingblock._TimeSteppingMethodName="cranck-nicolson";
                timesteppingblock._TimeSteppingMethod=TimeSteppingType::CrankNicolson;
                HasType=true;
            }
            else if(substr.find("cn")!=string::npos||
                    substr.find("CN")!=string::npos){
                timesteppingblock._TimeSteppingMethodName="cranck-nicolson";
                timesteppingblock._TimeSteppingMethod=TimeSteppingType::CrankNicolson;
                HasType=true;
            }
            else if(substr.find("alpha")!=string::npos||
                    substr.find("ALPHA")!=string::npos){
                timesteppingblock._TimeSteppingMethodName="alpha-method";
                timesteppingblock._TimeSteppingMethod=TimeSteppingType::AlphaMethod;
                HasType=true;
            }
            else if(substr.find("theta")!=string::npos||
                    substr.find("THETA")!=string::npos){
                timesteppingblock._TimeSteppingMethodName="theta-method";
                timesteppingblock._TimeSteppingMethod=TimeSteppingType::ThetaMethod;
                HasType=true;
            }
            else{
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: unsupported type= in [timestepping] block            !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        type=be[,cn,theta,alpha] is expected                 !!!   ***\n");
                Msg_AsFem_Exit();
            }
        }
        else if(str.find("endtime=")!=string::npos||
                str.find("ENDTIME=")!=string::npos){
            if(!HasType){
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: no type= found in [timestepping] block               !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        endtime= must be given after type= in [timestepping] !!!   ***\n");
                Msg_AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            numbers=SplitStrNum(substr);
            if(numbers.size()<1){
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: no endtime found in [timestepping] block             !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        endtime=real value is expected                       !!!   ***\n");
                Msg_AsFem_Exit();
            }
            else{
                if(numbers[0]<1.0e-10||numbers[0]>1.0e20){
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: invalid endtime found in [timestepping] block        !!!   ***\n");
                    PetscPrintf(PETSC_COMM_WORLD,"***        endtime=positive real value is expected              !!!   ***\n");
                    Msg_AsFem_Exit();
                }
                timesteppingblock._FinalTime=numbers[0];
                HasEndTime=true;
            }
        }
        else if(str.find("dt=")!=string::npos||
                str.find("DT=")!=string::npos){
            if(!HasType){
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: no type= found in [timestepping] block               !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        dt= must be given after type= in [timestepping]      !!!   ***\n");
                Msg_AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            numbers=SplitStrNum(substr);
            if(numbers.size()<1){
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: no dt= found in [timestepping] block                 !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        dt=real value is expected                            !!!   ***\n");
                Msg_AsFem_Exit();
            }
            else{
                if(numbers[0]<1.0e-10){
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: invalid dt= found in [timestepping] block            !!!   ***\n");
                    PetscPrintf(PETSC_COMM_WORLD,"***        dt=positive real value is expected                   !!!   ***\n");
                    Msg_AsFem_Exit();
                }
                timesteppingblock._dt0=numbers[0];
                HasDt=true;
            }
        }
        else if(str.find("adaptive=")!=string::npos||
                str.find("ADAPTIVE=")!=string::npos){
            if(!HasType){
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: no type= found in [timestepping] block               !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        adaptive= must be given after type= in [timestepping]!!!   ***\n");
                Msg_AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            if(substr.find("true")!=string::npos||
               substr.find("TRUE")!=string::npos){
                timesteppingblock._IsAdaptive=true;
            }
            else if(substr.find("false")!=string::npos||
                    substr.find("FALSE")!=string::npos){
                timesteppingblock._IsAdaptive=false;
            }
            else{
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: invalid adaptive= found in [timestepping] block      !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        adaptive=true[false] is expected                     !!!   ***\n");
                Msg_AsFem_Exit();
            }
        }
        else if(str.find("dtmax=")!=string::npos||
                str.find("DTMAX=")!=string::npos){
            if(!HasType){
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: no type= found in [timestepping] block               !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        dtmax= must be given after type= in [timestepping]   !!!   ***\n");
                Msg_AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            numbers=SplitStrNum(substr);
            if(numbers.size()<1){
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: no dtmax= found in [timestepping] block              !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        dtmax=real value is expected                         !!!   ***\n");
                Msg_AsFem_Exit();
            }
            else{
                if(numbers[0]<1.0e-12){
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: invalid dtmax= found in [timestepping] block         !!!   ***\n");
                    PetscPrintf(PETSC_COMM_WORLD,"***        dtmax=positive real value is expected                !!!   ***\n");
                    Msg_AsFem_Exit();
                }
                timesteppingblock._dtmax=numbers[0];
            }
        }
        else if(str.find("dtmin=")!=string::npos||
                str.find("DTMIN=")!=string::npos){
            if(!HasType){
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: no type= found in [timestepping] block               !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        dtmin= must be given after type= in [timestepping]   !!!   ***\n");
                Msg_AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            numbers=SplitStrNum(substr);
            if(numbers.size()<1){
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: no dtmin= found in [timestepping] block              !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        dtmin=real value is expected                         !!!   ***\n");
                Msg_AsFem_Exit();
            }
            else{
                if(numbers[0]<1.0e-12){
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: invalid dtmin= found in [timestepping] block         !!!   ***\n");
                    PetscPrintf(PETSC_COMM_WORLD,"***        dtmin=positive real value is expected                !!!   ***\n");
                    Msg_AsFem_Exit();
                }
                timesteppingblock._dtmin=numbers[0];
            }
        }
        else if(str.find("interval=")!=string::npos||
                str.find("INTERVAL=")!=string::npos){
            if(!HasType){
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: no type= found in [timestepping] block               !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        interval= must be given after type= in [timestepping]!!!   ***\n");
                Msg_AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            numbers=SplitStrNum(substr);
            if(numbers.size()<1){
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: no dt= found in [timestepping] block                 !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        dt=real value is expected                            !!!   ***\n");
                Msg_AsFem_Exit();
            }
            else{
                if(numbers[0]<1){
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: invalid interval= found in [timestepping] block      !!!   ***\n");
                    PetscPrintf(PETSC_COMM_WORLD,"***        interval=positive integer value is expected          !!!   ***\n");
                    Msg_AsFem_Exit();
                }
                timesteppingblock._interval=int(numbers[0]);
            }
        }
        else if(str.find("opts=")!=string::npos||
                str.find("OPTS=")!=string::npos){
            if(!HasType){
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: no type= found in [timestepping] block               !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        opts= must be given after type= in [timestepping]    !!!   ***\n");
                Msg_AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            numbers=SplitStrNum(substr);
            if(numbers.size()<1){
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: no opts= found in [timestepping] block               !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        opts=integer value is expected                       !!!   ***\n");
                Msg_AsFem_Exit();
            }
            else{
                if(numbers[0]<1){
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: invalid opts= found in [timestepping] block          !!!   ***\n");
                    PetscPrintf(PETSC_COMM_WORLD,"***        opts=positive integer value is expected              !!!   ***\n");
                    Msg_AsFem_Exit();
                }
                timesteppingblock._nOpts=int(numbers[0]);
            }
        }
        else if(str.find("cutfactor=")!=string::npos||
                str.find("CUTFACTOR=")!=string::npos){
            if(!HasType){
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: no type= found in [timestepping] block               !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        cutfactor= must be given after type in [timestepping]!!!   ***\n");
                Msg_AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            numbers=SplitStrNum(substr);
            if(numbers.size()<1){
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: no cutfactor= found in [timestepping] block          !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        cutfactor=real value is expected                     !!!   ***\n");
                Msg_AsFem_Exit();
            }
            else{
                if(numbers[0]<1.0e-12){
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: invalid cutfactor= found in [timestepping] block     !!!   ***\n");
                    PetscPrintf(PETSC_COMM_WORLD,"***        cutfactor=positive real value is expected            !!!   ***\n");
                    Msg_AsFem_Exit();
                }
                timesteppingblock._CutFactor=numbers[0];
            }
        }
        else if(str.find("growthfactor=")!=string::npos||
                str.find("GROWTHFACTOR=")!=string::npos){
            if(!HasType){
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: no growthfactor= found in [timestepping] block       !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        growthfactor= must be given after type=              !!!   ***\n");
                Msg_AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=RemoveStrSpace(substr);
            numbers=SplitStrNum(substr);
            if(numbers.size()<1){
                Msg_Input_LineError(linenum);
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: no growthfactor= found in [timestepping] block       !!!   ***\n");
                PetscPrintf(PETSC_COMM_WORLD,"***        growthfactor=real value is expected                  !!!   ***\n");
                Msg_AsFem_Exit();
            }
            else{
                if(numbers[0]<1.0e-12){
                    Msg_Input_LineError(linenum);
                    PetscPrintf(PETSC_COMM_WORLD,"*** Error: invalid growthfactor= found in [timestepping] block  !!!   ***\n");
                    PetscPrintf(PETSC_COMM_WORLD,"***        growthfactor=positive real value is expected         !!!   ***\n");
                    Msg_AsFem_Exit();
                }
                timesteppingblock._GrowthFactor=numbers[0];
            }
        }
        else if(str.find("[]")!=string::npos){
            Msg_Input_LineError(linenum);
            Msg_Input_BlockBracketNotComplete();
            Msg_AsFem_Exit();
        }
        else{
            Msg_Input_LineError(linenum);
            PetscPrintf(PETSC_COMM_WORLD,"*** Error: unknown option in [timestepping] block               !!!   ***\n");
            Msg_AsFem_Exit();
        }
        getline(in,str);linenum+=1;
    }
    if(!HasType){
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: type= not found in [timestepping] block              !!!   ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***        type=time-integration-method should be given         !!!   ***\n");
        Msg_AsFem_Exit();
    }

    if(!HasEndTime){
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: endtime= not found in [timestepping] block           !!!   ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***        endtime=final time should be given                   !!!   ***\n");
        Msg_AsFem_Exit();
    }

    if(!HasDt){
        PetscPrintf(PETSC_COMM_WORLD,"*** Warning: dt= not found in [timestepping] block              !!!   ***\n");
        PetscPrintf(PETSC_COMM_WORLD,"***          dt=default value will be                           !!!   ***\n");
    }

    return HasType;
}