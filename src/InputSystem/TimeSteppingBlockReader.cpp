//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.08.06
//+++ Purpose: Implement the reader for [timestepping] block
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/TimeSteppingBlockReader.h"


void TimeSteppingBlockReader::PrintHelper(){
    MessagePrinter::PrintStars(MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("The complete information for [timestepping] block:",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("[timestepping]",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("  type=be",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("  dt=1.0e-6",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("  adaptive=true",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("  optiters=4",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("  dtmin=1.0e-10",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("  dtmax=1.0e-2",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("  growthfactor=1.1",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("  cutfactor=0.85",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("[end]",MessageColor::BLUE);
    MessagePrinter::PrintStars(MessageColor::BLUE);
}
//**********************************************************
bool TimeSteppingBlockReader::ReadTimeSteppingBlock(ifstream &in, string str, int &linenum, TimeStepping &timestepping){
    bool HasType=false;
    vector<double> numbers;
    string namelist;
    // now str already contains [dofs]
    getline(in,str);linenum+=1;

    TimeSteppingBlock timesteppingBlock;

    timesteppingBlock.Init();// use the default value

    while(str.find("[end]")==string::npos&&
          str.find("[END]")==string::npos){
        if(StringUtils::IsCommentLine(str)||str.length()<1){
            getline(in,str);linenum+=1;
            str=StringUtils::StrToLower(str);
            str=StringUtils::RemoveStrSpace(str);
            continue;
        }
        if(str.find("type=")!=string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=StringUtils::RemoveStrSpace(substr);
            if(substr.find("helper")!=string::npos){
                PrintHelper();
                return false;
            }
            else if((substr.find("be")!=string::npos||substr.find("BE")!=string::npos)&&
               substr.length()==2){
                HasType=true;
                timesteppingBlock._TimeSteppingType=TimeSteppingType::BACKWARDEULER;
                timesteppingBlock._TimeSteppingTypeName="backward-euler";
            }
            else if((substr.find("cn")!=string::npos||substr.find("CN")!=string::npos)&&
               substr.length()==2){
                HasType=true;
                timesteppingBlock._TimeSteppingType=TimeSteppingType::CRANCKNICLSON;
                timesteppingBlock._TimeSteppingTypeName="crank-nicolson";
            }
            else if((substr.find("gl")!=string::npos||substr.find("GL")!=string::npos)&&
               substr.length()==2){
                HasType=true;
                timesteppingBlock._TimeSteppingType=TimeSteppingType::GL;
                timesteppingBlock._TimeSteppingTypeName="general-linear";
            }
            else if((substr.find("rosw")!=string::npos||substr.find("ROSW")!=string::npos)&&
               substr.length()==4){
                HasType=true;
                timesteppingBlock._TimeSteppingType=TimeSteppingType::ROSW;
                timesteppingBlock._TimeSteppingTypeName="rosenbrock-w";
            }
            else if((substr.find("bdf2")!=string::npos||substr.find("BDF2")!=string::npos)&&
               substr.length()==4){
                HasType=true;
                timesteppingBlock._TimeSteppingType=TimeSteppingType::BDF2;
                timesteppingBlock._TimeSteppingTypeName="back-difference-formula2";
            }
            else{
                MessagePrinter::PrintErrorInLineNumber(linenum);
                MessagePrinter::PrintErrorTxt("unsupported type in the [timestepping] block");
                MessagePrinter::AsFem_Exit();
            }
        }
        else if(str.find("optiters")!=string::npos||
            str.find("OptIters")!=string::npos||
            str.find("OPTITERS")!=string::npos){
            if(!HasType){
                MessagePrinter::PrintErrorInLineNumber(linenum);
                MessagePrinter::PrintErrorTxt("no 'type=' is found in the [timestepping] block, 'optiters=' must be given after 'type='");
                MessagePrinter::AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            numbers=StringUtils::SplitStrNum(substr);
            if(numbers.size()<1){
                MessagePrinter::PrintErrorInLineNumber(linenum);
                MessagePrinter::PrintErrorTxt("optiters= number can not be found in the [timestepping] block, optiters=integer is expected");
                MessagePrinter::AsFem_Exit();
            }
            else{
                if(int(numbers[0])<1){
                    MessagePrinter::PrintErrorInLineNumber(linenum);
                    MessagePrinter::PrintErrorTxt("invalid optiters number in the [timestepping] block, optiters=integer is expected",false);
                    MessagePrinter::AsFem_Exit();
                }
                timesteppingBlock._OptIters=int(numbers[0]);
            }
        }
        else if(str.find("dt=")!=string::npos||
                 str.find("Dt=")!=string::npos||
                 str.find("DT=")!=string::npos){
            if(!HasType){
                MessagePrinter::PrintErrorTxt("no 'type=' found in the [timestepping] block, dt= should be given after 'type='");
                MessagePrinter::AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            numbers=StringUtils::SplitStrNum(substr);
            if(numbers.size()<1){
                MessagePrinter::PrintErrorInLineNumber(linenum);
                MessagePrinter::PrintErrorTxt("no dt found in the [timestepping] block, dt=real should be given");
                MessagePrinter::AsFem_Exit();
            }
            else{
                if(numbers[0]<=0.0){
                    MessagePrinter::PrintErrorInLineNumber(linenum);
                    MessagePrinter::PrintErrorTxt("invalid dt found in [timestepping] block, dt=real should be given");
                    MessagePrinter::AsFem_Exit();
                }
                timesteppingBlock._Dt=numbers[0];
            }
        }
        else if(str.find("time=")!=string::npos||
                 str.find("Time=")!=string::npos||
                 str.find("TIME=")!=string::npos){
            if(!HasType){
                MessagePrinter::PrintErrorTxt("no 'type=' found in the [timestepping] block, time= should be given after 'type='");
                MessagePrinter::AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            numbers=StringUtils::SplitStrNum(substr);
            if(numbers.size()<1){
                MessagePrinter::PrintErrorInLineNumber(linenum);
                MessagePrinter::PrintErrorTxt("no time= found in the [timestepping] block, time=real should be given");
                MessagePrinter::AsFem_Exit();
            }
            else{
                if(numbers[0]<=0.0){
                    MessagePrinter::PrintErrorInLineNumber(linenum);
                    MessagePrinter::PrintErrorTxt("invalid time found in [timestepping] block, time=real should be given");
                    MessagePrinter::AsFem_Exit();
                }
                timesteppingBlock._FinalT=numbers[0];
            }
        }
        else if(str.find("dtmin=")!=string::npos||
                str.find("DtMin=")!=string::npos||
                str.find("DTMIN=")!=string::npos){
            if(!HasType){
                MessagePrinter::PrintErrorTxt("no 'type=' found in the [timestepping] block, dtmin= should be given after 'type='");
                MessagePrinter::AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            numbers=StringUtils::SplitStrNum(substr);
            if(numbers.size()<1){
                MessagePrinter::PrintErrorInLineNumber(linenum);
                MessagePrinter::PrintErrorTxt("no dtmin found in the [timestepping] block, dtmin=real should be given");
                MessagePrinter::AsFem_Exit();
            }
            else{
                if(numbers[0]<=0.0){
                    MessagePrinter::PrintErrorInLineNumber(linenum);
                    MessagePrinter::PrintErrorTxt("invalid dtmin found in [timestepping] block, dtmin=real should be given");
                    MessagePrinter::AsFem_Exit();
                }
                timesteppingBlock._DtMin=numbers[0];
            }
        }
        else if(str.find("dtmax=")!=string::npos||
                str.find("DtMax=")!=string::npos||
                str.find("DTMAX=")!=string::npos){
            if(!HasType){
                MessagePrinter::PrintErrorTxt("no 'type=' found in the [timestepping] block, dtmax= should be given after 'type='");
                MessagePrinter::AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            numbers=StringUtils::SplitStrNum(substr);
            if(numbers.size()<1){
                MessagePrinter::PrintErrorInLineNumber(linenum);
                MessagePrinter::PrintErrorTxt("no dtmax found in the [timestepping] block, dtmax=real should be given");
                MessagePrinter::AsFem_Exit();
            }
            else{
                if(numbers[0]<=0.0){
                    MessagePrinter::PrintErrorInLineNumber(linenum);
                    MessagePrinter::PrintErrorTxt("invalid dtmax found in [timestepping] block, dtmax=real should be given");
                    MessagePrinter::AsFem_Exit();
                }
                timesteppingBlock._DtMax=numbers[0];
            }
        }
        else if(str.find("growthfactor=")!=string::npos||
                 str.find("GrowthFactor=")!=string::npos||
                 str.find("GROWTHFACTOR=")!=string::npos){
            if(!HasType){
                MessagePrinter::PrintErrorTxt("no 'type=' found in the [timestepping] block, growthfactor= should be given after 'type='");
                MessagePrinter::AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            numbers=StringUtils::SplitStrNum(substr);
            if(numbers.size()<1){
                MessagePrinter::PrintErrorInLineNumber(linenum);
                MessagePrinter::PrintErrorTxt("no growthfactor found in the [timestepping] block, growthfactor=real should be given");
                MessagePrinter::AsFem_Exit();
            }
            else{
                if(numbers[0]<=0.0){
                    MessagePrinter::PrintErrorInLineNumber(linenum);
                    MessagePrinter::PrintErrorTxt("invalid growthfactor found in [timestepping] block, growthfactor=real should be given");
                    MessagePrinter::AsFem_Exit();
                }
                timesteppingBlock._GrowthFactor=numbers[0];
            }
        }
        else if(str.find("cutfactor=")!=string::npos||
                 str.find("CutFactor=")!=string::npos||
                 str.find("CUTFACTOR=")!=string::npos){
            if(!HasType){
                MessagePrinter::PrintErrorTxt("no 'type=' found in the [timestepping] block, cutfactor= should be given after 'type='");
                MessagePrinter::AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            numbers=StringUtils::SplitStrNum(substr);
            if(numbers.size()<1){
                MessagePrinter::PrintErrorInLineNumber(linenum);
                MessagePrinter::PrintErrorTxt("no cutfactor found in the [timestepping] block, cutfactor=real should be given");
                MessagePrinter::AsFem_Exit();
            }
            else{
                if(numbers[0]<=0.0){
                    MessagePrinter::PrintErrorInLineNumber(linenum);
                    MessagePrinter::PrintErrorTxt("invalid cutfactor found in [timestepping] block, cutfactor=real should be given");
                    MessagePrinter::AsFem_Exit();
                }
                timesteppingBlock._CutBackFactor=numbers[0];
            }
        }
        else if(str.find("adaptive=")!=string::npos||
                str.find("Adaptive=")!=string::npos||
                str.find("ADAPTIVE=")!=string::npos){
            if(!HasType){
                MessagePrinter::PrintErrorTxt("no 'type=' found in the [timestepping] block, adaptive= should be given after 'type='");
                MessagePrinter::AsFem_Exit();
            }
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            if(substr.find("true")!=string::npos||
               substr.find("True")!=string::npos||
               substr.find("TRUE")!=string::npos){
                timesteppingBlock._Adaptive=true;
            }
            else if(substr.find("false")!=string::npos||
               substr.find("False")!=string::npos||
               substr.find("FALSE")!=string::npos){
                timesteppingBlock._Adaptive=false;
            }
            else{
                MessagePrinter::PrintErrorInLineNumber(linenum);
                MessagePrinter::PrintErrorTxt("unknown option for adaptive= in [timestepping] block, true or false is expected");
                MessagePrinter::AsFem_Exit();
            }
        }
        else if(str.find("[]")!=string::npos){
            MessagePrinter::PrintErrorInLineNumber(linenum);
            MessagePrinter::PrintErrorTxt("the bracket pair is not complete in the [timestepping] block",false);
            MessagePrinter::AsFem_Exit();
        }
        else{
            MessagePrinter::PrintErrorInLineNumber(linenum);
            MessagePrinter::PrintErrorTxt("unknown option in [timestepping] block",false);
            MessagePrinter::AsFem_Exit();
        }
        getline(in,str);linenum+=1;
    }

    timestepping.SetOpitonsFromTimeSteppingBlock(timesteppingBlock);

    return HasType;
}
