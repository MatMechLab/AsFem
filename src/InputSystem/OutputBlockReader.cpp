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
//+++ Date   : 2021.08.05
//+++ Purpose: Implement the reader for [output] block
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/OutputBlockReader.h"


void OutputBlockReader::PrintHelper(){
    
    MessagePrinter::PrintStars(MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("The complete information for [output] block:",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("[output]",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("  type=vtk",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("  folder=foldername[default is empty]",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("  interval=5",MessageColor::BLUE);
    MessagePrinter::PrintNormalTxt("[end]",MessageColor::BLUE);
    MessagePrinter::PrintStars(MessageColor::BLUE);

}

//***************************************************
bool OutputBlockReader::ReadOutputBlock(ifstream &in,string str,int &linenum,OutputSystem &outputSystem){
    
    bool HasType;
    vector<double> numbers;
    OutputBlock outputblock;
    string msg;
    // now str already contains [dofs]
    getline(in,str);linenum+=1;
    str=StringUtils::StrToLower(str);
    str=StringUtils::RemoveStrSpace(str);
    numbers.clear();

    HasType=false;

    outputblock.Init();

    while(str.find("[end]")==string::npos&&
          str.find("[END]")==string::npos){
        if(StringUtils::IsCommentLine(str)||str.length()<1){
            getline(in,str);linenum+=1;
            str=StringUtils::StrToLower(str);
            continue;
        }
        str=StringUtils::RemoveStrSpace(str);

        if(str.find("helper")!=string::npos){
            PrintHelper();
            return false;
        }
        else if(str.find("type=")!=string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=StringUtils::RemoveStrSpace(substr);
            if(substr.find("vtu")!=string::npos&&substr.size()==3){
                   outputblock._OutputType=OutputType::VTU;
                   outputblock._OutputFormatName="vtu";
                   HasType=true;
            }
            else if(substr.find("vtk")!=string::npos&&substr.size()==3){
                   outputblock._OutputType=OutputType::VTK;
                   outputblock._OutputFormatName="vtk";
                   HasType=true;
            }
            else if(substr.find("csv")!=string::npos&&substr.size()==3){
                   outputblock._OutputType=OutputType::CSV;
                   outputblock._OutputFormatName="csv";
                   HasType=true;
            }
            else{
                HasType=false;
                MessagePrinter::PrintErrorInLineNumber(linenum);
                msg="unsupported output file type in the [output] block, type= vtu[vtk,csv] is expected";
                MessagePrinter::PrintErrorTxt(msg);
                MessagePrinter::AsFem_Exit();
            }
        }
        else if(str.find("interval=")!=string::npos){
            numbers=StringUtils::SplitStrNum(str);
            if(numbers.size()<1){
                MessagePrinter::PrintErrorInLineNumber(linenum);
                msg="unsupported output interval= option in the [output] block, option=integer is expected";
                MessagePrinter::PrintErrorTxt(msg);
                MessagePrinter::AsFem_Exit();
            }
            else{
                outputblock._Interval=static_cast<int>(numbers[0]);
            }
        }
        else if(str.find("folder=")!=string::npos){
            int i=str.find_first_of('=');
            string substr=str.substr(i+1,str.length());
            substr=StringUtils::RemoveStrSpace(substr);
            if(substr.size()<1){
                outputblock._OutputFolderName.clear();
            }
            else{
                outputblock._OutputFolderName=substr;
            }
        }
        else if(str.find("[]")!=string::npos){
            MessagePrinter::PrintErrorInLineNumber(linenum);
            msg="the bracket pair in the [output] block is not complete";
            MessagePrinter::PrintErrorTxt(msg);
            MessagePrinter::AsFem_Exit();
        }
        else{
            MessagePrinter::PrintErrorInLineNumber(linenum);
            msg="unknown option in the [output] block";
            MessagePrinter::PrintErrorTxt(msg);
            MessagePrinter::AsFem_Exit();
        }
        getline(in,str);linenum+=1;
        str=StringUtils::StrToLower(str);
    }
    HasType=true;
    outputSystem.InitFromOutputBlock(outputblock);
    outputSystem.SetInputFileName(_InputFileName);

    return HasType;
}
