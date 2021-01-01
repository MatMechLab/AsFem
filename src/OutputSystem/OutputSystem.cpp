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
//+++ Date   : 2020.07.12
//+++ Purpose: define the output system for AsFem, where all the 
//+++          results should be written out to the result file
//+++          by this class
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "OutputSystem/OutputSystem.h"


OutputSystem::OutputSystem(){
    _Interval=1;
    _OutputType=OutputType::VTU;
    _OutputTypeName="vtu";
    _OutputFolderName.clear();    
    _OutputFileName.clear();
    _InputFileName.clear();
    _CSVFieldNameList.clear();
}

void OutputSystem::Init(string inputfilename){
    _Interval=1;
    _OutputType=OutputType::VTU;
    _OutputTypeName="vtu";
    _OutputFolderName.clear();    
    _OutputFileName.clear();
    _InputFileName=inputfilename;
    _CSVFieldNameList.clear();
}

void OutputSystem::InitFromOutputBlock(OutputBlock &outputblock){
    _Interval=outputblock._Interval;
    _OutputType=outputblock._OutputType;
    _OutputTypeName=outputblock._OutputFormatName;
    _OutputFolderName=outputblock._OutputFolderName;
}

void OutputSystem::SetOutputType(OutputType outputtype){
    if(outputtype==OutputType::CSV){
        _OutputType=OutputType::CSV;
        _OutputTypeName="csv";
    }
    else if(outputtype==OutputType::VTK){
        _OutputType=OutputType::VTK;
        _OutputTypeName="vtk";
    }
    else if(outputtype==OutputType::VTU){
        _OutputType=OutputType::VTU;
        _OutputTypeName="vtu";
    }
}

//****************************************************
void OutputSystem::PrintInfo()const{
    MessagePrinter::PrintNormalTxt("Output system information summary:");
    MessagePrinter::PrintNormalTxt("  output file format ="+_OutputTypeName);
    MessagePrinter::PrintNormalTxt("  output interval="+to_string(_Interval));
    MessagePrinter::PrintDashLine();
}