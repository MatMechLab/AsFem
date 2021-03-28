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
//+++ Date   : 2021.02.22
//+++ Purpose: implement the post-process system in AsFem
//+++          this class can do the volume,surface,nodal integration
//+++          or other types post process, and write out the result
//+++          to the csv file
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Postprocess/Postprocess.h"

Postprocess::Postprocess(){
    _PostProcessBlockList.clear();
    _nPostProcessBlocks=0;
    _CSVFileName.clear();
    _VariableNameList.clear();
    _PPSValues.clear();
    _OutputInterval=1;
}
//************************************************************
void Postprocess::PrintPostprocessInfo()const{
    MessagePrinter::PrintNormalTxt("Postprocess system information summary:");
    MessagePrinter::PrintNormalTxt("  output interval="+to_string(GetOutputIntervalNum())
                                      +", number of postprocessor="+to_string(_nPostProcessBlocks));
    for(auto it:_PostProcessBlockList){
        it.PrintInfo();
    }
    MessagePrinter::PrintDashLine();
}
//*****************************************************************
void Postprocess::SetInputFileName(string inputfilename){
    _InputFileName=inputfilename;
    _CSVFileName=_InputFileName.substr(0,_InputFileName.size()-2)+".csv";
}
void Postprocess::AddPostprocessBlock(PostprocessBlock &postprocessblock){
    string msg;
    if(_PostProcessBlockList.size()<1){
        _PostProcessBlockList.push_back(postprocessblock);
        _nPostProcessBlocks=1;
    }
    else{
        bool IsBlockNameUnique;
        IsBlockNameUnique=true;
        for(const auto &it:_PostProcessBlockList){
            if(it._PPSBlockName==postprocessblock._PPSBlockName){
                IsBlockNameUnique=false;
                break;
            }
        }
        if(IsBlockNameUnique){
            _PostProcessBlockList.push_back(postprocessblock);
            _nPostProcessBlocks+=1;
        }
        else{
            msg="duplicated ["+postprocessblock._PPSBlockName+"] is found in [postprocess] sub block, the block name should be unique";
            MessagePrinter::PrintErrorTxt(msg);
            MessagePrinter::AsFem_Exit();
        }
    }
}
//**********************************************************
void Postprocess::InitPPSOutput(){
    if(_nPostProcessBlocks<1){
        MessagePrinter::PrintWarningTxt("no [postprocess] block is found, AsFem will not do any postprocess for you");
        return;
    }
    else{
        _VariableNameList.clear();
        _PPSValues.reserve(_nPostProcessBlocks);
        for(const auto &block:_PostProcessBlockList){
            _VariableNameList.push_back(block._PPSBlockName);
            _PPSValues.push_back(0.0);
        }
        MPI_Comm_rank(PETSC_COMM_WORLD, &_rank);
        if(_rank==0){
            ofstream out;
            out.open(_CSVFileName,ios::out);
            if (!out.is_open()){
                string str="can\'t create a new csv file(="+_CSVFileName+")!, please make sure you have write permission";
                MessagePrinter::PrintErrorTxt(str);
                MessagePrinter::AsFem_Exit();
            }
            out<<"time";
            for(const auto &it:_VariableNameList){
                out<<","<<it;
            }
            out<<endl;
            out.close();
        }
    }
}