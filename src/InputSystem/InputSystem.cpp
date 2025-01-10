//****************************************************************
//* This file is part of the AsFem framework
//* Advanced Simulation kit based on Finite Element Method (AsFem)
//* All rights reserved, Yang Bai/MM-Lab@CopyRight 2020-present
//* https://github.com/MatMechLab/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2022.05.07
//+++ Purpose: the input file reading system in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/InputSystem.h"

InputSystem::InputSystem(){
    m_InputFileName.clear();
    m_MeshFileName.clear();
    m_Json.clear();
    m_ReadOnly=false;
}
InputSystem::~InputSystem(){
    m_InputFileName.clear();
    m_MeshFileName.clear();
    m_Json.clear();
}
InputSystem::InputSystem(int args,char *argv[]){
    m_ReadOnly=false;
    if(args==1){
        // ./asfem or asfem
        m_InputFileName.clear();
        m_MeshFileName.clear();
        m_ReadOnly=false;
    }
    else if(args==3){
        // ./asfem -i input.json or asfem -i input.json
        m_ReadOnly=false;
        if(string(argv[1]).find("-i")!=string::npos){
            if(string(argv[2]).size()<5){
                MessagePrinter::printErrorTxt("invalid input file name after '-i', it must be xxx.json");
                MessagePrinter::exitAsFem();
            }
            m_InputFileName=argv[2];
            m_HasInputFile=true;
            if(m_InputFileName.compare(m_InputFileName.size()-5,5,".json")!=0){
                MessagePrinter::printErrorTxt("invalid input file name, your input file must have the extension '.json'");
                MessagePrinter::exitAsFem();
            }
        }
        else{
            m_HasInputFile=false;
            MessagePrinter::printErrorTxt("invalid command line args, the second args must be '-i'");
            MessagePrinter::exitAsFem();
        }
    }
    else{
        m_HasInputFile=false;
        if(string(argv[1]).find("-i")!=string::npos){
            m_InputFileName=argv[2];
            m_HasInputFile=true;
            if(m_InputFileName.compare(m_InputFileName.size()-5,5,".json")!=0){
                m_HasInputFile=false;
                MessagePrinter::printErrorTxt("invalid command line args, the second args must be '-i'");
                MessagePrinter::exitAsFem();
            }
        }
        else{
            m_HasInputFile=false;
            MessagePrinter::printErrorTxt("invalid command line args, the second args must be '-i'");
            MessagePrinter::exitAsFem();
        }
        for(int i=4-1;i<args;i++){
            if(string(argv[i]).find("--read-only")!=string::npos){
                m_ReadOnly=true;
            }
        }
    }
}
//*************************************************
void InputSystem::init(int args,char *argv[]){
    m_ReadOnly=false;
    if(args==1){
        // ./asfem or asfem
        m_InputFileName.clear();
        m_MeshFileName.clear();
        m_ReadOnly=false;
    }
    else if(args==3){
        // ./asfem -i input.json or asfem -i input.json
        if(string(argv[1]).find("-i")!=string::npos){
            if(string(argv[2]).size()<5){
                MessagePrinter::printErrorTxt("invalid input file name after '-i', it must be xxx.json");
                MessagePrinter::exitAsFem();
            }
            m_InputFileName=argv[2];
            m_HasInputFile=true;
            if(m_InputFileName.compare(m_InputFileName.size()-5,5,".json")!=0){
                MessagePrinter::printErrorTxt("invalid input file name, your input file must have the extension '.json'");
                MessagePrinter::exitAsFem();
            }
        }
        else{
            m_HasInputFile=false;
            MessagePrinter::printErrorTxt("invalid command line args, the second args must be '-i'");
            MessagePrinter::exitAsFem();
        }
    }
    else{
        m_HasInputFile=false;
        if(string(argv[1]).find("-i")!=string::npos){
            m_InputFileName=argv[2];
            m_HasInputFile=true;
            if(m_InputFileName.compare(m_InputFileName.size()-5,5,".json")!=0){
                m_HasInputFile=false;
                MessagePrinter::printErrorTxt("invalid command line args, the second args must be '-i'");
                MessagePrinter::exitAsFem();
            }
        }
        else{
            m_HasInputFile=false;
            MessagePrinter::printErrorTxt("invalid command line args, the second args must be '-i'");
            MessagePrinter::exitAsFem();
        }
        for(int i=4-1;i<args;i++){
            if(string(argv[i]).find("--read-only")!=string::npos){
                m_ReadOnly=true;
            }
        }
    }
}