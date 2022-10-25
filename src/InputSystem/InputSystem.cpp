//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
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
    m_inputfile_name.clear();
    m_meshfile_name.clear();
    m_json.clear();
    m_readonly=false;
}
InputSystem::~InputSystem(){
    m_inputfile_name.clear();
    m_meshfile_name.clear();
    m_json.clear();
}
InputSystem::InputSystem(int args,char *argv[]){
    m_readonly=false;
    if(args==1){
        // ./asfem or asfem
        m_inputfile_name.clear();
        m_meshfile_name.clear();
        m_readonly=false;
    }
    else if(args==3){
        // ./asfem -i input.json or asfem -i input.json
        m_readonly=false;
        if(string(argv[1]).find("-i")!=string::npos){
            if(string(argv[2]).size()<5){
                MessagePrinter::printErrorTxt("invalid input file name after '-i', it must be xxx.json");
                MessagePrinter::exitAsFem();
            }
            m_inputfile_name=argv[2];
            m_hasinputfile=true;
            if(m_inputfile_name.compare(m_inputfile_name.size()-5,5,".json")!=0){
                MessagePrinter::printErrorTxt("invalid input file name, your input file must have the extension '.json'");
                MessagePrinter::exitAsFem();
            }
        }
        else{
            m_hasinputfile=false;
            MessagePrinter::printErrorTxt("invalid command line args, the second args must be '-i'");
            MessagePrinter::exitAsFem();
        }
    }
    else{
        m_hasinputfile=false;
        if(string(argv[1]).find("-i")!=string::npos){
            m_inputfile_name=argv[2];
            m_hasinputfile=true;
            if(m_inputfile_name.compare(m_inputfile_name.size()-5,5,".json")!=0){
                m_hasinputfile=false;
                MessagePrinter::printErrorTxt("invalid command line args, the second args must be '-i'");
                MessagePrinter::exitAsFem();
            }
        }
        else{
            m_hasinputfile=false;
            MessagePrinter::printErrorTxt("invalid command line args, the second args must be '-i'");
            MessagePrinter::exitAsFem();
        }
        for(int i=4-1;i<args;i++){
            if(string(argv[i]).find("--read-only")!=string::npos){
                m_readonly=true;
            }
        }
    }
}
//*************************************************
void InputSystem::init(int args,char *argv[]){
    m_readonly=false;
    if(args==1){
        // ./asfem or asfem
        m_inputfile_name.clear();
        m_meshfile_name.clear();
        m_readonly=false;
    }
    else if(args==3){
        // ./asfem -i input.json or asfem -i input.json
        if(string(argv[1]).find("-i")!=string::npos){
            if(string(argv[2]).size()<5){
                MessagePrinter::printErrorTxt("invalid input file name after '-i', it must be xxx.json");
                MessagePrinter::exitAsFem();
            }
            m_inputfile_name=argv[2];
            m_hasinputfile=true;
            if(m_inputfile_name.compare(m_inputfile_name.size()-5,5,".json")!=0){
                MessagePrinter::printErrorTxt("invalid input file name, your input file must have the extension '.json'");
                MessagePrinter::exitAsFem();
            }
        }
        else{
            m_hasinputfile=false;
            MessagePrinter::printErrorTxt("invalid command line args, the second args must be '-i'");
            MessagePrinter::exitAsFem();
        }
    }
    else{
        m_hasinputfile=false;
        if(string(argv[1]).find("-i")!=string::npos){
            m_inputfile_name=argv[2];
            m_hasinputfile=true;
            if(m_inputfile_name.compare(m_inputfile_name.size()-5,5,".json")!=0){
                m_hasinputfile=false;
                MessagePrinter::printErrorTxt("invalid command line args, the second args must be '-i'");
                MessagePrinter::exitAsFem();
            }
        }
        else{
            m_hasinputfile=false;
            MessagePrinter::printErrorTxt("invalid command line args, the second args must be '-i'");
            MessagePrinter::exitAsFem();
        }
        for(int i=4-1;i<args;i++){
            if(string(argv[i]).find("--read-only")!=string::npos){
                m_readonly=true;
            }
        }
    }
}