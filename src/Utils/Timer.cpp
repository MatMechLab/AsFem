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
//+++ Date   : 2022.05.09
//+++ Purpose: timer for time elapse counting
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Utils/Timer.h"

Timer::Timer(){
    MPI_Comm_rank(PETSC_COMM_WORLD,&m_rank);
    if(m_rank==0){
        m_start_timer=std::chrono::high_resolution_clock::now();
        m_end_timer=std::chrono::high_resolution_clock::now();
        m_duration_seconds=0.0;
    }
}
void Timer::resetTimer(){
    MPI_Comm_rank(PETSC_COMM_WORLD,&m_rank);
    if(m_rank==0){
        m_start_timer=std::chrono::high_resolution_clock::now();
        m_end_timer=std::chrono::high_resolution_clock::now();
        m_duration_seconds=0.0;
    }
}

void Timer::startTimer(){
    MPI_Comm_rank(PETSC_COMM_WORLD,&m_rank);
    if(m_rank==0){
        m_start_timer=std::chrono::high_resolution_clock::now();
    }
}
void Timer::endTimer(){
    MPI_Comm_rank(PETSC_COMM_WORLD,&m_rank);
    if(m_rank==0){
        m_end_timer=std::chrono::high_resolution_clock::now();
    }
}
void Timer::printElapseTime(const string &instr,const bool &flag)const{
    char buff[16];
    string str;
    snprintf(buff,16,"%14.5e",getDurationInSecond());
    str=buff;
    if(flag) MessagePrinter::printStars();
    MessagePrinter::printNormalTxt(instr+", elapsed time="+str+" [s]");
    MessagePrinter::printStars();
}