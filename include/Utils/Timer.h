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

#pragma once

#include <iostream>
#include <string>
#include <chrono>

/**
 * for AsFem's header
 */
#include "Utils/MessagePrinter.h"

using std::string;

/**
 * this class implement the timer for time elapse estimation
 */
class Timer{
public:
    /**
     * constructor
     */
    Timer();
    /**
     * start the time counting
     */
    void startTimer();
    /**
     * end the time counting
     */
    void endTimer();
    /**
     * reset the timer
     */
    void resetTimer();
    /**
     * print out the elapsed time
     * @param str string, default one is empty string
     * @param flag true for double star lines, false for single star line(bottom)
     */
    void printElapseTime(const string &str="",const bool &flag=true)const;

    //*****************************************************
    //*** general gettings
    //*****************************************************
    /**
     * get the duration time between start and end timer in second
     */
    inline double getDurationInSecond()const{
        if(m_rank==0){
            return std::chrono::duration_cast<std::chrono::milliseconds>(m_end_timer-m_start_timer).count()/1.0e3;
        }
        return 0.0;
    }
    /**
     * get the duration time between start and end timer in minute
     */
    inline double getDurationInMinute()const{
        if(m_rank==0){
            return std::chrono::duration_cast<std::chrono::milliseconds>(m_end_timer-m_start_timer).count()/1.0e3/60.0;
        }
        return 0.0;
    }

private:
    std::chrono::high_resolution_clock::time_point m_start_timer;/**< the start timer */
    std::chrono::high_resolution_clock::time_point m_end_timer;/**< the end timer */
    double m_duration_seconds;/**< duration time in seconds */
    int m_rank;/**< rank number */

};