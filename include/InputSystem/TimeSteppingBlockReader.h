//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.08.06
//+++ Purpose: Implement the reader for [timestepping] block
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "InputSystem/SingleBlockReader.h"

#include "TimeStepping/TimeStepping.h"

/**
 * This class implement the reader for [timestepping] block
 */
class TimeSteppingBlockReader:public SingleBlockReader{
public:
    /**
     * Implement the helper function for [timestepping] block, if users dont know what to do, then
     * 'type=helper' can show you the basic information for this block
     */
    virtual void PrintHelper() override;

    /**
     * This function responsible for reading the [timestepping] block, it could contains several sub block<br>
     * The basic structure of this block should look like: <br>
     * <pre>
     * [timestepping]
     *   type=be,cn,...
     *   optiters=3
     *   dt=1.0e-5
     *   dtmin=1.0e-10
     *   dtmax=1.0e-2
     *   time=1.0e1
     *   growthfactor=1.1
     *   cutfactor=0.85
     *   adaptive=true
     * [end]
     * </pre>
     * @param in the ifstream for input file reading
     * @param str the string variable constains '[mates]' line
     * @param lasteendlinenum the line number of the last '[end]' block of [maates]
     * @param linenum the current line number, which should be update during the file reading
     * @param timestepping the time stepping system, which store the params we read and set up the time integration method for transient analysis
     */
    bool ReadTimeSteppingBlock(ifstream &in,string str,int &linenum,TimeStepping &timestepping);


};
