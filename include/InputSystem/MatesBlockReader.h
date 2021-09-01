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
//+++ Purpose: Implement the reader for [mates] block
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "InputSystem/SingleBlockReader.h"

#include "MateSystem/MateSystem.h"

/**
 * This class responsible for the '[mates]' block reading
 */
class MatesBlockReader:public SingleBlockReader{
public:
    /**
     * Implement the helper function for [mates] block, if users dont know what to do, then
     * 'type=helper' can show you the basic information for this block
     */
    virtual void PrintHelper() override;

    /**
     * This function responsible for reading the [mates] block, it could contains several sub block<br>
     * The basic structure of this block should look like: <br>
     * <pre>
     * [mates]
     *   [mate-1]
     *     type=mate-type-1 
     *     params=val1 val2 
     *   [end]
     *   [mate-2]
     *     type=mate-type-20
     *     params=val3 val4
     *   [end]
     * [end]
     * </pre>
     * @param in the ifstream for input file reading
     * @param str the string variable constains '[mates]' line
     * @param lasteendlinenum the line number of the last '[end]' block of [maates]
     * @param linenum the current line number, which should be update during the file reading
     * @param mateSystem the material system, which store the params we read and set up the material properties calculation
     */ 
    bool ReadMatesBlock(ifstream &in,string str,const int &lastendlinenum,int &linenum,MateSystem &mateSystem);

};
