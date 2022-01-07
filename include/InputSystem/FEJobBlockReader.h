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
//+++ Date   : 2021.08.05
//+++ Purpose: Implement the reader for [job] block
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "InputSystem/SingleBlockReader.h"
#include "FEProblem/FEJobBlock.h"

class FEJobBlockReader:public SingleBlockReader{
public:
    /**
     * Implement the helper function for [job] block, if users dont know what to do, then
     * 'type=helper' can show you the basic information for this block
     */
    virtual void PrintHelper() override;


    /**
     * This function responsible for reading the [mates] block, it could contains several sub block<br>
     * The basic structure of this block should look like: <br>
     * <pre>
     * [mates]
     *   type=static,transient
     *   debug=true,false,dep
     * [end]
     * </pre>
     * @param in the ifstream for input file reading
     * @param str the string variable constains '[mates]' line
     * @param linenum the current line number, which should be update during the file reading
     * @param feJobBlock the fe analysis job block class, which store the params we read and set up the FE analysis job
     */
    bool ReadFEJobBlock(ifstream &in,string str,int &linenum,FEJobBlock &feJobBlock);

};
