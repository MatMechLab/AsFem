//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.08.06
//+++ Purpose: Implement the reader for [postprocess] block
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "InputSystem/SingleBlockReader.h"

#include "Postprocess/Postprocess.h"
#include "DofHandler/DofHandler.h"

/**
 * This class implement the reader for [postprocess] block
 */
class PostprocessBlockReader:public SingleBlockReader{
public:
    /**
     * Implement the helper function for [mates] block, if users dont know what to do, then
     * 'type=helper' can show you the basic information for this block
     */
    virtual void PrintHelper() override;

    /**
     * This function responsible for reading the [postprocess] block, it could contains several sub block<br>
     * The basic structure of this block should look like: <br>
     * <pre>
     * [postprocess]
     *   [processname-1]
     *     type=postprocess-type-1
     *     dof=dof1
     *     projvariable=name
     *     scalarmate=name
     *     vectormate=name
     *     rank2mate=name
     *     side=side-name
     *     domain=domain-name
     *     elmtid=id
     *     nodeid=id
     *     iindex=id
     *     jindex=id
     *     component=component
     *   [end]
     *   ...
     * [end]
     * </pre>
     * @param in the ifstream for input file reading
     * @param str the string variable constains '[mates]' line
     * @param lasteendlinenum the line number of the last '[end]' block of [maates]
     * @param linenum the current line number, which should be update during the file reading
     * @param postprocess the postprocess system, which store the params we read and set up the postprocess properties calculation
     * @param dofHandler DofHandler class, which is responsible for the DoFs related calculation
     */
    bool ReadPostprocessBlock(ifstream &in,string str,const int &lastendlinenum,int &linenum,Postprocess &postprocess,DofHandler &dofHandler);
};
