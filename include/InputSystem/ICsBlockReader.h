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
//+++ Purpose: Implement the reader for [ics] block
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "InputSystem/SingleBlockReader.h"

#include "DofHandler/DofHandler.h"
#include "ICSystem/ICSystem.h"

/**
 * This class implement the reader for [ics] block in the input file
 */
class ICsBlockReader:public SingleBlockReader{
public:
    /**
     * This function can print out the basic usage of [ics] block.
     * If users dont know what to do , 'type=helper' can show you some demos
     */
    virtual void PrintHelper() override;

    /**
     * This function will read the [ics] block. The basic structure of [ics] block
     * should look like: <br>
     * <pre>
     * [ics]
     *   [ic-1]
     *     type=ic-type-1
     *     dof=dof-name-1
     *     params=val1 val2
     *     domain=domain-name-1
     *   [end]
     *   [ic-2]
     *     type=ic-type-2
     *     dof=dof-name-1
     *     params=val3 val4
     *     domain=domain-name-2
     *   [end]
     *   ...
     * [end]
     * </pre>
     * @param in ifstream for input file
     * @param str the string which contains '[elmts]' line 
     * @param lastendlinenum the line number of the the last '[end]' bracket in '[elmts]' block
     * @param linenum the current line number, which will be update during the file reading
     * @param icSystem the initial condition system class, which is used to set up the basic information in ICs 
     * @param dofHandler the dofHandler class, which is used to check the validation of the applied DoFs
     */ 
    bool ReadICBlock(ifstream &in,string str,const int &lastendlinenum,int &linenum,ICSystem &icSystem,DofHandler &dofHandler);

};
