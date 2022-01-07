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
//+++ Purpose: Implement the reader for [bcs] block
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "InputSystem/SingleBlockReader.h"

#include "BCSystem/BCSystem.h"

/**
 * This class implement the reader for '[bcs]' block, this block could contain
 * multiple sub block
 */
class BCsBlockReader:public SingleBlockReader{

public:
    /**
     * Print the help information for '[bcs]' block,
     * if users dont know what to do, 'type=helper' will
     * show your the basic information
     */
    virtual void PrintHelper() override;

    /**
     * This function can read the [bcs] block information<br>
     * The basic structure of the [bcs] block should look like:
     * <pre>
     * [bcs]
     *   [mybc1]
     *     type=dirichlet,neumann,user1,user2,...
     *     dof=dof-name-1
     *     value=boundary-value-1
     *     boundary=boundary-name-1
     *   [end]
     *   [mybc2]
     *     type=dirichlet, neumann, user1, user2, ...
     *     dof=dof-name-2
     *     value=boundary-value-2
     *     boundary=boundary-name-2
     *   [end]
     *   ...
     * [end]
     * </pre>
     * @param in ifstream for the file reading
     * @param str the string contain '[bcs]' line
     * @param lastendlinenum the last line number of [bcs] block
     * @param linenum the current line number 
     * @param bcSystem the bc class, which will set up the boundary system and the related dofs name, and so on...
     */
    bool ReadBCBlock(ifstream &in,string str,const int &lastendlinenum,int &linenum,BCSystem &bcSystem,DofHandler &dofHandler);


};
