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
//+++ Date   : 2021.07.13
//+++ Purpose: Implement the reader for [dofs] block
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "InputSystem/SingleBlockReader.h"

#include "DofHandler/DofHandler.h"

/**
 * This class implement the reader for the [dofs] block in our input file
 */
class DofsBlockReader:public SingleBlockReader{
public:
    /**
     * Print out the complete information for the [dofs] block <br>
     * If you dont know what to do, just using 'type=helper' <br>
     */
    virtual void PrintHelper() override;

    /**
     * Read the [dofs] block, if everything works fine, it should return True value<br>
     * The block should looks like: <br>
     * <pre>
     * [dofs]
     * names=dof1_name dof2_name
     * [end]
     * </pre>
     * @param in ifstream for the file reading
     * @param str the string which constains '[dofs]' line
     * @param linenum the current line number
     * @param dofHandler the dofHandler class, which is used to set up the dofs name we read
     */
    bool ReadDofsBlock(ifstream &in,string str,int &linenum,DofHandler &dofHandler);

};
