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
//+++ Purpose: Implement the reader for [qpoint] block
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "InputSystem/SingleBlockReader.h"

#include "FE/FE.h"

/**
 * This class implement the reader for [qpoint] block
 */
class QPointBlockReader:public SingleBlockReader{
public:
    /**
     * This function can print out the [elmts] block 
     * template for users
     */
    virtual void PrintHelper() override;

    /**
     * This function will read the [qpoint] block,
     * For example: <br>
     * <pre>
     * [qpoint]
     *   type=gauss,[gausslobatto]
     *   order=3
     *   bcorder=1
     * [end]
     * </pre>
     * @param in ifstream for input file
     * @param str the string which contains '[elmts]' line 
     * @param linenum the current line number, which will be update during the file reading
     * @param fe the FE space class, which is used to set up the basic information of the Guass point integration method  
     */
    bool ReadQPointBlock(ifstream &in,string str,int &linenum,FE &fe);


};
