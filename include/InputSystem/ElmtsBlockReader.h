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
//+++ Date   : 2021.07.14
//+++ Purpose: Implement the reader for [elmts] block
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "InputSystem/SingleBlockReader.h"

#include "DofHandler/DofHandler.h"
#include "ElmtSystem/ElmtSystem.h"

/**
 * The class responsible for the [elmts] block reader <br>
 * It should noted that, the [elmts] block could have multiple sub-blocks <br>
 */
class ElmtsBlockReader:public SingleBlockReader{
public:
    /**
     * This function can print out the [elmts] block 
     * template for users
     */
    virtual void PrintHelper() override;

    /**
     * This function will read the [elmts] block,
     * the [elmts] should contain several sub-block
     * for different elements(PDEs/ODEs)
     * It should noted that, the [elmts] block could have multiple sub-blocks <br>
     * For example: <br>
     * <pre>
     * [elmts]
     *   [elmt-1]
     *     type=mechanics
     *     dofs=disp_x disp_y
     *     mate=myelasticmate
     *     domain=mydomain1
     *   [end]
     *   [elmt-2]
     *     type=diffusion
     *     dofs=c
     *     mate=mymate2
     *     domain=mydomain2
     *   [end]
     * [end]
     *</pre>
     * @param in ifstream for input file
     * @param str the string which contains '[elmts]' line 
     * @param lastendlinenum the line number of the the last '[end]' bracket in '[elmts]' block
     * @param linenum the current line number, which will be update during the file reading
     * @param elmtSystem the elmenet system class, wich is used to set up the configuration of elments during the input file reading
     * @param dofHandler the dofHandler class, which is used to check whether the dof's name is valid or not 
     */
    bool ReadElmtsBlock(ifstream &in,string str,const int &lastendlinenum,int &linenum,ElmtSystem &elmtSystem,DofHandler &dofHandler);

};
