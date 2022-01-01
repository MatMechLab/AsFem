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
//+++ Purpose: Implement the reader for [projection] block
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "InputSystem/SingleBlockReader.h"

#include "SolutionSystem/SolutionSystem.h"

class ProjectionBlockReader:public SingleBlockReader{
public:
    /**
     * Implement the helper function for [projection] block, if users dont know what to do, then
     * 'type=helper' can show you the basic information for this block
     */
    virtual void PrintHelper() override;

    /**
     * This function responsible for reading the [projection] block
     * The basic structure of this block should look like: <br>
     * <pre>
     * [projection]
     *   type=mate-type-1
     *   name=projection-variable-name
     *   scalarmate=scalar-material-name
     *   vectormate=vector-material-name
     *   rank2mate=rank2-material-name
     *   rank4mate=rank4-material-name
     * [end]
     * </pre>
     * @param in the ifstream for input file reading
     * @param str the string variable constains '[projection]' line
     * @param linenum the current line number, which should be update during the file reading
     * @param solution the solution system, which store the properties we want to project to nodal points
     */
    bool ReadProjectionBlock(ifstream &in,string str,int &linenum,SolutionSystem &solutionSystem);
};
