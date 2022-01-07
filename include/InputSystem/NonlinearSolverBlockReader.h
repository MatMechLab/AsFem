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
//+++ Purpose: Implement the reader for [nonlinearsolver] block
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "InputSystem/SingleBlockReader.h"

#include "NonlinearSolver/NonlinearSolver.h"

/**
 * This class implement the reader for [nonlinearsolver] block
 */
class NonlinearSolverBlockReader:public SingleBlockReader{
public:
    /**
     * Implement the helper function for [nonlinearsolver] block, if users dont know what to do, then
     * 'type=helper' can show you the basic information for this block
     */
    virtual void PrintHelper() override;

    /**
     * This function responsible for reading the [mates] block, it could contains several sub block<br>
     * The basic structure of this block should look like: <br>
     * <pre>
     * [nonlinearsolver]
     *   type=nr
     *   maxiters=maximum-iteration-number
     *   r_rel_tol=relative-tolerance-of-residual
     *   r_abs_tol=absolute-tolerance-of-residual
     *   stol=tolerance-of-delta-U
     *   solver=superlu,mumps
     * [end]
     * </pre>
     * @param in the ifstream for input file reading
     * @param str the string variable constains '[mates]' line
     * @param linenum the current line number, which should be update during the file reading
     * @param nonlinearSolver the nonlinear solver system, which store the params we read and set up the nonlinear solver
     */
    bool ReadNonlinearSolverBlock(ifstream &in,string str,int &linenum,NonlinearSolver &nonlinearSolver);

private:
    NonlinearSolverBlock _nonlinearSolverBlock;

};
