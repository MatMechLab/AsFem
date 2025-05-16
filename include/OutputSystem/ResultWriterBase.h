//****************************************************************
//* This file is part of the AsFem framework
//* Advanced Simulation kit based on Finite Element Method (AsFem)
//* All rights reserved, Yang Bai/MM-Lab@CopyRight 2020-present
//* https://github.com/MatMechLab/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2022.08.22
//+++ Purpose: Defines the abstract class for result output
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once


#include "petsc.h"
#include "Utils/MessagePrinter.h"

/**
 * For AsFem's built-in classes
 */
#include "FECell/FECell.h"
#include "DofHandler/DofHandler.h"
#include "SolutionSystem/SolutionSystem.h"
#include "ProjectionSystem/ProjectionSystem.h"



/**
 * This class defines the abstract class for result output
 */
class ResultWriterBase{
public:
    /**
     * save result to different files according to the output format
     * @param t_filename the string name of result file
     * @param t_fecell the fe cell class
     * @param t_dofHandler the dofhandler class
     * @param t_solution the solution class
     * @param t_projection the projection class
     */
    virtual void saveResults(const string &t_filename,
                             const FECell &t_fecell,
                             const DofHandler &t_dofHandler,
                             SolutionSystem &t_solution,
                             ProjectionSystem &t_projection)=0;

};