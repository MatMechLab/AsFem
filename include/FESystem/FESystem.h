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
//+++ Date   : 2020.11.29
//+++ Purpose: Implement the general tasks of FEM calculation for
//+++          both the bulk and interface system.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#pragma once

#include "petsc.h"
/**
 * For AsFem's own header files
 */
#include "FESystem/BulkFESystem.h"


/**
 * This class implements the general FEM calculation actions for 
 * both the bulk and interface element
 */
class FESystem:public BulkFESystem{
public:
    /**
     * constructor
     */
    FESystem();

};