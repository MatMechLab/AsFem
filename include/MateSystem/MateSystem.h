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
//+++ Date   : 2020.07.10
//+++ Purpose: Implement the materials system for AsFem.
//+++          This class offer some built-in material models,
//+++          as well as the User-Defined-Material (umat) models
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

/**
 * For AsFem's own headers
 */
#include "MateSystem/BulkMateSystem.h"

/**
 * This class implement the calculation of user-defined materials for both
 * the bulk domain and interface domain
 */
class MateSystem:public BulkMateSystem{
public:
    /**
     * constructor
     */
    MateSystem();

    /**
     * print out the mate system info
     */
    void printMateSystemInfo()const;

};