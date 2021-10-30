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
//+++ Date   : 2021.09.03
//+++ Purpose: Implement the user defined material properties 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "MateSystem/BulkMaterialBase.h"

/**
 * This class calculate the user-defined-material(umat4)
 */
class User4Material:public BulkMaterialBase{
public:
    /**
     * Initialze the material properties
     */
    virtual void InitMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, Materials &Mate) override;

    /**
     * Calculate the material properties
     */
    virtual void ComputeMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, const Materials &MateOld, Materials &Mate) override;

};
