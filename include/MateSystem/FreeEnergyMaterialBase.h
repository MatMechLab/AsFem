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
//+++ Date   : 2021.04.04
//+++ Purpose: Define the base material for free energy materials
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#pragma once

#include "MateSystem/BulkMaterialBase.h"

/**
 * The base class for free energy materials, where the calculation of free energy and its 1st+2nd order derivatives are required
 */
class FreeEnergyMaterialBase: public BulkMaterialBase{
protected:
    /**
     * Compute the general free energy
     * @param InputParams the material parameters read from the input file
     * @param elmtsoln the solution vector of current element
     * @param F the free energy value vector
     */
    virtual void ComputeF(const vector<double> &InputParams,const LocalElmtSolution &elmtsoln,vector<double> &F)=0;
    
    /**
     * Compute the general free energy's 1st order derivative (chemical potentials)
     * @param InputParams the material parameters read from the input file
     * @param elmtsoln the solution vector of current element
     * @param dF the free energy's 1st derivatives (\f$\mu\f$)
     */
    virtual void ComputedFdU(const vector<double> &InputParams,const LocalElmtSolution &elmtsoln,vector<double> &dF)=0;
    
    /**
     * Compute the general free energy's 2nd order derivative (\f$\partial\mu/\partial c\f$)
     * @param InputParams the material parameters read from the input file
     * @param elmtsoln the solution vector of current element
     * @param d2F the free energy's 1st derivatives (\f$d\mu/dc\f$)
     */
    virtual void Computed2FdU2(const vector<double> &InputParams,const LocalElmtSolution &elmtsoln,vector<double> &d2F)=0;

};
