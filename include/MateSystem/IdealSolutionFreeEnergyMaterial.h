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
//+++ Date   : 2021.08.31
//+++ Purpose: Calculate the free energy, chemical potential and its
//+++          derivatives of ideal solution case
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "MateSystem/FreeEnergyMaterialBase.h"

/**
 * This class calculate the free energy of ideal solution
 */
class IdealSolutionFreeEnergyMaterial: public FreeEnergyMaterialBase{
public:
    /**
     * constructor
     */
    IdealSolutionFreeEnergyMaterial();
    
    
    virtual void InitMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, Materials &Mate) override;


    virtual void ComputeMaterialProperties(const vector<double> &InputParams, const LocalElmtInfo &elmtinfo, const LocalElmtSolution &elmtsoln, const Materials &MateOld, Materials &Mate) override;


private:
    /**
     * Compute the general free energy
     * @param InputParams the material parameters read from the input file
     * @param elmtsoln the solution vector of current element
     * @param F the free energy value vector
     */
    virtual void ComputeF(const vector<double> &InputParams,const LocalElmtSolution &elmtsoln,vector<double> &F) override;
    

    /**
     * Compute the general free energy's 1st order derivative (chemical potentials)
     * @param InputParams the material parameters read from the input file
     * @param elmtsoln the solution vector of current element
     * @param dF the free energy's 1st derivatives (\f$\mu\f$)
     */
    virtual void ComputedFdU(const vector<double> &InputParams,const LocalElmtSolution &elmtsoln,vector<double> &dF) override;
    
    /**
     * Compute the general free energy's 2nd order derivative (\f$\partial\mu/\partial c\f$)
     * @param InputParams the material parameters read from the input file
     * @param elmtsoln the solution vector of current element
     * @param d2F the free energy's 1st derivatives (\f$d\mu/dc\f$)
     */
    virtual void Computed2FdU2(const vector<double> &InputParams,const LocalElmtSolution &elmtsoln,vector<double> &d2F) override;

    

private:
    double c;/**< local concentration*/
    vector<double> _F,_dFdc,_d2Fdc2;/**< local array for F and its derivatives*/

};

