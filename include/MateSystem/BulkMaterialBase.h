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
//+++ Date   : 2021.04.04
//+++ Purpose: Define the base material abstract class, all the
//+++          materials should inherit this class!
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once


#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

// For AsFem's own header files
#include "MateSystem/Materials.h"
#include "ElmtSystem/LocalElmtData.h"

using namespace std;

/**
 * The abstract bulk material base class, which defines the common functions for 
 * material property calculation
 */
class BulkMaterialBase{
public:
    /**
     * Initial the preset material properties, if you don't need the history information of some materials, then you can avoid calling this function
     * @param InputParams the input material parameters read from the input file
     * @param elmtinfo the data structure for the local element information
     * @param elmtsoln the solution like U and V of the local element
     * @param Mate the materials to be initialized
     */
    virtual void InitMaterialProperties(const vector<double> &InputParams,
                                        const LocalElmtInfo &elmtinfo,
                                        const LocalElmtSolution &elmtsoln,
                                        Materials &Mate)=0;

    /**
     * Compute the material property accroding to your model
     * @param InputParams the input material parameters read from the input file
     * @param elmtinfo the data structure for the local element information
     * @param elmtsoln the solution like U and V of the local element
     * @param MateOld the materials from previous step
     * @param Mate the materials to be calculated
     */
    virtual void ComputeMaterialProperties(const vector<double> &InputParams,
                                           const LocalElmtInfo &elmtinfo,
                                           const LocalElmtSolution &elmtsoln,
                                           const Materials &MateOld,Materials &Mate)=0;


};
