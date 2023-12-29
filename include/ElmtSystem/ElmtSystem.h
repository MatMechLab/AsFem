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
//+++ Date   : 2022.05.12
//+++ Purpose: the element system in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once


/**
 * for AsFem's headers
 */
#include "ElmtSystem/BulkElmtSystem.h"


/**
 * This class implement the calculation of bulk and interface elements
 */
class ElmtSystem:public BulkElmtSystem{
public:
    /**
     * constructor
     */
    ElmtSystem();

    /**
     * print the information of element system
     */
    void printElmtSystemInfo()const;

};