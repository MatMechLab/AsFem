//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
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