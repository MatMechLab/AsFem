//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.07.10
//+++ Purpose: Implement the materials system for AsFem
//+++          this class offer some built-in material models
//+++          as well as the User-Defined-Material (umat) models
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once



//**********************************
//*** For AsFem's own header file
//**********************************
#include "MateSystem/BulkMateSystem.h"

/**
 * This class implement the calculation of user-defined materials for both
 * the bulk domain and interface domain
 */
class MateSystem:public BulkMateSystem{
public:
    MateSystem();

    void PrintMateSystemInfo()const{PrintBulkMateSystemInfo();}
};
