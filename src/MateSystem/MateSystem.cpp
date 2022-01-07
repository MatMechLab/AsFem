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
//+++ Date   : 2020.07.10
//+++ Purpose: Implement the materials system for AsFem
//+++          this class offer some built-in material models
//+++          as well as the User-Defined-Material (umat) models
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MateSystem/MateSystem.h"

MateSystem::MateSystem(){
    _nBulkMateBlocks=0;
    _BulkMateBlockList.clear();

    _Materials.Clean();

    _MaterialsOld.Clean();

}
