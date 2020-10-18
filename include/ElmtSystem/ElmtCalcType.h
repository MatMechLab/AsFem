//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.07.18
//+++ Purpose: Define the basic calculation job in element system
//+++          the action like compute residual and compute jacobian
//+++          should be list here
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

enum class BulkElmtCalcType{
    ComputeResidual,
    ComputeJacobian,
    ComputeProjectionValues,
    InitHistoryValues,
    UpdateHistoryValues
};

