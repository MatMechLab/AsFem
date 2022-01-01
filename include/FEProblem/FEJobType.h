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
//+++ Date   : 2020.12.28
//+++ Purpose: Define the analysis type in AsFem, i.e. static, transient
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

enum class FEJobType{
    STATIC,
    TRANSIENT
};