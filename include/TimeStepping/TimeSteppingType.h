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
//+++ Date   : 2020.12.29
//+++ Purpose: Define time stepping type for AsFem, i.e.
//+++          backward-euler, cranck-nicolson,...
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#pragma once

enum class TimeSteppingType{
    STATIC,
    BACKWARDEULER,
    CRANCKNICLSON,
    BDF2
};
