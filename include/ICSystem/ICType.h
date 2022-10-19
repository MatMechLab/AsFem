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
//+++ Date   : 2020.07.10
//+++ Purpose: Define the initial condition type in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

/**
 * This enum class list out all supported initial conditions in AsFem
 */
enum class ICType{
    NULLIC,
    CONSTIC,
    RANDOMIC,
    RECTANGLEIC,
    CUBICIC,
    CIRCLEIC,
    SMOOTHCIRCLEIC,
    SPHERICALIC,
    USER1IC,
    USER2IC,
    USER3IC,
    USER4IC,
    USER5IC,
    USER6IC,
    USER7IC,
    USER8IC,
    USER9IC,
    USER10IC
};