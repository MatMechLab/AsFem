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
//+++ Date   : 2022.04.20
//+++ Purpose: this class defines the basic mesh type in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

/**
 * the enum class for the different type of mesh supported in AsFem
 */
enum class MeshType{
    NULLTYPE,
    // for 0d case
    POINT,
    // for 1d case
    EDGE2,
    EDGE3,
    EDGE4,
    EDGE5,
    // for 2d case
    TRI3,
    TRI6,
    QUAD4,
    QUAD8,
    QUAD9,
    // for 3d case
    TET4,
    TET10,
    HEX8,
    HEX20,
    HEX27
};