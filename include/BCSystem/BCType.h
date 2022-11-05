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
//+++ Date   : 2020.07.01
//+++ Purpose: Define the boundary condition type in AsFem
//+++          For example: Dirichlet BC, Neumann BC, and so on
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

enum class BCType{
    NULLBC,
    DIRICHLETBC,
    CYCLICDIRICHLETBC,
    USER1DIRICHLETBC,
    USER2DIRICHLETBC,
    USER3DIRICHLETBC,
    USER4DIRICHLETBC,
    USER5DIRICHLETBC,
    POISSON2DBENCHMARKBC,
    ROTATEDDIRICHLETBC,
    //**************************
    NEUMANNBC,
    PRESSUREBC,
    TRACTIONBC,
    NODALDIRICHLETBC,
    NODALNEUMANNBC,
    NODALFORCEBC,
    NODALFLUXBC,
    // for user-defined boundary conditions
    USER1BC,
    USER2BC,
    USER3BC,
    USER4BC,
    USER5BC,
    USER6BC,
    USER7BC,
    USER8BC,
    USER9BC,
    USER10BC
};
