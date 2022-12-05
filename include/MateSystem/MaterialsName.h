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
//+++ Date   : 2020.12.20
//+++ Purpose: Define the built-in material type's short name in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <map>
#include <vector>
#include <string>

#include "MathUtils/Vector3d.h"
#include "MathUtils/Rank2Tensor.h"
#include "MathUtils/Rank4Tensor.h"

using std::map;
using std::string;
using std::vector;
//*******************************************
//*** for type redefine, in order to use its short name
typedef map<string,bool>        BooleanMateType;/**< for boolean materials, this material cant be exported */
typedef map<string,double>      ScalarMateType;/**< for scalar materials */
typedef map<string,Vector3d>    VectorMateType;/**< for vector materials */
typedef map<string,Rank2Tensor> Rank2MateType;/**< for rank-2 materials */
typedef map<string,Rank4Tensor> Rank4MateType;/**< for rank-4 materials */