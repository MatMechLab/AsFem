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
//+++ Date   : 2020.12.20
//+++ Purpose: Define the built-in material type's short name in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <map>
#include <vector>
#include <string>

#include "Utils/Vector3d.h"
#include "Utils/RankTwoTensor.h"
#include "Utils/RankFourTensor.h"

using namespace std;
//*******************************************
//*** for type redfine in order to use its short name
typedef map<string,double>         ScalarMateType;
typedef map<string,Vector3d>       VectorMateType;
typedef map<string,RankTwoTensor>  Rank2MateType;
typedef map<string,RankFourTensor> Rank4MateType;

