//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "MateSystem/MateSystem.h"

void MateSystem::InitMateSystem(){
    _ScalarMaterials.resize(50,0.0);
    _VectorMaterials.resize(50,Vector3d(0.0));
    _Rank2Materials.resize(10,RankTwoTensor(0.0));
    _Rank4Materials.resize(5,RankFourTensor(0.0));
}