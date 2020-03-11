//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_BCTYPE_H
#define ASFEM_BCTYPE_H

enum class BCType{
    NullBC,
    PresetBC,
    DirichletBC,
    NeumannBC,
    PressureBC,
    User1BC,
    User2BC,
    User3BC,
    User4BC,
    User5BC,
    User6BC,
    User7BC,
    User8BC,
    User9BC,
    User10BC
};

#endif //ASFEM_BCTYPE_h