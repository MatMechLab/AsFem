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
//+++ Date   : 2022.05.13
//+++ Purpose: defines the shape function type in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

/**
 * The enum clas defines the basic shape function type in AsFem.
 * It must be mentioned that, in most of cases, you should not change
 * it to user1~user5, it should be decided by the type of your mesh, namely the 'default' type.
 */
enum class ShapeFunType{
    DEFAULT,
    USER1,
    USER2,
    USER3,
    USER4,
    USER5
};