//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.07.12
//+++ Purpose: implement the general FE space for FEM calculation
//+++          in AsFem, here one can use:
//+++            1) gauss integration 
//+++            2) shape functions for different mesh
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once


#include "FE/QPoint.h"


class FE:public QPoint{
public:
    FE();

private:
    
};