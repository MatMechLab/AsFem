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
//+++ Date   : 2020.07.01
//+++ Purpose: Define the element system in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <iomanip>

#include "Utils/MessagePrinter.h"


#include "ElmtSystem/ElmtBlock.h"
#include "ElmtSystem/BulkElmtSystem.h"


using namespace std;

class ElmtSystem:public BulkElmtSystem{
public:
    ElmtSystem();

    void PrintElmtSystemInfo()const{PrintBulkElmtInfo();}
    //*******************************************
    //*** some getting funs
    //*******************************************


};