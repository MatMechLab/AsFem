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
//+++ Date   : 2020.06.28
//+++ Purpose: Implement a general message printer for AsFem
//+++          This class is general to provide almost all the
//+++          message print in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <string>
#include <iomanip>
#include <vector>

#include "petsc.h"

using namespace std;


class MessagePrinter{
public:
    MessagePrinter();

    static void PrintTxt(string str);
    static void PrintShortTxt(string str);
    static void PrintLongTxt(string str);
    
    static void PrintWelcomeTxt(string str);
    static void PrintStars();

private:
    static const int _nWords=75;
    vector<string> SplitStr2Vec(string str);
};