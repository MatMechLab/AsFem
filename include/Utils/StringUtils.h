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
//+++ Date   : 2020.06.30
//+++ Purpose: Implement a general class for string manipulate
//+++          For example, string cases convert and the split..
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <string>
#include <algorithm>
#include <sstream>
#include <vector>

using std::vector;
using std::string;

/**
 * This class implements the general manipulation of strings.
 */
class StringUtils{
public:
    /**
     * constructor
     */
    StringUtils();

    /**
     * convert the string to lower case
     * @param instr the input string
     */
    static string strToLowerCase(string instr);

    /**
     * convert the string to upper case
     * @param instr the input string
     */
    static string strToUpperCase(string instr);

    /**
     * remove the empty space in the input str(not modified)
     * @param instr input string, its value will be unchanged
     */
    static string removeStrSpace(string instr);
    /**
     * split the float number from a given string
     * @param instr the input string
     */
    static vector<double> splitStrNum(string instr);

    /**
     * check if an expression is valid for a time-dependent case
     */
    static bool isValidTimeDependentExpression(string instr);

};