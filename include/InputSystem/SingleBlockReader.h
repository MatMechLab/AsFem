//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.07.08
//+++ Purpose: Define the base class for single block reader
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>


#include "petsc.h"

/**
 * include AsFem's own header files
 */
#include "Utils/StringUtils.h"
#include "Utils/MessagePrinter.h"
#include "Utils/MessageColor.h"

using namespace std;


/**
 * Define the base class for single block in our input file,
 * each reader for different single block should inherit this base class
 */
class SingleBlockReader{
protected:
    /**
     * Each single block reader should re-implement this helper function to 
     * displace the complete block information in user's terminal
     * It should be called when 'type=helper' is assigned!
     */
    virtual void PrintHelper()=0;

};


