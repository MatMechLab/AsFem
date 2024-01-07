//****************************************************************
//* This file is part of the AsFem framework
//* Advanced Simulation kit based on Finite Element Method (AsFem)
//* All rights reserved, Yang Bai/MM-Lab@CopyRight 2020-present
//* https://github.com/MatMechLab/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author  : Yang Bai
//+++ Date    : 2023.12.30
//+++ Function: finite element cell management class
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once


#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>


#include "Utils/MessagePrinter.h"
#include "FECell/FECellData.h"


using std::vector;
using std::map;
using std::make_pair;
using std::string;

/**
 * This class defines and manages the basic/single finite element cell data structure,
 * users can use this class for either the elemental loop or system vector/matrix assemble.
*/
class FECell{
public:
    /**
     * Constructor
    */
    FECell();
    /**
     * Constructor
    */
    FECell(const FECell &a);

    void setMeshInfo(const int &nx,const double &xmin,const double &xmax,const MeshType &meshtype);


    /**
     * Print out the summary information of FECell class, i.e., mpi-based mesh distribution, elements num, nodes num, etc.
    */
    void printSummaryInfo()const;

private:
    FECellData CellData;/**< the FE cell data structure of the whole system */


};