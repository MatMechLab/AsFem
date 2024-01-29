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

    /**
     * setup mesh geometry info for 1d case
     * @param nx mesh number along x-axis
     * @param xmin xmin value of the domain
     * @param xmax xmax value of the domain
     * @param meshtype the given mesh type
    */
    void setMeshInfo(const int &nx,const double &xmin,const double &xmax,const MeshType &meshtype);
    /**
     * setup mesh geometry info for 2d case
     * @param nx mesh number along x-axis
     * @param ny mesh number along y-axis
     * @param xmin xmin value of the domain
     * @param xmax xmax value of the domain
     * @param ymin ymin value of the domain
     * @param ymax ymax value of the domain
     * @param meshtype the given mesh type
    */
    void setMeshInfo(const int &nx,const int &ny,
                     const double &xmin,const double &xmax,
                     const double &ymin,const double &ymax,const MeshType &meshtype);
    /**
     * setup mesh geometry info for 3d case
     * @param nx mesh number along x-axis
     * @param ny mesh number along y-axis
     * @param xmin xmin value of the domain
     * @param xmax xmax value of the domain
     * @param ymin ymin value of the domain
     * @param ymax ymax value of the domain
     * @param meshtype the given mesh type
    */
    void setMeshInfo(const int &nx,const int &ny,const int &nz,
                     const double &xmin,const double &xmax,
                     const double &ymin,const double &ymax,
                     const double &zmin,const double &zmax,const MeshType &meshtype);


    /**
     * Get the reference of fecell data
    */
    inline FECellData& getCellDataRef(){return m_CellData;}

    /**
     * Save the FECell mesh to vtu file
    */
    void saveFECell2VTUFile()const;
    /**
     * Print out the summary information of FECell class, i.e., mpi-based mesh distribution, elements num, nodes num, etc.
    */
    void printSummaryInfo()const;

private:
    FECellData m_CellData;/**< the FE cell data structure of the whole system */


};