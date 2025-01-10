//****************************************************************
//* This file is part of the AsFem framework
//* Advanced Simulation kit based on Finite Element Method (AsFem)
//* All rights reserved, Yang Bai/MM-Lab@CopyRight 2020-present
//* https://github.com/MatMechLab/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2022.06.05
//+++ Purpose: this class offers the functions and management of
//+++          shape function class and qpoint class for FEM calc
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "FE/QPoint.h"
#include "FE/ShapeFun.h"

#include "FECell/FECell.h"

/**
 * This class implement the general function and management of both shape functions and qpoints classes
 */
class FE{
public:
    /**
     * constructor
     */
    FE();
    /**
     * initialize the FE class, this option assumes that no info is given in your input file
     * @param t_FECell the fe cell class for FE initializing
     */
    void initdefault(const FECell &t_FECell);

    /**
     * init the FE class with preset info(defined in your input file)
     * @param t_FECell the fe cell class for FE initializing
     */
    void init(const FECell &t_FECell);

    /**
     * get the maximum dim of FE space
     */
    inline int getMaxDim()const{return m_MaxDim;}
    /**
     * get the mini dim of FE space
     */
    inline int getMinDim()const{return m_MinDim;}

    /**
     * print out the summary info of FE space
     */
    void printFEInfo()const;
    /**
     * release allocated memory
     */
    void releaseMemory();

public:
    ShapeFun m_BulkShp;/**< shapefun for the bulk element */
    ShapeFun m_SurfaceShp;/**< shapefun for the surface element (for 3d case) */
    ShapeFun m_LineShp;/**< shapefun for the line element (for 2d case) */

    QPoint m_BulkQpoints;/**< gauss integration points for the bulk element */
    QPoint m_SurfaceQpoints;/**< gauss integration points for the surface element */
    QPoint m_LineQpoints;/**< gauss integration points for the line element */

private:
    int m_MaxDim;/**< the max dimension of fe space */
    int m_MinDim;/**< the min dimension of fe space */

};