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
//+++ Date   : 2022.06.05
//+++ Purpose: this class offers the functions and management of
//+++          shape function class and qpoint class for FEM calc
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "FE/QPoint.h"
#include "FE/ShapeFun.h"

#include "Mesh/Mesh.h"

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
     * @param t_mesh the mesh class for FE initializing
     */
    void initdefault(const Mesh &t_mesh);

    /**
     * init the FE class with preset info(defined in your input file)
     * @param t_mesh the mesh class for FE initializing
     */
    void init(const Mesh &t_mesh);

    /**
     * get the maximum dim of FE space
     */
    inline int getMaxDim()const{return m_maxdim;}
    /**
     * get the mini dim of FE space
     */
    inline int getMinDim()const{return m_mindim;}

    /**
     * print out the summary info of FE space
     */
    void printFEInfo()const;
    /**
     * release allocated memory
     */
    void releaseMemory();

public:
    ShapeFun m_bulk_shp;/**< shapefun for the bulk element */
    ShapeFun m_surface_shp;/**< shapefun for the surface element (for 3d case) */
    ShapeFun m_line_shp;/**< shapefun for the line element (for 2d case) */

    QPoint m_bulk_qpoints;/**< gauss integration points for the bulk element */
    QPoint m_surface_qpoints;/**< gauss integration points for the surface element */
    QPoint m_line_qpoints;/**< gauss integration points for the line element */

private:
    int m_maxdim;/**< the max dimension of fe space */
    int m_mindim;/**< the min dimension of fe space */

};