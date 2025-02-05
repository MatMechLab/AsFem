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
//+++ Date   : 2022.08.22
//+++ Purpose: Defines the data structure for projection
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "MathUtils/Vector.h"

/**
 * This structure defines the basic data for projection
 */
struct ProjectionData{
    int m_ScalarProjMateNum;/**< number of scalar material to be projected */
    int m_VectorProjMateNum;/**< number of vector material to be projected */
    int m_Rank2ProjMateNum;/**< number of rank2tensor material to be projected */
    int m_Rank4ProjMateNum;/**< number of rank4tensor material to be projected */

    vector<Vector> m_ScalarProjMateVecList;/**< the list stores the projected scalar materials vector */
    vector<Vector> m_VectorProjMateVecList;/**< the list stores the projected vector materials vector */
    vector<Vector> m_Rank2ProjMateVecList;/**< the list stores the projected rank-2 materials vector */
    vector<Vector> m_Rank4ProjMateVecList;/**< the list stores the projected rank-4 materials vector */

    vector<string> m_ScalarProjMateNameList;/**< vector for the name of scalar materials to be projected */
    vector<string> m_VectorProjMateNamelist;/**< vector for the name of vector materials to be projected */
    vector<string> m_Rank2ProjMateNameList;/**< vector for the name of rank-2 tensor materials to be projected */
    vector<string> m_Rank4ProjMateNameList;/**< vector for the name of rank-4 tensor materials to be projected */
};