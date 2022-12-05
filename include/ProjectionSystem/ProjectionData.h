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
//+++ Date   : 2022.08.22
//+++ Purpose: Defines the data structure for projection
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "MathUtils/Vector.h"

/**
 * This structure defines the basic data for projection
 */
struct ProjectionData{
    int m_scalarmate_num;/**< number of scalar material to be projected */
    int m_vectormate_num;/**< number of vector material to be projected */
    int m_rank2mate_num;/**< number of rank2tensor material to be projected */
    int m_rank4mate_num;/**< number of rank4tensor material to be projected */
    Vector m_proj_scalarmate_vec;/**< the vector stores the projected scalar materials */
    Vector m_proj_vectormate_vec;/**< the vector stores the projected vector materials */
    Vector m_proj_rank2mate_vec;/**< the vector stores the projected rank-2 materials */
    Vector m_proj_rank4mate_vec;/**< the vector stores the projected rank-4 materials */

    vector<string> m_scalarmate_namelist;/**< vector for the name of scalar materials to be projected */
    vector<string> m_vectormate_namelist;/**< vector for the name of vector materials to be projected */
    vector<string> m_rank2mate_namelist;/**< vector for the name of rank-2 tensor materials to be projected */
    vector<string> m_rank4mate_namelist;/**< vector for the name of rank-4 tensor materials to be projected */
};