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
//+++ Date   : 2024.07.10
//+++ Purpose: this class stores the coordinates of nodes of 
//+++          a single element
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FECell/Nodes.h"

Nodes::Nodes(){
    m_coordinates.clear();m_size=0;
}
Nodes::Nodes(const int &n){
    m_size=n;m_coordinates.resize(n*3,0.0);
}
Nodes::Nodes(const Nodes &nodes){
    m_size=nodes.m_size;
    m_coordinates.clear();
    for(const auto &it:nodes.m_coordinates) m_coordinates.push_back(it);
}
   
Nodes::~Nodes(){
    m_size=0;
    m_coordinates.clear();
}