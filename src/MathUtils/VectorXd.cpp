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
//+++ Date   : 2020.10.18
//+++ Purpose: Define the general Vector array in AsFem.
//+++          We mainly use this for the calculation of residual.
//+++          If one wants to use Eigen's VectorXd, please use
//+++          Eigen::VectorXd, which is different with ours !!!
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "MathUtils/VectorXd.h"

VectorXd::VectorXd(){
    m_vals.clear();
    m_m=0;
}
VectorXd::VectorXd(const VectorXd &a){
    m_vals=a.m_vals;
    m_m=a.m_m;
}
VectorXd::VectorXd(const int &m){
    m_m=m;
    m_vals.resize(m,0.0);
}
VectorXd::VectorXd(const int &m,const double &val){
    m_m=m;
    m_vals.resize(m,val);
}