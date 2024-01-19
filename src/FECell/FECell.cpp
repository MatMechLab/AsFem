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
//+++ Date    : 2024.01.19
//+++ Function: finite element cell management class
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FECell/FECell.h"

FECell::FECell(){}

// setup 1d mesh info
void FECell::setMeshInfo(const int &nx,const double &xmin,const double &xmax,const MeshType &meshtype){
    m_CellData.Nx=nx;
    m_CellData.Xmin=xmin;
    m_CellData.Xmax=xmax;
    m_CellData.BulkElmtMeshType=meshtype;
}
// setup 2d mesh info
void FECell::setMeshInfo(const int &nx,const int &ny,
                         const double &xmin,const double &xmax,
                         const double &ymin,const double &ymax,const MeshType &meshtype){
    m_CellData.Nx=nx;
    m_CellData.Ny=ny;
    m_CellData.Xmin=xmin;
    m_CellData.Xmax=xmax;
    m_CellData.Ymin=ymin;
    m_CellData.Ymax=ymax;
    m_CellData.BulkElmtMeshType=meshtype;
}
// setup 3d mesh info
void FECell::setMeshInfo(const int &nx,const int &ny,const int &nz,
                         const double &xmin,const double &xmax,
                         const double &ymin,const double &ymax,
                         const double &zmin,const double &zmax,const MeshType &meshtype){
    m_CellData.Nx=nx;
    m_CellData.Ny=ny;
    m_CellData.Nz=nz;
    m_CellData.Xmin=xmin;
    m_CellData.Xmax=xmax;
    m_CellData.Ymin=ymin;
    m_CellData.Ymax=ymax;
    m_CellData.Zmin=zmin;
    m_CellData.Zmax=zmax;
    m_CellData.BulkElmtMeshType=meshtype;
}
void FECell::printSummaryInfo()const{
    
}