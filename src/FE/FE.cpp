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
//+++ Date   : 2022.06.12
//+++ Purpose: this class offers the functions and management of
//+++          shape function class and qpoint class for FEM calc
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/FE.h"

FE::FE(){
    m_MaxDim=0;
    m_MinDim=0;
}

void FE::initdefault(const FECell &t_FECell){
    if(t_FECell.getFECellMaxDim()==1){
        // for 1d mesh
        m_MaxDim=1;
        m_MinDim=0;
        // for shape functions
        m_BulkShp.setMeshType(t_FECell.getFECellBulkElmtMeshType());
        m_BulkShp.init();

        // for gauss points
        m_BulkQpoints.setDim(t_FECell.getFECellMaxDim());
        m_BulkQpoints.setMeshType(t_FECell.getFECellBulkElmtMeshType());
        m_BulkQpoints.setQPointType(QPointType::GAUSSLEGENDRE);
        m_BulkQpoints.setOrder(t_FECell.getFECellBulkMeshOrder()+1);
        m_BulkQpoints.createQPoints();
    }
    else if(t_FECell.getFECellMaxDim()==2){
        // for 2d mesh
        m_MaxDim=2;
        m_MinDim=1;
        // for shape functions
        m_BulkShp.setMeshType(t_FECell.getFECellBulkElmtMeshType());
        m_BulkShp.init();

        m_LineShp.setMeshType(t_FECell.getFECellLineElmtMeshType());
        m_LineShp.init();

        // for gauss points
        m_BulkQpoints.setDim(t_FECell.getFECellMaxDim());
        m_BulkQpoints.setMeshType(t_FECell.getFECellBulkElmtMeshType());
        m_BulkQpoints.setQPointType(QPointType::GAUSSLEGENDRE);
        m_BulkQpoints.setOrder(t_FECell.getFECellBulkMeshOrder()+1);
        m_BulkQpoints.createQPoints();
        //
        m_LineQpoints.setDim(m_MinDim);
        m_LineQpoints.setMeshType(t_FECell.getFECellBulkElmtMeshType());
        m_LineQpoints.setQPointType(QPointType::GAUSSLEGENDRE);
        m_LineQpoints.setOrder(t_FECell.getFECellBulkMeshOrder()+1);
        m_LineQpoints.createQPoints();
    }
    else if(t_FECell.getFECellMaxDim()==3){
        // for 3d mesh
        m_MaxDim=3;
        m_MinDim=t_FECell.getFeCellMinDim();
        // for shape functions
        m_BulkShp.setMeshType(t_FECell.getFECellBulkElmtMeshType());
        m_BulkShp.init();

        m_SurfaceShp.setMeshType(t_FECell.getFECellSurfElmtMeshType());
        m_SurfaceShp.init();

        m_LineShp.setMeshType(t_FECell.getFECellLineElmtMeshType());
        m_LineShp.init();

        // for gauss points
        m_BulkQpoints.setDim(t_FECell.getFECellMaxDim());
        m_BulkQpoints.setMeshType(t_FECell.getFECellBulkElmtMeshType());
        m_BulkQpoints.setQPointType(QPointType::GAUSSLEGENDRE);
        m_BulkQpoints.setOrder(t_FECell.getFECellBulkMeshOrder()+1);
        m_BulkQpoints.createQPoints();
        //
        m_SurfaceQpoints.setDim(2);
        m_SurfaceQpoints.setMeshType(t_FECell.getFECellSurfElmtMeshType());
        m_SurfaceQpoints.setQPointType(QPointType::GAUSSLEGENDRE);
        m_SurfaceQpoints.setOrder(t_FECell.getFECellBulkMeshOrder()+1);
        m_SurfaceQpoints.createQPoints();
        //
        m_LineQpoints.setDim(1);
        m_LineQpoints.setMeshType(t_FECell.getFECellLineElmtMeshType());
        m_LineQpoints.setQPointType(QPointType::GAUSSLEGENDRE);
        m_LineQpoints.setOrder(t_FECell.getFECellBulkMeshOrder()+1);
        m_LineQpoints.createQPoints();
    }
    else{
        MessagePrinter::printErrorTxt("unsupported dim(="+to_string(t_FECell.getFECellMaxDim())+") for FE initializing");
        MessagePrinter::exitAsFem();
    }
}
void FE::init(const FECell &t_FECell){
    // this should only be used on user defines shp and qp block in their input file
    // here the different shp and qpoint have already been defined/given in your input file !!!
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    m_MaxDim=t_FECell.getFECellMaxDim();
    m_MinDim=t_FECell.getFeCellMinDim();
    if(t_FECell.getFECellMaxDim()==1){
        // for 1d case
        m_BulkShp.init();
        m_BulkQpoints.createQPoints();
    }
    else if(t_FECell.getFECellMaxDim()==2){
        // for 2d case
        m_BulkShp.init();
        m_BulkQpoints.createQPoints();

        m_LineShp.init();
        m_LineQpoints.createQPoints();
    }
    else if(t_FECell.getFECellMaxDim()==3){
        // for 3d case
        m_BulkShp.init();
        m_BulkQpoints.createQPoints();

        m_SurfaceShp.init();
        m_SurfaceQpoints.createQPoints();

        m_LineShp.init();
        m_LineQpoints.createQPoints();
    }
    else{
        MessagePrinter::printErrorTxt("unsupported dim(="+to_string(t_FECell.getFECellMaxDim())+") for FE initializing");
        MessagePrinter::exitAsFem();
    }
}

void FE::printFEInfo()const{
    if(getMaxDim()==1){
        MessagePrinter::printNormalTxt("Qpoint info for bulk elmts");
        m_BulkQpoints.printQPointsInfo();
    }
    else if(getMaxDim()==2){
        MessagePrinter::printNormalTxt("Qpoint info for bulk elmts");
        m_BulkQpoints.printQPointsInfo();

        MessagePrinter::printNormalTxt("Qpoint info for line elmts");
        m_LineQpoints.printQPointsInfo();
    }
    else if(getMaxDim()==3){
        MessagePrinter::printNormalTxt("Qpoint info for bulk elmts");
        m_BulkQpoints.printQPointsInfo();

        MessagePrinter::printNormalTxt("Qpoint info for surface elmts");
        m_SurfaceQpoints.printQPointsInfo();

        if(getMinDim()==1){
            MessagePrinter::printNormalTxt("Qpoint info for line elmts");
            m_LineQpoints.printQPointsInfo();
        }
    }
    MessagePrinter::printStars();
}

void FE::releaseMemory(){
    m_BulkShp.releaseMemory();
    m_LineShp.releaseMemory();
    m_SurfaceShp.releaseMemory();

    m_BulkQpoints.releaseMemory();
    m_LineQpoints.releaseMemory();
    m_SurfaceQpoints.releaseMemory();
}