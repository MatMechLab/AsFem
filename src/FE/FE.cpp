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
    m_maxdim=0;
    m_mindim=0;
}

void FE::initdefault(const FECell &t_fecell){
    if(t_fecell.getFECellMaxDim()==1){
        // for 1d mesh
        m_maxdim=1;
        m_mindim=0;
        // for shape functions
        m_bulk_shp.setMeshType(t_fecell.getFECellBulkElmtMeshType());
        m_bulk_shp.init();

        // for gauss points
        m_bulk_qpoints.setDim(t_fecell.getFECellMaxDim());
        m_bulk_qpoints.setMeshType(t_fecell.getFECellBulkElmtMeshType());
        m_bulk_qpoints.setQPointType(QPointType::GAUSSLEGENDRE);
        m_bulk_qpoints.setOrder(t_fecell.getFECellBulkMeshOrder()+1);
        m_bulk_qpoints.createQPoints();
    }
    else if(t_fecell.getFECellMaxDim()==2){
        // for 2d mesh
        m_maxdim=2;
        m_mindim=1;
        // for shape functions
        m_bulk_shp.setMeshType(t_fecell.getFECellBulkElmtMeshType());
        m_bulk_shp.init();

        m_line_shp.setMeshType(t_fecell.getFECellLineElmtMeshType());
        m_line_shp.init();

        // for gauss points
        m_bulk_qpoints.setDim(t_fecell.getFECellMaxDim());
        m_bulk_qpoints.setMeshType(t_fecell.getFECellBulkElmtMeshType());
        m_bulk_qpoints.setQPointType(QPointType::GAUSSLEGENDRE);
        m_bulk_qpoints.setOrder(t_fecell.getFECellBulkMeshOrder()+1);
        m_bulk_qpoints.createQPoints();
        //
        m_line_qpoints.setDim(m_mindim);
        m_line_qpoints.setMeshType(t_fecell.getFECellBulkElmtMeshType());
        m_line_qpoints.setQPointType(QPointType::GAUSSLEGENDRE);
        m_line_qpoints.setOrder(t_fecell.getFECellBulkMeshOrder()+1);
        m_line_qpoints.createQPoints();
    }
    else if(t_fecell.getFECellMaxDim()==3){
        // for 3d mesh
        m_maxdim=3;
        m_mindim=t_fecell.getFeCellMinDim();
        // for shape functions
        m_bulk_shp.setMeshType(t_fecell.getFECellBulkElmtMeshType());
        m_bulk_shp.init();

        m_surface_shp.setMeshType(t_fecell.getFECellSurfElmtMeshType());
        m_surface_shp.init();

        m_line_shp.setMeshType(t_fecell.getFECellLineElmtMeshType());
        m_line_shp.init();

        // for gauss points
        m_bulk_qpoints.setDim(t_fecell.getFECellMaxDim());
        m_bulk_qpoints.setMeshType(t_fecell.getFECellBulkElmtMeshType());
        m_bulk_qpoints.setQPointType(QPointType::GAUSSLEGENDRE);
        m_bulk_qpoints.setOrder(t_fecell.getFECellBulkMeshOrder()+1);
        m_bulk_qpoints.createQPoints();
        //
        m_surface_qpoints.setDim(2);
        m_surface_qpoints.setMeshType(t_fecell.getFECellSurfElmtMeshType());
        m_surface_qpoints.setQPointType(QPointType::GAUSSLEGENDRE);
        m_surface_qpoints.setOrder(t_fecell.getFECellBulkMeshOrder()+1);
        m_surface_qpoints.createQPoints();
        //
        m_line_qpoints.setDim(1);
        m_line_qpoints.setMeshType(t_fecell.getFECellLineElmtMeshType());
        m_line_qpoints.setQPointType(QPointType::GAUSSLEGENDRE);
        m_line_qpoints.setOrder(t_fecell.getFECellBulkMeshOrder()+1);
        m_line_qpoints.createQPoints();
    }
    else{
        MessagePrinter::printErrorTxt("unsupported dim(="+to_string(t_fecell.getFECellMaxDim())+") for FE initializing");
        MessagePrinter::exitAsFem();
    }
}
void FE::init(const FECell &t_fecell){
    // this should only be used on user defines shp and qp block in their input file
    // here the different shp and qpoint have already been defined/given in your input file !!!
    m_maxdim=t_fecell.getFECellMaxDim();
    m_mindim=t_fecell.getFeCellMinDim();
    if(t_fecell.getFECellMaxDim()==1){
        // for 1d case
        m_bulk_shp.init();
        m_bulk_qpoints.createQPoints();
    }
    else if(t_fecell.getFECellMaxDim()==2){
        // for 2d case
        m_bulk_shp.init();
        m_bulk_qpoints.createQPoints();

        m_line_shp.init();
        m_line_qpoints.createQPoints();
    }
    else if(t_fecell.getFECellMaxDim()==3){
        // for 3d case
        m_bulk_shp.init();
        m_bulk_qpoints.createQPoints();

        m_surface_shp.init();
        m_surface_qpoints.createQPoints();

        m_line_shp.init();
        m_line_qpoints.createQPoints();
    }
    else{
        MessagePrinter::printErrorTxt("unsupported dim(="+to_string(t_fecell.getFECellMaxDim())+") for FE initializing");
        MessagePrinter::exitAsFem();
    }
}

void FE::printFEInfo()const{
    if(getMaxDim()==1){
        MessagePrinter::printNormalTxt("Qpoint info for bulk elmts");
        m_bulk_qpoints.printQPointsInfo();
    }
    else if(getMaxDim()==2){
        MessagePrinter::printNormalTxt("Qpoint info for bulk elmts");
        m_bulk_qpoints.printQPointsInfo();

        MessagePrinter::printNormalTxt("Qpoint info for line elmts");
        m_line_qpoints.printQPointsInfo();
    }
    else if(getMaxDim()==3){
        MessagePrinter::printNormalTxt("Qpoint info for bulk elmts");
        m_bulk_qpoints.printQPointsInfo();

        MessagePrinter::printNormalTxt("Qpoint info for surface elmts");
        m_surface_qpoints.printQPointsInfo();

        if(getMinDim()==1){
            MessagePrinter::printNormalTxt("Qpoint info for line elmts");
            m_line_qpoints.printQPointsInfo();
        }
    }
    MessagePrinter::printStars();
}

void FE::releaseMemory(){
    m_bulk_shp.releaseMemory();
    m_line_shp.releaseMemory();
    m_surface_shp.releaseMemory();

    m_bulk_qpoints.releaseMemory();
    m_line_qpoints.releaseMemory();
    m_surface_qpoints.releaseMemory();
}