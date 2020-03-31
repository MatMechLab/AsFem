//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "FE/FE.h"

void FE::InitFE(Mesh &mesh){

    // SetOrder(mesh.GetMeshOrder());
    // SetBCOrder(mesh.GetMeshOrder());
    SetDim(mesh.GetMaxDim());
    SetDimMin(mesh.GetMinDim());

    _shp_bulk=ShapeFun(mesh.GetDim(),mesh.GetBulkMeshType());
    _shp_bulk.PreCalc();

    _qp_bulk=QPoint(mesh.GetDim(),GetOrder());
    _qp_bulk.SetQPointType(_QPointType);
    _qp_bulk.CreateQPoints(mesh.GetBulkMeshType());
    

    _nodes=Nodes(mesh.GetNodesNumPerBulkElmt());
    if(mesh.GetDim()==3){
        _shp_surface=ShapeFun(2,mesh.GetSurfaceMeshType());
        _shp_surface.PreCalc();

        _shp_line=ShapeFun(1,mesh.GetLineMeshType());
        _shp_line.PreCalc();

        _qp_surface=QPoint(2,GetBCOrder());
        _qp_surface.SetQPointType(_QPointType);
        _qp_surface.CreateQPoints(mesh.GetSurfaceMeshType());

        _qp_line=QPoint(1,GetBCOrder());
        _qp_line.SetQPointType(_QPointType);
        _qp_line.CreateQPoints(mesh.GetLineMeshType());

        _surface_nodes=Nodes(mesh.GetNodesNumPerSurfaceElmt());
        _line_nodes=Nodes(mesh.GetNodesNumPerLineElmt());
    }
    else if(mesh.GetDim()==2){
        _shp_line=ShapeFun(1,mesh.GetLineMeshType());
        _shp_line.PreCalc();

        _qp_line=QPoint(1,GetBCOrder());
        _qp_line.SetQPointType(_QPointType);
        _qp_line.CreateQPoints(mesh.GetLineMeshType());

        _line_nodes=Nodes(mesh.GetNodesNumPerLineElmt());
    }
}