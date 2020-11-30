//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.07.12
//+++ Purpose: implement the general FE space for FEM calculation
//+++          in AsFem, here one can use:
//+++            1) gauss integration 
//+++            2) shape functions for different mesh
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/FE.h"

FE::FE(){
    _nDim=1;_nMinDim=0;
    _HasDimSet=false;
    _IsInit=false;
    _BulkQPoint.SetQPointType(QPointType::GAUSSLEGENDRE);
    _SurfaceQPoint.SetQPointType(QPointType::GAUSSLEGENDRE);
    _LineQPoint.SetQPointType(QPointType::GAUSSLEGENDRE);
    _nBulkQpOrder=1;_nBCQpOrder=1;
}
//**********************************************
void FE::SetQPointType(QPointType qptype){
    _LineQPoint.SetQPointType(qptype);
    _SurfaceQPoint.SetQPointType(qptype);
    _BulkQPoint.SetQPointType(qptype);
}
void FE::SetBulkQpOrder(int order){
    _nBulkQpOrder=order;
    _LineQPoint.SetQPointOrder(order);
    _SurfaceQPoint.SetQPointOrder(order);
    _BulkQPoint.SetQPointOrder(order);
}
void FE::SetBCQpOrder(int order){
    _nBCQpOrder=order;
    _LineQPoint.SetQPointOrder(order);
    _SurfaceQPoint.SetQPointOrder(order);
}
//******************************************************
void FE::CreateQPoints(Mesh &mesh){
    if(_HasDimSet){
        if(GetDim()==1){
            _BulkQPoint.SetDim(1);
            _BulkQPoint.CreateQpoints(mesh.GetBulkMeshBulkElmtType());
        }
        else if(GetDim()==2){
            _BulkQPoint.SetDim(2);
            _BulkQPoint.CreateQpoints(mesh.GetBulkMeshBulkElmtType());

            _LineQPoint.SetDim(1);
            _LineQPoint.CreateQpoints(mesh.GetBulkMeshLineElmtType());
        }
        else if(GetDim()==3){
            _BulkQPoint.SetDim(3);
            _BulkQPoint.CreateQpoints(mesh.GetBulkMeshBulkElmtType());

            _SurfaceQPoint.SetDim(2);
            _SurfaceQPoint.CreateQpoints(mesh.GetBulkMeshSurfaceElmtType());

            _LineQPoint.SetDim(1);
            _LineQPoint.CreateQpoints(mesh.GetBulkMeshLineElmtType());
        }
    }
    else{
        MessagePrinter::PrintErrorTxt("can\'t create qpoints for FE space, the dim has not been given yet");
        MessagePrinter::AsFem_Exit();
    }
}
//**************************************************************************
//*** for shape function related functions
//**************************************************************************
void FE::CreateShapeFuns(Mesh &mesh){
    SetDim(mesh.GetBulkMeshDim());
    SetMinDim(mesh.GetBulkMeshMinDim());

    _BulkShp=ShapeFun(mesh.GetBulkMeshDim(),mesh.GetBulkMeshBulkElmtType());
    _BulkShp.PreCalc();

    _BulkNodes=Nodes(mesh.GetBulkMeshNodesNumPerBulkElmt());
    if(GetDim()==3){
        _SurfaceShp=ShapeFun(2,mesh.GetBulkMeshSurfaceElmtType());
        _SurfaceShp.PreCalc();

        _LineShp=ShapeFun(1,mesh.GetBulkMeshLineElmtType());
        _LineShp.PreCalc();

        _SurfaceNodes=Nodes(mesh.GetBulkMeshNodesNumPerSurfaceElmt());
        _LineNodes=Nodes(mesh.GetBulkMeshNodesNumPerLineElmt());
    }
    else if(GetDim()==2){
        _LineShp=ShapeFun(1,mesh.GetBulkMeshLineElmtType());
        _LineShp.PreCalc();

        _LineNodes=Nodes(mesh.GetBulkMeshNodesNumPerLineElmt());
    }
}
//**************************************************************************
void FE::InitFE(Mesh &mesh){
    CreateQPoints(mesh);
    CreateShapeFuns(mesh);
}
//*******************************************
void FE::PrintFEInfo()const{
    string msg;
    MessagePrinter::PrintNormalTxt("Summary information of qpoint");
    if(_BulkQPoint.GetQpPointType()==QPointType::GAUSSLEGENDRE){
        msg="  qpoint type= Gauss-Legendre, dim="+to_string(GetDim());
    }
    else if(_BulkQPoint.GetQpPointType()==QPointType::GAUSSLOBATTO){
        msg="  qpoint type= Gauss-Lobatto, dim="+to_string(GetDim());
    }
    MessagePrinter::PrintNormalTxt(msg);

    if(GetDim()==1){
        msg="  for bulk: order="+to_string(_BulkQPoint.GetQpOrder())
           +", num of qpoints="+to_string(_BulkQPoint.GetQpPointsNum());
        MessagePrinter::PrintNormalTxt(msg);
    }
    else if(GetDim()==2){
        msg="  for boun: order="+to_string(_BulkQPoint.GetQpOrder())
           +", num of qpoints="+to_string(_BulkQPoint.GetQpPointsNum());
        MessagePrinter::PrintNormalTxt(msg);

        msg="  for boun: order="+to_string(_LineQPoint.GetQpOrder())
           +", num of qpoints="+to_string(_LineQPoint.GetQpPointsNum());
        MessagePrinter::PrintNormalTxt(msg);
    }
    else if(GetDim()==3){
        msg="  for bulk: order="+to_string(_BulkQPoint.GetQpOrder())
           +", num of qpoints="+to_string(_BulkQPoint.GetQpPointsNum());
        MessagePrinter::PrintNormalTxt(msg);

        msg="  for boun: order="+to_string(_SurfaceQPoint.GetQpOrder())
           +", num of qpoints="+to_string(_SurfaceQPoint.GetQpPointsNum());
        MessagePrinter::PrintNormalTxt(msg);
    }
    MessagePrinter::PrintDashLine();
}