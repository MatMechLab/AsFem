//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author  : Yang Bai
//+++ Date    : 2021.02.18
//+++ Reviewer: Xiaoyuan @ 2021.08.28
//+++ Purpose : implement Abaqus' inp mesh file importer
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <algorithm>
#include "Mesh/MeshIOBase.h"
#include "Mesh/MeshType.h"

#include "Utils/MessagePrinter.h"
#include "Utils/StringUtils.h"

class AbaqusIO:public MeshIOBase{
public:
    AbaqusIO(){
        _MeshFileName.clear();
        _HasSetMeshFileName=false;
    }

    virtual bool ReadMeshFromFile(Mesh &mesh) override;
    virtual void SetMeshFileName(string filename)override{_MeshFileName=filename;_HasSetMeshFileName=true;}
    virtual string GetMeshFileName()const override{return _MeshFileName;}

private:
    int GetElmtNodesNumFromInp() const;
    int GetElmtDimFromInp() const;
    int GetElmtOrderFromInp()const;
    int GetElmtVTKCellTypeFromInp()const;
    MeshType GetElmtMeshTypeFromInp()const;
    string GetElmtMeshTypeNameFromInp()const;

    int GetSubElmtNodesNumFromInp() const;
    int GetSubElmtDimFromInp() const;
    int GetSubElmtOrderFromInpElmt()const;
    MeshType GetSubElmtMeshTypeFromInp()const;

    int GetSubSubElmtNodesNumFromInp() const;
    MeshType GetSubSubElmtMeshTypeFromInp()const;


    int GetNodesNumFromInp()const;
    int GetElmtsNumFromInp()const;
    int GetNsetsNumFromInp()const;
    int GetElsetsNumFromInp()const;
    int GetSurfacesNumFromInp()const;
    int GetSurfaceElmtsNumFromInp()const;
    int GetSurfaceEdgeIDViaSurfaceNameFromInp(string surfacesetname)const;
    vector<int> GetSurfaceElmtIDViaSurfaceNameFromInp(string surfacesetname)const;


    ifstream _in;

    int _nMaxDim=-1,_nMinDim=4;
    int _nPhysicGroups=0;
    int _nNodeSetPhysicalGroups=0;
    int _nNodes=0,_nElmts=0;
    int _nBulkElmts=0,_nSurfaceElmts=0,_nLineElmts=0;

    double _Xmax,_Xmin,_Ymax,_Ymin,_Zmax,_Zmin;

    int _nNodesPerBulkElmt=-1;
    int _nNodesPerLineElmt=0;
    int _nNodesPerSurfaceElmt=0;
    int _nOrder=1;

};
