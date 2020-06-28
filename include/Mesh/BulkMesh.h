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
//+++ Date   : 2020.06.28
//+++ Purpose: Define bulk mesh class for the bulk domain(except 
//+++          the interface domain)
//+++          This class offers the functions to get as well as
//+++          create the mesh and the related geometric info
//+++          In order to import mesh from external file, one should
//+++          call the meshio class
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include<iostream>
#include<iomanip>
#include<vector>
#include<cmath>
#include<set>

#include "Mesh/MeshType.h"


using namespace std;


class BulkMesh{
public:
    BulkMesh();

    void CreateMesh();
    void SaveMesh();
    //************************************************************
    //*** for the basic settings
    //************************************************************
    void SetDim(const int &ndim) {_nMaxDim=ndim;_nMinDim=ndim-1;}
    void SetMaxDim(const int &ndim) {_nMaxDim=ndim;}
    void SetMinDim(const int &ndim) {_nMinDim=ndim;}
    //*** for element numbers along different axis
    void SetNx(const int &nx){_Nx=nx;}
    void SetNy(const int &ny){_Ny=ny;}
    void SetNz(const int &nz){_Nz=nz;}
    //*** for element geometry settings
    void SetXmin(const double &xmin){_Xmin=xmin;}
    void SetXmax(const double &xmax){_Xmax=xmax;}
    void SetYmin(const double &ymin){_Ymin=ymin;}
    void SetYmax(const double &ymax){_Ymax=ymax;}
    void SetZmin(const double &zmin){_Zmin=zmin;}
    void SetZmax(const double &zmax){_Zmax=zmax;}
    //*** for mesh type setting
    void SetMeshType(MeshType type){_BulkMeshType=type;}

    //************************************************************
    //*** for the basic getting functions
    //************************************************************
    //*** for dimention of the bulk mesh
    inline int GetBulkMeshDim()const{return _nMaxDim;}
    inline int GetBulkMeshMinDim()const{return _nMinDim;}
    //*** for nodes num of the bulk mesh
    inline int GetBulkMeshNodesNum()const{return _nNodes;}
    //*** for elmts num of the bulk mesh
    inline int GetNx()const{return _Nx;}
    inline int GetNy()const{return _Ny;}
    inline int GetNz()const{return _Nz;}
    inline int GetBulkMeshElmtsNum()const{return _nElmts;}
    inline int GetBulkMeshBulkElmtsNum()const{return _nBulkElmts;}
    inline int GetBulkMeshSurfaceElmtsNum()const{return _nSurfaceElmts;}
    inline int GetBulkMeshLineElmtsNum()const{return _nLineElmts;}
    //*** for mesh order
    inline int GetBulkMeshOrder()const{return _nOrder;}
    //*** for mesh type
    inline MeshType GetBulkMeshBulkElmtType()const{return _BulkMeshType;}
    inline MeshType GetBulkMeshSurfaceElmtType()const{return _SurfaceMeshType;}
    inline MeshType GetBulkMeshLineElmtType()const{return _LineMeshType;}
    //*** for node's coordinate
    inline double GetIthNodeJthCoord(const int &i,const int &j)const{return _NodeCoords[(i-1)*3+j-1];}
    inline void   GetIthNodeCoords(const int &i,double (&coords)[4])const{
        coords[0]=0.0;
        coords[1]=_NodeCoords[(i-1)*3+1-1];
        coords[2]=_NodeCoords[(i-1)*3+2-1];
        coords[3]=_NodeCoords[(i-1)*3+3-1];
    }
    inline void   GetIthNodeCoords(const int &i,vector<double> &coords)const{
        coords[0]=0.0;
        coords[1]=_NodeCoords[(i-1)*3+1-1];
        coords[2]=_NodeCoords[(i-1)*3+2-1];
        coords[3]=_NodeCoords[(i-1)*3+3-1];
    }
    inline vector<double> GetIthNodeCoords(const int &i)const{
        vector<double> coords(3,0.0);
        coords[1]=_NodeCoords[(i-1)*3+1-1];
        coords[2]=_NodeCoords[(i-1)*3+2-1];
        coords[3]=_NodeCoords[(i-1)*3+3-1];
        return coords;
    }
    //*** for elmt connectivity information
    inline int  GetIthElmtNodesNum(const int &i)const{return _ElmtConn[i-1][0];}
    inline int  GetIthElmtJthNodeID(const int &i,const int &j)const{return _ElmtConn[i-1][j];}
    inline void GetIthElmtNodeIDs(const int &i,vector<int> &elConn)const{
        for(int j=1;j<=_ElmtConn[i-1][0];j++) elConn[j-1]=_ElmtConn[i-1][j];
    }
    inline vector<int> GetIthElmtNodeIDs(const int &i)const{
        vector<int> temp(_ElmtConn[i-1][0],0);
        for(int j=1;j<=_ElmtConn[i-1][0];j++) temp[j-1]=_ElmtConn[i-1][j];
        return temp;
    }
    //*** for physical group information
    inline int GetPhysicalGroupNum()const{return _nPhysicalGroups;}
    inline int GetIthPhysicalID(const int &i)const{return _PhysicalGroupIDList[i-1];}
    inline int GetIthPhysicalDim(const int &i)const{return _PhysicalGroupDimList[i-1];}
    inline string GetIthPhysicalName(const int &i)const{return _PhysicalGroupNameList[i-1];}
    inline string GetPhysicalNameViaID(const int &phyid)const{
        for(const auto &it:_PhysicalGroupID2NameList){
            if(it.first==phyid){
                return it.second;
            }
        }
        return "";
    }
    inline int GetPhysicalIDViaName(string phyname)const{
        for(const auto &it:_PhysicalGroupName2IDList){
            if(it.first==phyname){
                return it.second;
            }
        }
        return -1;
    }
    inline int GetNodesNumPerElmtViaPhysicalName(string phyname)const{
        for(const auto &it:_PhysicalGroupName2NodesNumPerElmtList){
            if(it.first==phyname){
                return it.second;
            }
        }
        return -1;
    }


    //************************************************************
    //*** for mesh information printer
    //************************************************************
    void PrintBulkMesh()const;
    void PrintBulkMeshDetails()const;

protected:
    //************************************************************
    //*** for the basic information of our mesh
    //************************************************************
    int _nNodes,_nElmts;
    int _nBulkElmts,_nSurfaceElmts,_nLineElmts;
    int _nMaxDim,_nMinDim;
    int _Nx,_Ny,_Nz;
    double _Xmax,_Xmin,_Ymax,_Ymin,_Zmin,_Zmax;
    int _nOrder;
    MeshType _BulkMeshType,_SurfaceMeshType,_LineMeshType;
    vector<double>      _NodeCoords;// store all the nodes/controlpts' coordinate
    vector<vector<int>> _ElmtConn;  // store all the element(bulk+surface+line+node elements)
    vector<double>      _ElmtVolume;// it could be: volume(3D), area(2D) or length(1D)


    //************************************************************
    //*** for the basic physical group information
    //************************************************************
    int _nPhysicalGroups;
    vector<string>                   _PhysicalGroupNameList;
    vector<int>                      _PhysicalGroupIDList;
    vector<int>                      _PhysicalGroupDimList;
    vector<pair<int,string>>         _PhysicalGroupID2NameList;
    vector<pair<string,int>>         _PhysicalGroupName2IDList;
    vector<pair<string,int>>         _PhysicalGroupName2NodesNumPerElmtList;
    vector<pair<string,vector<int>>> _PhysicalName2ElmtIDsList;
    vector<pair<string,vector<int>>> _PhysicalName2NodeIDsList;

};