//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_MESH_H
#define ASFEM_MESH_H


#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <algorithm>

#include "petsc.h"

// For AsFem's own header file
#include "MessagePrinter/MessagePrinter.h"
#include "MeshType.h"
#include "Nodes.h"

#include "ElmtSystem/ElmtSystem.h"

using namespace std;

class ElmtSystem;

class Mesh{
public:
    Mesh();
    void Reset();

    void SetElmtInfo(ElmtSystem &elmtSystem);

    //*****************************************
    //*** some common setting functions
    //*****************************************
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
    //*** for save mesh
    void SetMeshFileName(string name){_MeshFileName=name;}
    //*** set mesh type
    void SetMeshType(string meshtype);
    //***  for mesh mode(true for built-in, false for external one)
    void SetMeshMode(bool flag){
        _IsBuiltInMesh=flag;
        if(flag){
            _IsGmshMesh=false;
            _IsAbaqusMesh=false;
        }
    }
    void SetMeshGmshMode(bool flag){
        _IsGmshMesh=flag;
        if(flag) {
            _IsBuiltInMesh=false;
            _IsAbaqusMesh=false;
        }
    }
    void SetMeshAbaqusMode(bool flag){
        _IsAbaqusMesh=flag;
        if(flag){
            _IsBuiltInMesh=false;
            _IsGmshMesh=false;
        }
    }
    //*** set msh file name for gmsh import
    void SetMshFileName(string filename){_GmshFileName=filename;_AbaqusFileName=filename;}
    //*** set ith element's volume
    void SetIthElmtVolume(const int &i,const double &volume){
        _ElmtVolume[i-1]=volume;
    }
    //*** set i-th bulk elmt volume
    void SetIthBulkElmtVolume(const int &i,const double &volume){
        _ElmtVolume[i+_nElmts-_nBulkElmts-1]=volume;
    }
    //*** if one set the volume for each element, the you should update the total volume
    void SetTotalVolume(const double &volume){_TotalVolume=volume;}

    //*****************************************
    //*** some common getting functions
    //*****************************************
    inline int GetDim()const{return _nMaxDim;}
    inline int GetMaxDim()const{return _nMaxDim;}
    inline int GetMinDim()const{return _nMinDim;}
    //*** for element numbers
    inline int GetNx()const{return _Nx;}
    inline int GetNy()const{return _Ny;}
    inline int GetNz()const{return _Nz;}
    //*** for geometry
    inline double GetXmin()const{return _Xmin;}
    inline double GetXmax()const{return _Xmax;}
    inline double GetYmin()const{return _Ymin;}
    inline double GetYmax()const{return _Ymax;}
    inline double GetZmin()const{return _Zmin;}
    inline double GetZmax()const{return _Zmax;}
    //*** for mesh type
    inline MeshType GetBulkMeshType()const{return _BulkMeshType;}
    inline MeshType GetSurfaceMeshType()const{return _SurfaceMeshType;}
    inline MeshType GetLineMeshType()const{return _LineMeshType;}
    inline PetscInt GetMeshOrder()const{return _nOrder;}
    //*** for nodes number and element numbers
    inline int GetNodesNum()const{return _nNodes;}
    inline int GetNodesNumPerBulkElmt()const{return _nNodesPerBulkElmt;}
    inline int GetNodesNumPerLineElmt()const{return _nNodesPerLineElmt;}
    inline int GetNodesNumPerSurfaceElmt()const{return _nNodesPerSurfaceElmt;}
    inline int GetElmtsNum()const{return _nElmts;}
    inline int GetBulkElmtsNum()const{return _nBulkElmts;}
    inline int GetBCElmtsNum()const{return _nElmts-_nBulkElmts;}

    //*** for bc elements information
    inline int GetDimViaPhyName(string bcname)const{
        for(auto it:_PhysicNameToDimList){
            if(it.first==bcname){
                return it.second;
            }
        }
        return -1;
    }
    inline int GetBCElmtNodesNumViaPhyName(string bcname)const{
        for(auto it:_PhysicNameToElmtNodesNumList){
            if(it.first==bcname){
                return it.second;
            }
        }
        return -1;
    }

    //**** for physical groups
    inline int GetPhysicalGroupNum()const{return _nPhysicGroups;}
    //**** for i-th element information
    inline int GetIthBulkElmtVTKType(const int &i)const{
        return _ElmtVTKCellType[i+GetBCElmtsNum()-1];
    }
    inline ElmtType GetIthBulkElmtElmtType(const int &i)const{
        return _MeshElmtInfoList[i-1].first;
    }
    inline MateType GetIthBulkElmtMateType(const int &i)const{
        return _MeshElmtInfoList[i-1].second.first;
    }
    inline int GetIthBulkElmtMateIndex(const int &i)const{
        return _MeshElmtInfoList[i-1].second.second;
    }
    inline int GetIthElmtNodesNum(const int &i)const{
        return _ElmtConn[i-1][0];
    }
    inline int GetIthBulkElmtNodesNum(const int &i)const{
        return _ElmtConn[i-1+_nElmts-_nBulkElmts][0];
    }
    inline int GetElmtsNumViaPhyName(string phyname)const{
        for(auto it:_PhysicNameToElmtIndexSet){
            if(it.first==phyname){
                return PetscInt(it.second.size());
            }
        }
        return 0;
    }
    inline int GetIthElmtIndexViaPhyName(string phyname,const int &i)const{
        for(auto it:_PhysicNameToElmtIndexSet){
            if(it.first==phyname){
                return it.second[i-1];
            }
        }
        return 0;
    }
    inline int GetIthElmtNodesNumViaPhyName(string phyname,const int &i)const{
        for(auto it:_PhysicNameToElmtIndexSet){
            if(it.first==phyname){
                return GetIthElmtNodesNum(it.second[i-1]);
            }
        }
        return 0;
    }
    //*** for i-th element connectivity
    inline void GetIthBulkElmtConn(const int &i,vector<int> &elconn)const{
        for(int j=1;j<=_ElmtConn[i-1+GetBCElmtsNum()][0];++j){
            elconn[j-1]=_ElmtConn[i-1+GetBCElmtsNum()][j];
        }
    }
    vector<int> GetIthBulkElmtDofsIndex(const int &i)const{
        return _BulkElmtDofIndexList[i-1];
    }
    inline int GetIthElmtJthConn(const int &i,const int &j)const{
        return _ElmtConn[i-1][j];
    }
    void GetIthElmtNodes(const int &e,Nodes &nodes)const{
        for(int i=1;i<=GetIthElmtNodesNum(e);++i){
            nodes(i,0)=GetIthNodeJthCoord(GetIthElmtJthConn(e,i),0);
            nodes(i,1)=GetIthNodeJthCoord(GetIthElmtJthConn(e,i),1);
            nodes(i,2)=GetIthNodeJthCoord(GetIthElmtJthConn(e,i),2);
            nodes(i,3)=GetIthNodeJthCoord(GetIthElmtJthConn(e,i),3);
        }
    }
    void GetIthBulkElmtNodes(const int &e,Nodes &nodes)const{
        for(int i=1;i<=GetIthBulkElmtNodesNum(e);++i){
            nodes(i,0)=GetIthNodeJthCoord(GetIthBulkElmtJthConn(e,i),0);
            nodes(i,1)=GetIthNodeJthCoord(GetIthBulkElmtJthConn(e,i),1);
            nodes(i,2)=GetIthNodeJthCoord(GetIthBulkElmtJthConn(e,i),2);
            nodes(i,3)=GetIthNodeJthCoord(GetIthBulkElmtJthConn(e,i),3);
        }
    }
    inline int GetIthBulkElmtJthConn(const int &i,const int &j)const{
        return _ElmtConn[i-1+_nElmts-_nBulkElmts][j];
    }
    inline int GetIthElmtVTKCellType(const int &i)const{
        return _ElmtVTKCellType[i-1];
    }
    inline int GetIthBulkElmtVTKCellType(const int &i)const{
        return _ElmtVTKCellType[i-1+_nElmts-_nBulkElmts];
    }

    //*****************************************************
    //*** get information for nodes and connectivity
    //*****************************************************
    inline double GetIthNodeJthCoord(const int &i,const int &j)const{
        return _NodeCoords[(i-1)*4+j];
    }
    //*** for elemental 'volume'
    inline double GetTotalVolume()const{return _TotalVolume;}
    inline double GetIthElmtVolume(const int &i)const{
        return _ElmtVolume[i-1];
    }

public:
    //*****************************************
    //*** For mesh generation
    //*****************************************
    //TODO: currently, this can only be execute by the master rank, change it to general case
    bool CreateMesh();
    void SaveMesh();
    
private:
    bool Create1DMesh();
    bool Create2DMesh();
    bool Create3DMesh();
    //******************************************
    //*** For gmsh related thinigs
    //******************************************
    bool ReadMeshFromGmsh();
    int GetNodesNumViaGmshElmtType(int elmttype) const;
    int GetSurfaceElmtTypeViaGmshBulkElmtType(int elmttype) const;
    int GetElmtDimViaGmshElmtType(int elmttype) const;
    int GetElmtOrderViaGmshElmtType(int elmttype) const;
    string GetElmtNameViaGmshElmtType(int elmttype) const;
    MeshType GetElmtTypeViaGmshElmtType(int elmttype) const;
    MeshType GetBCElmtTypeViaGmshElmtType(int elmttype) const;
    int GetElmtSurfaceNumsViaGmshElmtType(int elmttype) const;
    vector<int> GetIthElmtJthSurfaceConn(const int &elmttype,const int &e,const int &j)const;
    int GetElmtVTKCellTypeViaGmshElmtType(int elmttype) const;
    void ModifyElmtConnViaGmshElmtType(int elmttype,vector<int> &conn) const;

    //******************************************
    //*** For mesh from Abaqus's inp file
    //******************************************
    bool ReadMeshFromAbaqus();
    int GetDimFromInpMeshTypeName(string meshtypename) const;
    int GetVTKCellTypeFormInpMeshTypeName(string meshtypename) const;
    int GetElmtOrderViaInpElmtTypeName(string meshtypename)const;
    int GetAbaqusNodesNumFromInp(string filename) const;
    int GetAbaqusElmtsNumFromInp(string filename) const;
    int GetAbaqusBCElmtsNumFromInp(string filename,int nNodesPerBCElmt) const;
    int GetNodeSetsNumFromInp(string filename) const;
    int GetElmtSetsNumFromInp(string filename) const;
    int GetAbaqusBCElmtNodesNumFromInp(string meshtypename) const;
    int GetElmtNodesNumFromInpElmtName(string meshtypename) const;
    int GetSurfaceElmtNodesNumFromInptElmtName(string meshtypename)const;
    int GetLineElmtNodesNumFromInputElmtName(string meshtypename)const;

    MeshType GetMeshTypeViaAbaqusMeshName(string meshtypename) const;
    MeshType GetBCMeshTypeViaAbaqusMeshName(string meshtypename) const;
    

    string GetElmtTypeNameFromInp(string filename) const;
    vector<string> GetNodeSetsNameFromInp(string filename) const;
    vector<string> GetElmtSetsNameFromInp(string filename) const;
    vector<vector<int>> GetNodeIndexSetsFromInp(string filename) const;
    vector<vector<int>> GetElmtIndexSetsFromInp(string filename) const;

    vector<int> GetNodeIndexVecFromInpNodeSetName(string filename,string nodesetname)const;
    vector<int> GetElmtIndexVecFromInpNodeSetName(string filename,string elmtsetname)const;
    //******************************************
    //*** For print mesh information
    //******************************************
public:
    void PrintMeshInfo() const;
    void PrintMeshDetailInfo() const;
private:
    //************************************************
    //*** For basic mesh information, i.e. node coords and connectivity
    //************************************************
    bool _IsMeshCreated;
    string _MeshFileName;
    string _GmshFileName,_AbaqusFileName;// for gmsh import
    vector<int> _ElmtVTKCellType;PetscInt _BulkElmtVTKCellType;
    vector<double> _NodeCoords;// store all the node coornidates
    vector<vector<int>> _ElmtConn;// store all the elements' connectivity
    vector<vector<int>> _BulkElmtConn;// store the bulk element's connectivity
                                           // if 3D, then surface element conn is not empty
                                           // if 2D, bulk element conn is already the surface element conn
                                           // thus only line element conn is not empty
                                           // if 1D, then only bulk element conn is not empty
    vector<vector<int>> _LineElmtConn,_SurfaceElmtConn;
    MeshType _BulkMeshType,_SurfaceMeshType,_LineMeshType;
    bool _IsBuiltInMesh=true,_IsGmshMesh=false,_IsAbaqusMesh=false;

    //****************************************
    //*** for gmsh related information
    //****************************************
    vector<int> _PhyGroupDimVec;
    vector<pair<int,string>> _PhyGroupArray;// id<--->physical name pair
    vector<int> _MeshUniGeoID,_MeshUniPhyID,_MeshUniPhyDim;
    vector<int> _ElmtDimVec,_ElmtTypeVec,_ElmtPhyIDVec,_ElmtGeoIDVec;
    
    //*************************************
    //*** For each mesh's uel and umat information
    vector<pair<ElmtType,pair<MateType,int>>> _MeshElmtInfoList;
    vector<vector<int>> _BulkElmtDofIndexList;
private:
    //******************************************
    //*** For mesh's geometry information
    //******************************************
    PetscInt _Nx,_Ny,_Nz;
    PetscReal _Xmin,_Xmax,_Ymin,_Ymax,_Zmin,_Zmax;
    PetscInt _nMaxDim,_nMinDim;
    PetscInt _nElmts,_nBulkElmts,_nNodes;
    PetscInt _nNodesPerBulkElmt,_nNodesPerSurfaceElmt,_nNodesPerLineElmt;
    PetscInt _nOrder;
    vector<PetscReal> _ElmtVolume;
    PetscReal _TotalVolume;
private:
    //*********************************************
    //*** For mesh's physical group management
    //*********************************************
    vector<string> _PhysicGroupNameList;
    vector<pair<string,int>> _PhysicNameToIDList;
    vector<pair<string,int>> _PhysicNameToDimList;
    vector<pair<string,int>> _PhysicNameToElmtNodesNumList;
    vector<pair<int,string>> _PhysicIDToNameList;
    vector<pair<string,vector<int>>> _PhysicNameToElmtIndexSet;
    int _nPhysicGroups;
    string _BulkMeshTypeName;

    
};

#endif // ASFEM_MESH_H