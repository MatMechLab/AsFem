//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2019
//* https://github.com/walkandthinker/AsFem
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
    void SetDim(const PetscInt &ndim) {_nMaxDim=ndim;_nMinDim=ndim-1;}
    void SetMaxDim(const PetscInt &ndim) {_nMaxDim=ndim;}
    void SetMinDim(const PetscInt &ndim) {_nMinDim=ndim;}
    //*** for element numbers along different axis
    void SetNx(const PetscInt &nx){_Nx=nx;}
    void SetNy(const PetscInt &ny){_Ny=ny;}
    void SetNz(const PetscInt &nz){_Nz=nz;}
    //*** for element geometry settings
    void SetXmin(const PetscReal &xmin){_Xmin=xmin;}
    void SetXmax(const PetscReal &xmax){_Xmax=xmax;}
    void SetYmin(const PetscReal &ymin){_Ymin=ymin;}
    void SetYmax(const PetscReal &ymax){_Ymax=ymax;}
    void SetZmin(const PetscReal &zmin){_Zmin=zmin;}
    void SetZmax(const PetscReal &zmax){_Zmax=zmax;}
    //*** for mesh type setting
    void SetMeshType(MeshType type){_BulkMeshType=type;}
    //*** for save mesh
    void SetMeshFileName(string name){_MeshFileName=name;}
    //*** set mesh type
    void SetMeshType(string meshtype);
    //***  for mesh mode(true for built-in, false for external one)
    void SetMeshMode(bool flag){_IsBuiltInMesh=flag;if(flag)_IsGmshMesh=false;}
    void SetMeshGmshMode(bool flag){_IsGmshMesh=flag;if(flag) _IsBuiltInMesh=false;}
    //*** set msh file name for gmsh import
    void SetMshFileName(string filename){_GmshFileName=filename;}
    //*** set ith element's volume
    void SetIthElmtVolume(const PetscInt &i,const PetscReal &volume){
        _ElmtVolume[i-1]=volume;
    }
    //*** set i-th bulk elmt volume
    void SetIthBulkElmtVolume(const PetscInt &i,const PetscReal &volume){
        _ElmtVolume[i+_nElmts-_nBulkElmts-1]=volume;
    }
    //*** if one set the volume for each element, the you should update the total volume
    void SetTotalVolume(const PetscReal &volume){_TotalVolume=volume;}

    //*****************************************
    //*** some common getting functions
    //*****************************************
    inline PetscInt GetDim()const{return _nMaxDim;}
    inline PetscInt GetMaxDim()const{return _nMaxDim;}
    inline PetscInt GetMinDim()const{return _nMinDim;}
    //*** for element numbers
    inline PetscInt GetNx()const{return _Nx;}
    inline PetscInt GetNy()const{return _Ny;}
    inline PetscInt GetNz()const{return _Nz;}
    //*** for geometry
    inline PetscReal GetXmin()const{return _Xmin;}
    inline PetscReal GetXmax()const{return _Xmax;}
    inline PetscReal GetYmin()const{return _Ymin;}
    inline PetscReal GetYmax()const{return _Ymax;}
    inline PetscReal GetZmin()const{return _Zmin;}
    inline PetscReal GetZmax()const{return _Zmax;}
    //*** for mesh type
    inline MeshType GetBulkMeshType()const{return _BulkMeshType;}
    inline MeshType GetSurfaceMeshType()const{return _SurfaceMeshType;}
    inline MeshType GetLineMeshType()const{return _LineMeshType;}
    inline PetscInt GetMeshOrder()const{return _nOrder;}
    //*** for nodes number and element numbers
    inline PetscInt GetNodesNum()const{return _nNodes;}
    inline PetscInt GetNodesNumPerBulkElmt()const{return _nNodesPerBulkElmt;}
    inline PetscInt GetNodesNumPerLineElmt()const{return _nNodesPerLineElmt;}
    inline PetscInt GetNodesNumPerSurfaceElmt()const{return _nNodesPerSurfaceElmt;}
    inline PetscInt GetElmtsNum()const{return _nElmts;}
    inline PetscInt GetBulkElmtsNum()const{return _nBulkElmts;}
    inline PetscInt GetBCElmtsNum()const{return _nElmts-_nBulkElmts;}

    //*** for bc elements information
    inline PetscInt GetDimViaPhyName(string bcname)const{
        for(auto it:_PhysicNameToDimList){
            if(it.first==bcname){
                return it.second;
            }
        }
        return -1;
    }
    inline PetscInt GetBCElmtNodesNumViaPhyName(string bcname)const{
        for(auto it:_PhysicNameToElmtNodesNumList){
            if(it.first==bcname){
                return it.second;
            }
        }
        return -1;
    }

    //**** for physical groups
    inline PetscInt GetPhysicalGroupNum()const{return _nPhysicGroups;}
    //**** for i-th element information
    inline PetscInt GetIthBulkElmtVTKType(const PetscInt &i)const{
        return _ElmtVTKCellType[i+GetBCElmtsNum()-1];
    }
    inline ElmtType GetIthBulkElmtElmtType(const PetscInt &i)const{
        return _MeshElmtInfoList[i-1].first;
    }
    inline MateType GetIthBulkElmtMateType(const PetscInt &i)const{
        return _MeshElmtInfoList[i-1].second.first;
    }
    inline int GetIthBulkElmtMateIndex(const PetscInt &i)const{
        return _MeshElmtInfoList[i-1].second.second;
    }
    inline PetscInt GetIthElmtNodesNum(const PetscInt &i)const{
        return _ElmtConn[i-1][0];
    }
    inline PetscInt GetIthBulkElmtNodesNum(const PetscInt &i)const{
        return _ElmtConn[i-1+_nElmts-_nBulkElmts][0];
    }
    inline PetscInt GetElmtsNumViaPhyName(string phyname)const{
        for(auto it:_PhysicNameToElmtIndexSet){
            if(it.first==phyname){
                return PetscInt(it.second.size());
            }
        }
        return 0;
    }
    inline PetscInt GetIthElmtIndexViaPhyName(string phyname,const PetscInt &i)const{
        for(auto it:_PhysicNameToElmtIndexSet){
            if(it.first==phyname){
                return it.second[i-1];
            }
        }
        return 0;
    }
    inline PetscInt GetIthElmtNodesNumViaPhyName(string phyname,const PetscInt &i)const{
        for(auto it:_PhysicNameToElmtIndexSet){
            if(it.first==phyname){
                return GetIthElmtNodesNum(it.second[i-1]);
            }
        }
        return 0;
    }
    //*** for i-th element connectivity
    inline void GetIthBulkElmtConn(const PetscInt &i,vector<PetscInt> &elconn)const{
        for(PetscInt j=1;j<=_ElmtConn[i-1+GetBCElmtsNum()][0];++j){
            elconn[j-1]=_ElmtConn[i-1+GetBCElmtsNum()][j];
        }
    }
    vector<PetscInt> GetIthBulkElmtDofsIndex(const PetscInt &i)const{
        return _BulkElmtDofIndexList[i-1];
    }
    inline PetscInt GetIthElmtJthConn(const PetscInt &i,const PetscInt &j)const{
        return _ElmtConn[i-1][j];
    }
    void GetIthElmtNodes(const PetscInt &e,Nodes &nodes)const{
        for(int i=1;i<=GetIthElmtNodesNum(e);++i){
            nodes(i,0)=GetIthNodeJthCoord(GetIthElmtJthConn(e,i),0);
            nodes(i,1)=GetIthNodeJthCoord(GetIthElmtJthConn(e,i),1);
            nodes(i,2)=GetIthNodeJthCoord(GetIthElmtJthConn(e,i),2);
            nodes(i,3)=GetIthNodeJthCoord(GetIthElmtJthConn(e,i),3);
        }
    }
    void GetIthBulkElmtNodes(const PetscInt &e,Nodes &nodes)const{
        for(int i=1;i<=GetIthBulkElmtNodesNum(e);++i){
            nodes(i,0)=GetIthNodeJthCoord(GetIthBulkElmtJthConn(e,i),0);
            nodes(i,1)=GetIthNodeJthCoord(GetIthBulkElmtJthConn(e,i),1);
            nodes(i,2)=GetIthNodeJthCoord(GetIthBulkElmtJthConn(e,i),2);
            nodes(i,3)=GetIthNodeJthCoord(GetIthBulkElmtJthConn(e,i),3);
        }
    }
    inline PetscInt GetIthBulkElmtJthConn(const PetscInt &i,const PetscInt &j)const{
        return _ElmtConn[i-1+_nElmts-_nBulkElmts][j];
    }
    inline PetscInt GetIthElmtVTKCellType(const PetscInt &i)const{
        return _ElmtVTKCellType[i-1];
    }
    inline PetscInt GetIthBulkElmtVTKCellType(const PetscInt &i)const{
        return _ElmtVTKCellType[i-1+_nElmts-_nBulkElmts];
    }

    //*****************************************************
    //*** get information for nodes and connectivity
    //*****************************************************
    inline PetscReal GetIthNodeJthCoord(const PetscInt &i,const PetscInt &j)const{
        return _NodeCoords[(i-1)*4+j];
    }
    //*** for elemental 'volume'
    inline PetscReal GetTotalVolume()const{return _TotalVolume;}
    inline PetscReal GetIthElmtVolume(const PetscInt &i)const{
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
    string _GmshFileName;// for gmsh import
    vector<PetscInt> _ElmtVTKCellType;PetscInt _BulkElmtVTKCellType;
    vector<PetscReal> _NodeCoords;// store all the node coornidates
    vector<vector<PetscInt>> _ElmtConn;// store all the elements' connectivity
    vector<vector<PetscInt>> _BulkElmtConn;// store the bulk element's connectivity
                                           // if 3D, then surface element conn is not empty
                                           // if 2D, bulk element conn is already the surface element conn
                                           // thus only line element conn is not empty
                                           // if 1D, then only bulk element conn is not empty
    vector<vector<PetscInt>> _LineElmtConn,_SurfaceElmtConn;
    MeshType _BulkMeshType,_SurfaceMeshType,_LineMeshType;
    bool _IsBuiltInMesh,_IsGmshMesh;

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
    vector<vector<PetscInt>> _BulkElmtDofIndexList;
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
    vector<pair<string,PetscInt>> _PhysicNameToIDList;
    vector<pair<string,PetscInt>> _PhysicNameToDimList;
    vector<pair<string,PetscInt>> _PhysicNameToElmtNodesNumList;
    vector<pair<PetscInt,string>> _PhysicIDToNameList;
    vector<pair<string,vector<PetscInt>>> _PhysicNameToElmtIndexSet;
    PetscInt _nPhysicGroups;
    string _BulkMeshTypeName;

    
};

#endif // ASFEM_MESH_H