#ifndef ASFEM_MESH_H
#define ASFEM_MESH_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <set>


#include "Nodes.h"
#include "MeshTypeDefine.h"
#include "MaterialSystem/MaterialType.h"
#include "ElmtSystem/ElmtSystem.h"
#include "ElmtSystem/ElmtType.h"


using namespace std;

class Nodes;
class MeshTypeDefine;

class Mesh
{
public:
    Mesh();

    void CreateMesh();
    void SetMeshElmtInfo(ElmtSystem &elmtSystem);
    void SaveMesh();

private:
    void Create1DMesh();
    void Create2DMesh();
    void Create3DMesh();

// settings for mesh class
public:
    void SetDim(int dim) {_nDim=dim;}
    void SetNx(int nx) {_Nx=nx;}
    void SetNy(int ny) {_Ny=ny;}
    void SetNz(int nz) {_Nz=nz;}
    void SetXmin(double xmin) {_Xmin=xmin;}
    void SetXmax(double xmax) {_Xmax=xmax;}
    void SetYmin(double ymin) {_Ymin=ymin;}
    void SetYmax(double ymax) {_Ymax=ymax;}
    void SetZmin(double zmin) {_Zmin=zmin;}
    void SetZmax(double zmax) {_Zmax=zmax;}
    void SetMeshType(string meshtype) {_MeshType=meshtype;}
    void SetMeshMode(bool flag) {_IsBuiltInMesh=flag;}
    void SetMeshFileName(string filename) {_MeshFileName=filename;}
    void SetMeshOrder(int norder) {_nOrders=norder;}

    // Get information for mesh class
    inline MeshType GetBulkElmtType() const {return _BulkElmtType;}
    inline MeshType GetSurfaceElmtType() const {return _SurfaceElmtType;}
    inline MeshType GetLineElmtType() const {return _LineElmtType;}

    inline int GetDim() const {return _nDim;}
    inline int GetDimMin() const {return _nDimMin;}

    inline double GetXmin() const{return _Xmin;}
    inline double GetXmax() const{return _Xmax;}
    inline double GetYmin() const{return _Ymin;}
    inline double GetYmax() const{return _Ymax;}
    inline double GetZmin() const{return _Zmin;}
    inline double GetZmax() const{return _Zmax;}

    inline int GetNodesNum() const {return _nNodes;}
    inline int GetMeshOrder() const {return _nOrders;}

    inline int GetNodesNumPerBulkElmt() const {return _nNodesPerBulkElmt;}
    inline int GetNodesNumPerSurfaceElmt() const {return _nNodesPerSurfaceElmt;}
    inline int GetNodesNumPerLineElmt() const {return _nNodesPerLineElmt;}

    inline int GetElmtsNum() const {return _nElmts;}
    inline int GetBulkElmtsNum() const {return _nBulkElmts;}
    inline int GetBCElmtsNum() const {return _nElmts-_nBulkElmts;}
    inline string GetMeshType() const {return _MeshType;}

    inline ElmtType GetIthBulkElmtElmtID(const int &e) const{
        return _MeshElmtInfo[e-1].first;
    }
    inline MaterialType GetIthBulkElmtMateID(const int &e)const{
        return _MeshElmtInfo[e-1].second.first;
    }
    inline int GetIthBulkElmtMateIDIndex(const int &e){
        return _MeshElmtInfo[e-1].second.second;
    }
    inline int GetElmtsNumViaPhyName(string phyname) const
    {
        for(auto it:_PhyNameToElmtIndexList)
        {
            if(it.first==phyname)
            {
                return int(it.second.size());
            }
        }
        return 0;
    }
    inline int GetElmtIndexViaPhyName(string phyname,const int &e) const{
        for(auto it:_PhyNameToElmtIndexList)
        {
            if(it.first==phyname)
            {
                return it.second[e-1];
            }
        }
        return 0;
    }
    inline int GetPhyGroupsNum() const {return _nPhyGroups;}
    inline vector<string> GetPhyGroupList() const {return _PhyGroupNameList;}
    inline vector<string> GetBulkGroupList() const {return _BulkPhyGroupNameList;}
    inline vector<string> GetBounGroupList() const {return _BounPhyGroupNameList;}

    inline double GetIthNodeJthCoord(const int &i,const int &j) const {return _NodeCoords[(i-1)*4+j];}
    inline int GetIthElmtNodesNum(const int &e) const {return _ElmtConn[e-1][0];}
    inline int GetIthBulkElmtNodesNum(const int &e) const {return _ElmtConn[e-1+_nElmts-_nBulkElmts][0];}
    inline int GetIthBulkElmtJthConn(const int &e,const int &j) const {return _ElmtConn[e-1+_nElmts-_nBulkElmts][j];}
    inline int GetIthElmtJthConn(const int &i,const int &j) const {return _ElmtConn[i-1][j];}

    inline int GetIthElmtVTKType(int e) const {return _CellVTKType[e-1];}
    inline int GetBulkElmtVTKType() const {return _VTKCellType;}
    inline int GetIthBulkElmtVTKType(int e) const {return _CellVTKType[e-1+_nElmts-_nBulkElmts];}

    inline void GetIthBulkElmtConn(const int &e,vector<int> &elConn)const{
        for(int i=1;i<=GetIthBulkElmtNodesNum(e);++i){
            elConn[i-1]=GetIthBulkElmtJthConn(e,i);
        }
    }
    inline void GetIthElmtNodes(const int &e,Nodes &nodes) const
    {
        for(int i=1;i<=GetIthElmtNodesNum(e);++i)
        {
            nodes(i,0)=GetIthNodeJthCoord(GetIthElmtJthConn(e,i),0);
            nodes(i,1)=GetIthNodeJthCoord(GetIthElmtJthConn(e,i),1);
            nodes(i,2)=GetIthNodeJthCoord(GetIthElmtJthConn(e,i),2);
            nodes(i,3)=GetIthNodeJthCoord(GetIthElmtJthConn(e,i),3);
        }
    }
    inline void GetIthBulkElmtNodes(const int &e,Nodes &nodes) const
    {
        for(int i=1;i<=GetIthBulkElmtNodesNum(e);++i)
        {
            nodes(i,0)=GetIthNodeJthCoord(GetIthBulkElmtJthConn(e,i),0);
            nodes(i,1)=GetIthNodeJthCoord(GetIthBulkElmtJthConn(e,i),1);
            nodes(i,2)=GetIthNodeJthCoord(GetIthBulkElmtJthConn(e,i),2);
            nodes(i,3)=GetIthNodeJthCoord(GetIthBulkElmtJthConn(e,i),3);
        }
    }
    inline vector<int> GetIthBulkElmtDofsIndex(const int &e)const{
        return _BulkElmtDofIndexList[e-1];
    }

    // print out mesh information
    void PrintMeshInfo() const;
    void PrintMeshDetailInfo() const;
// settings for gmsh
public:
    void SetMshFileName(string filename) {_GmshFileName=filename;_IsBuiltInMesh=false;_UseGmsh=true;}
private:
    bool ReadMeshFromGmsh();
    void GmshReaderForOldVersion();
    void GmshReaderForNewVersion();
private:
    int GetNodesNumViaGmshElmtType(int elmttype) const;
    int GetSurfaceElmtTypeViaGmshBulkElmtType(int elmttype) const;
    int GetElmtDimViaGmshElmtType(int elmttype) const;
    string GetElmtNameViaGmshElmtType(int elmttype) const;
    MeshType GetElmtTypeViaGmshElmtType(int elmttype) const;
    MeshType GetBCElmtTypeViaGmshElmtType(int elmttype) const;
    int GetElmtSurfaceNumsViaGmshElmtType(int elmttype) const;
    vector<int> GetIthElmtJthSurfaceConn(const int &elmttype,const int &e,const int &j)const;
    int GetElmtVTKCellTypeViaGmshElmtType(int elmttype) const;
    void ModifyElmtConnViaGmshElmtType(int elmttype,vector<int> &conn) const;
    

private:
    int _nOrders=1;
    int _Nx,_Ny,_Nz,_nDim,_nDimMin;
    double _Xmin,_Xmax,_Ymin,_Ymax,_Zmin,_Zmax;
    string _MeshType;
    bool _IsBuiltInMesh=true;
    bool _SaveMesh=false;
    string _MeshFileName;

// for gmsh information
private:
    string _GmshFileName="";
    bool _UseGmsh=false;
    vector<int> _PhyGroupDimVec;
    vector<pair<int,string>> _PhyGroupArray;// id<--->physical name pair
    vector<int> _MeshUniGeoID,_MeshUniPhyID,_MeshUniPhyDim;
    vector<int> _ElmtDimVec,_ElmtTypeVec,_ElmtPhyIDVec,_ElmtGeoIDVec;

    // Information for mesh (for both the built-in mesh and gmsh)
    bool _IsMeshCreated;
    int _nNodes,_nElmts,_nBulkElmts,_nNodesPerElmt;

    int _nNodesPerBulkElmt=0,_nNodesPerSurfaceElmt=0,_nNodesPerLineElmt=0;
    int _nPhyGroups;
    int _VTKCellType;
    vector<double> _NodeCoords;
    vector<vector<int>> _ElmtConn;
    //vector<int> _BulkElmtConn;
    vector<string> _PhyGroupNameList,_BulkPhyGroupNameList,_BounPhyGroupNameList;
    vector<pair<string,int>> _PhyNameToIDList;
    vector<pair<int,string>> _PhyIDToNameList;
    vector<pair<string,vector<int>>> _PhyNameToElmtIndexList;
    vector<int> _CellVTKType;

    vector<pair<ElmtType,pair<MaterialType,int>>> _MeshElmtInfo;// store the [elmt] block info for 
                                                      // first array store iuel information
                                                      // second array store imate information
    vector<vector<int>> _BulkElmtDofIndexList;// store the active dof index(i.e. :e=1-->1,2;e=2-->2,3)

    MeshType _BulkElmtType,_SurfaceElmtType,_LineElmtType;
};

#endif 