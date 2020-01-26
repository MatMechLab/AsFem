#include "Mesh/Mesh.h"

Mesh::Mesh()
{
    _nOrders=1;
    _Nx=2;_Ny=2;_Nz=2;_nDim=2;
    _Xmin=0.0;_Xmax=1.0;
    _Ymin=0.0;_Ymax=1.0;
    _Zmin=0.0;_Zmax=1.0;
    _MeshType="quad4";
    _IsBuiltInMesh=true;


    _nNodesPerBulkElmt=0;_nNodesPerSurfaceElmt=0;_nNodesPerLineElmt=0;

    // for mesh output
    _SaveMesh=false;
    _MeshFileName="";
    // gmsh information
    _GmshFileName="";
    _UseGmsh=false;

    // mesh information
    _nNodes=0;_nElmts=0;
    _nPhyGroups=0;
    _NodeCoords.clear();
    _ElmtConn.clear();
    _PhyGroupNameList.clear();
    _PhyNameToIDList.clear();
    _PhyIDToNameList.clear();
    _PhyNameToElmtIndexList.clear();

    _MeshElmtInfo.clear();

}