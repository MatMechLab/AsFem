#include "Mesh/Mesh.h"

void Mesh::Create1DMesh()
{
    _nNodes=0;_nElmts=0;_nBulkElmts=0;
    _nPhyGroups=0;
    _NodeCoords.clear();
    _ElmtConn.clear();
    _PhyGroupNameList.clear();
    _PhyNameToIDList.clear();
    _PhyIDToNameList.clear();
    _PhyNameToElmtIndexList.clear();
    _CellVTKType.clear();

    _IsMeshCreated=false;

    _nNodesPerSurfaceElmt=0;_nNodesPerLineElmt=0;

    //**************************************
    //*** Start to create mesh
    //**************************************
    int i,j,e,P;
    _VTKCellType=4;
    if(_MeshType=="edge2")
    {
        P=1;
        _VTKCellType=3;
        _BulkElmtType=MeshType::EDGE2;
        _nNodesPerBulkElmt=2;
        _nOrders=1;
    }
    if(_MeshType=="edge3") 
    {
        P=2;
        _BulkElmtType=MeshType::EDGE3;
        _nNodesPerBulkElmt=3;
        _nOrders=2;
    }
    if(_MeshType=="edge4")
    {
        P=3;
        _BulkElmtType=MeshType::EDGE4;
        _nNodesPerBulkElmt=4;
        _nOrders=3;
    }

    _nBulkElmts=_Nx;
    _nElmts=_nBulkElmts+2;// line element with 2 point element
    _nNodes=_nBulkElmts*P+1;
    _nNodesPerElmt=P+1;
    _nDimMin=0;
    

    double dx=(_Xmax-_Xmin)/(_nNodes-1);


    
    _ElmtConn.resize(_nElmts,vector<int>(0));
    _CellVTKType.resize(_nElmts,0);
    _NodeCoords.resize(4*_nNodes,0.0);
    for(i=0;i<_nNodes;i++)
    {
        _NodeCoords[i*4+0]=1.0;
        _NodeCoords[i*4+1]=_Xmin+i*dx;
        _NodeCoords[i*4+2]=0.0;
        _NodeCoords[i*4+3]=0.0;
    }
    

    // generate the element connectivity information
    _ElmtConn[0].push_back(1);// first component is length
    _ElmtConn[0].push_back(1);// second one is the node index(start from 1, not zero!!!)
    _CellVTKType[1-1]=1;

    _ElmtConn[1].push_back(1);
    _ElmtConn[1].push_back(_nNodes);
    _CellVTKType[2-1]=1;

    

    vector<int> left,right;
    left.clear();left.push_back(1);// node index
    right.clear();right.push_back(2);
    _PhyNameToElmtIndexList.push_back(make_pair("left",left));
    _PhyNameToElmtIndexList.push_back(make_pair("right",right));

    vector<int> tempconn;
    tempconn.clear();
    for(e=0;e<_nBulkElmts;++e)
    {
        _ElmtConn[2+e].clear();
        _ElmtConn[2+e].push_back(_nNodesPerElmt);
        tempconn.push_back(2+e+1);
        for(j=1;j<=_nNodesPerElmt;++j)
        {
            _ElmtConn[2+e].push_back(e*P+j);
        }
        if(_nNodesPerElmt==2)
        {
            _CellVTKType[2+e]=3;
        }
        else
        {
            _CellVTKType[2+e]=4;
        }
    }
    _PhyNameToElmtIndexList.push_back(make_pair("alldomain",tempconn));

    _BulkPhyGroupNameList.clear();
    _BulkPhyGroupNameList.push_back("alldomain");

    _BounPhyGroupNameList.clear();
    _BounPhyGroupNameList.push_back("left");
    _BounPhyGroupNameList.push_back("right");

    
    _nPhyGroups=3;
    _PhyGroupNameList.push_back("left");
    _PhyGroupNameList.push_back("right");
    _PhyGroupNameList.push_back("alldomain");

    _PhyNameToIDList.push_back(make_pair("left",1));
    _PhyNameToIDList.push_back(make_pair("right",2));
    _PhyNameToIDList.push_back(make_pair("alldomain",3));

    _PhyIDToNameList.push_back(make_pair(1,"left"));
    _PhyIDToNameList.push_back(make_pair(2,"right"));
    _PhyIDToNameList.push_back(make_pair(3,"alldomain"));

    _IsMeshCreated=true;


}