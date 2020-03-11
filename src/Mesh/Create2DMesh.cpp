//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************


#include "Mesh/Mesh.h"

bool Mesh::Create2DMesh(){
    _nNodes=0;_nElmts=0;_nBulkElmts=0;
    _nMinDim=1;_nMaxDim=2;
    _nPhysicGroups=0;
    _NodeCoords.clear();
    _ElmtConn.clear();
    _PhysicGroupNameList.clear();
    _PhysicIDToNameList.clear();
    _PhysicNameToIDList.clear();
    _PhysicNameToElmtIndexSet.clear();
    _ElmtVTKCellType.clear();

    _nNodesPerBulkElmt=0;_nNodesPerSurfaceElmt=0;_nNodesPerLineElmt=0;

    vector<int> tempconn;

    _IsMeshCreated=false;
    //****************************************
    //*** Start to generate mesh
    //****************************************
    double dx,dy;
    int e,i,j,k;
    int i1,i2,i3,i4,i5,i6,i7,i8,i9;
    int nNodesPerBCElmt;

    nNodesPerBCElmt=2;

    if(_BulkMeshType==MeshType::QUAD4){
        _nOrder=1;
        dx=(_Xmax-_Xmin)/(_Nx);
        dy=(_Ymax-_Ymin)/(_Ny);
        _nBulkElmts=_Nx*_Ny;
        _nElmts=_nBulkElmts+2*(_Nx+_Ny);
        _nNodes=(_Nx+1)*(_Ny+1);
        _nNodesPerBulkElmt=4;
        _nNodesPerLineElmt=2;
        nNodesPerBCElmt=2;


        _BulkMeshType=MeshType::QUAD4;
        _LineMeshType=MeshType::EDGE2;

        _nNodesPerBulkElmt=4;
        _nNodesPerLineElmt=2;

        _BulkElmtVTKCellType=9;

        _NodeCoords.resize(4*_nNodes,0.0);
        for(j=1;j<=_Ny+1;++j){
            for(i=1;i<=_Nx+1;++i){
                k=(j-1)*(_Nx+1)+i;
                _NodeCoords[(k-1)*4+0]=1.0;
                _NodeCoords[(k-1)*4+1]=_Xmin+(i-1)*dx;
                _NodeCoords[(k-1)*4+2]=_Ymin+(j-1)*dy;
                _NodeCoords[(k-1)*4+3]=0.0;
            }
        }
        // Connectivity information for bulk element
        _ElmtConn.resize(_nElmts,vector<int>(0));
        _ElmtVTKCellType.resize(_nElmts,0);
        _ElmtVolume.resize(_nElmts,0.0);
        tempconn.clear();
        for(j=1;j<=_Ny;++j){
            for(i=1;i<=_Nx;++i){
                e=(j-1)*_Nx+i;
                i1=(j-1)*(_Nx+1)+i;
                i2=i1+1;
                i3=i2+_Nx+1;
                i4=i3-1;

                _ElmtConn[e-1+_nElmts-_nBulkElmts].clear();
                _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(_nNodesPerBulkElmt);
                _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i1);
                _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i2);
                _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i3);
                _ElmtConn[e-1+_nElmts-_nBulkElmts].push_back(i4);
                _ElmtVTKCellType[e-1+_nElmts-_nBulkElmts]=_BulkElmtVTKCellType;

                _ElmtVolume[e-1+_nElmts-_nBulkElmts]=dx*dy;

                tempconn.push_back(e+_nElmts-_nBulkElmts);
            }
        }
    }
    else if(_BulkMeshType==MeshType::QUAD8){
        _nOrder=2;
        dx=(_Xmax-_Xmin)/(2.0*_Nx);
        dy=(_Ymax-_Ymin)/(2.0*_Ny);

        _nBulkElmts=_Nx*_Ny;
        _nElmts=_nBulkElmts+2*(_Nx+_Ny);
        _nNodes=(2*_Nx+1)*(2*_Ny+1)-_nBulkElmts;
        _nNodesPerBulkElmt=8;
        _nNodesPerLineElmt=2;
        nNodesPerBCElmt=3;


        _BulkMeshType=MeshType::QUAD8;
        _LineMeshType=MeshType::EDGE3;

        _nNodesPerBulkElmt=8;
        _nNodesPerLineElmt=3;

        _BulkElmtVTKCellType=23;



        // Create node
        _NodeCoords.resize(_nNodes*4,0.0);
        for(j=1;j<=_Ny;++j){
            // for bottom line of each element
            for(i=1;i<=2*_Nx+1;i++){
                k=(j-1)*(2*_Nx+1+_Nx+1)+i;

                _NodeCoords[(k-1)*4  ]=1.0;
                _NodeCoords[(k-1)*4+1]=_Xmin+(i-1)*dx;
                _NodeCoords[(k-1)*4+2]=_Ymin+(j-1)*2*dy;
                _NodeCoords[(k-1)*4+3]=0.0;
            }
            // for middle line of each element
            for(i=1;i<=_Nx+1;i++){
                k=(j-1)*(2*_Nx+1+_Nx+1)+2*_Nx+1+i;

                _NodeCoords[(k-1)*4  ]=1.0;
                _NodeCoords[(k-1)*4+1]=_Xmin+(i-1)*2*dx;
                _NodeCoords[(k-1)*4+2]=_Ymin+(j-1)*2*dy+dy;
                _NodeCoords[(k-1)*4+3]=0.0;
            }
        }
        // for the last top line
        j=_Ny+1;
        for(i=1;i<=2*_Nx+1;i++){
            k=(j-1)*(2*_Nx+1+_Nx+1)+i;

            _NodeCoords[(k-1)*4  ]=1.0;
            _NodeCoords[(k-1)*4+1]=_Xmin+(i-1)*dx;
            _NodeCoords[(k-1)*4+2]=_Ymin+(j-1)*2*dy;
            _NodeCoords[(k-1)*4+3]=0.0;
        }
        // Create connectivity matrix
        _ElmtConn.resize(_nElmts,vector<int>(0));
        _ElmtVTKCellType.resize(_nElmts,0);
        _ElmtVolume.resize(_nElmts,0.0);
        tempconn.clear();
        for(j=1;j<=_Ny;j++){
            for(i=1;i<=_Nx;i++){
                e=(j-1)*_Nx+i;
                i1=(j-1)*(2*_Nx+1+_Nx+1)+2*i-1;
                i2=i1+2;
                i3=i2+(2*_Nx+1+_Nx+1);
                i4=i3-2;

                i5=i1+1;
                i6=i2+(2*_Nx+1)-i;
                i7=i3-1;
                i8=i1+(2*_Nx+1)-(i-1);


                _ElmtConn[(e-1)+_nElmts-_nBulkElmts].clear();
                _ElmtConn[(e-1)+_nElmts-_nBulkElmts].push_back(_nNodesPerBulkElmt);
                _ElmtConn[(e-1)+_nElmts-_nBulkElmts].push_back(i1);
                _ElmtConn[(e-1)+_nElmts-_nBulkElmts].push_back(i2);
                _ElmtConn[(e-1)+_nElmts-_nBulkElmts].push_back(i3);
                _ElmtConn[(e-1)+_nElmts-_nBulkElmts].push_back(i4);
                _ElmtConn[(e-1)+_nElmts-_nBulkElmts].push_back(i5);
                _ElmtConn[(e-1)+_nElmts-_nBulkElmts].push_back(i6);
                _ElmtConn[(e-1)+_nElmts-_nBulkElmts].push_back(i7);
                _ElmtConn[(e-1)+_nElmts-_nBulkElmts].push_back(i8);

                _ElmtVTKCellType[e-1+_nElmts-_nBulkElmts]=_BulkElmtVTKCellType;

                tempconn.push_back(e+_nElmts-_nBulkElmts);

                _ElmtVolume[e-1+_nElmts-_nBulkElmts]=2*dx*2*dy;
            }
        }
    }
    else if(_BulkMeshType==MeshType::QUAD9){
        _nOrder=2;
        dx=(_Xmax-_Xmin)/(2.0*_Nx);
        dy=(_Ymax-_Ymin)/(2.0*_Ny);

        _nBulkElmts=_Nx*_Ny;
        _nElmts=_nBulkElmts+2*(_Nx+_Ny);
        _nNodes=(2*_Nx+1)*(2*_Ny+1);
        _nNodesPerBulkElmt=9;
        nNodesPerBCElmt=3;

        _BulkMeshType=MeshType::QUAD9;
        _LineMeshType=MeshType::EDGE3;

        _nNodesPerBulkElmt=9;
        _nNodesPerLineElmt=3;

        _BulkElmtVTKCellType=28;

        _NodeCoords.resize(_nNodes*4,0.0);
        for(j=1;j<=2*_Ny+1;j++){
            for(i=1;i<=2*_Nx+1;i++){
                k=(j-1)*(2*_Nx+1)+i;
                _NodeCoords[(k-1)*4  ]=1.0;
                _NodeCoords[(k-1)*4+1]=_Xmin+(i-1)*dx;
                _NodeCoords[(k-1)*4+2]=_Ymin+(j-1)*dy;
                _NodeCoords[(k-1)*4+3]=0.0;
            }
        }
        // Create Connectivity matrix
        _ElmtConn.resize(_nElmts,vector<int>(0));
        _ElmtVTKCellType.resize(_nElmts,0);
        _ElmtVolume.resize(_nElmts,0.0);
        tempconn.clear();
        for(j=1;j<=_Ny;j++){
            for(i=1;i<=_Nx;i++){
                e=(j-1)*_Nx+i;
                i1=(j-1)*2*(2*_Nx+1)+2*i-1;
                i2=i1+2;
                i3=i2+2*(2*_Nx+1);
                i4=i3-2;

                i5=i1+1;
                i6=i2+(2*_Nx+1);
                i7=i3-1;
                i8=i1+(2*_Nx+1);
                i9=i8+1;


                _ElmtConn[(e-1)+_nElmts-_nBulkElmts].clear();
                _ElmtConn[(e-1)+_nElmts-_nBulkElmts].push_back(_nNodesPerBulkElmt);
                _ElmtConn[(e-1)+_nElmts-_nBulkElmts].push_back(i1);
                _ElmtConn[(e-1)+_nElmts-_nBulkElmts].push_back(i2);
                _ElmtConn[(e-1)+_nElmts-_nBulkElmts].push_back(i3);
                _ElmtConn[(e-1)+_nElmts-_nBulkElmts].push_back(i4);
                _ElmtConn[(e-1)+_nElmts-_nBulkElmts].push_back(i5);
                _ElmtConn[(e-1)+_nElmts-_nBulkElmts].push_back(i6);
                _ElmtConn[(e-1)+_nElmts-_nBulkElmts].push_back(i7);
                _ElmtConn[(e-1)+_nElmts-_nBulkElmts].push_back(i8);
                _ElmtConn[(e-1)+_nElmts-_nBulkElmts].push_back(i9);

                _ElmtVTKCellType[e-1+_nElmts-_nBulkElmts]=_BulkElmtVTKCellType;

                _ElmtVolume[e-1+_nElmts-_nBulkElmts]=2*dx*2*dy;

                tempconn.push_back(e+_nElmts-_nBulkElmts);
            }
        }
    }
    //********************************************************+
    // split the boundary element and bulk element, generate the physical group information
    //********************************************************
    _nPhysicGroups=4+1;
    _PhysicGroupNameList.clear();
    _PhysicGroupNameList.push_back("left");
    _PhysicGroupNameList.push_back("right");
    _PhysicGroupNameList.push_back("bottom");
    _PhysicGroupNameList.push_back("top");
    _PhysicGroupNameList.push_back("alldomain");


    //******************************
    _PhysicNameToIDList.clear();
    _PhysicNameToIDList.push_back(make_pair("left",1));
    _PhysicNameToIDList.push_back(make_pair("right",2));
    _PhysicNameToIDList.push_back(make_pair("bottom",3));
    _PhysicNameToIDList.push_back(make_pair("top",4));
    _PhysicNameToIDList.push_back(make_pair("alldomain",5));

    _PhysicIDToNameList.clear();
    _PhysicIDToNameList.push_back(make_pair(1,"left"));
    _PhysicIDToNameList.push_back(make_pair(2,"right"));
    _PhysicIDToNameList.push_back(make_pair(3,"bottom"));
    _PhysicIDToNameList.push_back(make_pair(4,"top"));
    _PhysicIDToNameList.push_back(make_pair(5,"alldomain"));

    // now re-loop all the boundary element to store them into different bc vectors
    vector<int> left,right,bottom,top;

    _PhysicNameToElmtIndexSet.clear();
    

    int nBCElmts=0;
    // for leftside
    left.clear();
    i=1;
    for(j=1;j<=_Ny;j++){
        e=(j-1)*_Nx+i;
        _ElmtConn[nBCElmts].clear();
        if(nNodesPerBCElmt==2){
            // quad4 case
            _ElmtConn[nBCElmts].push_back(2);
            _ElmtConn[nBCElmts].push_back(GetIthBulkElmtJthConn(e,4));
            _ElmtConn[nBCElmts].push_back(GetIthBulkElmtJthConn(e,1));
            _ElmtVTKCellType[nBCElmts]=3;
        }
        else{
            _ElmtConn[nBCElmts].push_back(3);
            _ElmtConn[nBCElmts].push_back(GetIthBulkElmtJthConn(e,4));
            _ElmtConn[nBCElmts].push_back(GetIthBulkElmtJthConn(e,8));
            _ElmtConn[nBCElmts].push_back(GetIthBulkElmtJthConn(e,1));
            _ElmtVTKCellType[nBCElmts]=4;
        }
        nBCElmts+=1;
        left.push_back(nBCElmts);
    }
    
    // For right side
    right.clear();
    i=_Nx;
    for(j=1;j<=_Ny;j++){
        _ElmtConn[nBCElmts].clear();
        e=(j-1)*_Nx+i;
        if(nNodesPerBCElmt==2){
            // quad4 case
            _ElmtConn[nBCElmts].push_back(2);
            _ElmtConn[nBCElmts].push_back(GetIthBulkElmtJthConn(e,2));
            _ElmtConn[nBCElmts].push_back(GetIthBulkElmtJthConn(e,3));
            _ElmtVTKCellType[nBCElmts]=3;
        }
        else{
            _ElmtConn[nBCElmts].push_back(3);
            _ElmtConn[nBCElmts].push_back(GetIthBulkElmtJthConn(e,2));
            _ElmtConn[nBCElmts].push_back(GetIthBulkElmtJthConn(e,6));
            _ElmtConn[nBCElmts].push_back(GetIthBulkElmtJthConn(e,3));
            _ElmtVTKCellType[nBCElmts]=4;
        }
        nBCElmts+=1;
        right.push_back(nBCElmts);
    }
    // For bottom edge
    bottom.clear();
    j=1;
    for(i=1;i<=_Nx;i++){
        e=(j-1)*_Nx+i;
        _ElmtConn[nBCElmts].clear();
        if(nNodesPerBCElmt==2){
            // quad4 case
            _ElmtConn[nBCElmts].push_back(2);
            _ElmtConn[nBCElmts].push_back(GetIthBulkElmtJthConn(e,1));
            _ElmtConn[nBCElmts].push_back(GetIthBulkElmtJthConn(e,2));
            _ElmtVTKCellType[nBCElmts]=3;
        }
        else{
            _ElmtConn[nBCElmts].push_back(3);
            _ElmtConn[nBCElmts].push_back(GetIthBulkElmtJthConn(e,1));
            _ElmtConn[nBCElmts].push_back(GetIthBulkElmtJthConn(e,5));
            _ElmtConn[nBCElmts].push_back(GetIthBulkElmtJthConn(e,2));
            _ElmtVTKCellType[nBCElmts]=4;
        }
        nBCElmts+=1;
        bottom.push_back(nBCElmts);
    }
    // For top edge
    top.clear();
    j=_Ny;
    for(i=1;i<=_Nx;i++){
        e=(j-1)*_Nx+i;
        _ElmtConn[nBCElmts].clear();
        if(nNodesPerBCElmt==2){
            // quad4 case
            _ElmtConn[nBCElmts].push_back(2);
            _ElmtConn[nBCElmts].push_back(GetIthBulkElmtJthConn(e,3));
            _ElmtConn[nBCElmts].push_back(GetIthBulkElmtJthConn(e,4));
            _ElmtVTKCellType[nBCElmts]=3;
        }
        else{
            _ElmtConn[nBCElmts].push_back(3);
            _ElmtConn[nBCElmts].push_back(GetIthBulkElmtJthConn(e,3));
            _ElmtConn[nBCElmts].push_back(GetIthBulkElmtJthConn(e,7));
            _ElmtConn[nBCElmts].push_back(GetIthBulkElmtJthConn(e,4));
            _ElmtVTKCellType[nBCElmts]=4;
        }
        nBCElmts+=1;
        top.push_back(nBCElmts);
    }

    _PhysicNameToElmtIndexSet.push_back(make_pair("left",left));
    _PhysicNameToElmtIndexSet.push_back(make_pair("right",right));
    _PhysicNameToElmtIndexSet.push_back(make_pair("bottom",bottom));
    _PhysicNameToElmtIndexSet.push_back(make_pair("top",top));
    _PhysicNameToElmtIndexSet.push_back(make_pair("alldomain",tempconn));

    //*********************************
    _PhysicNameToDimList.clear();
    _PhysicNameToDimList.push_back(make_pair("left",1));
    _PhysicNameToDimList.push_back(make_pair("right",1));
    _PhysicNameToDimList.push_back(make_pair("bottom",1));
    _PhysicNameToDimList.push_back(make_pair("top",1));
    _PhysicNameToDimList.push_back(make_pair("alldomain",2));

    _PhysicNameToElmtNodesNumList.clear();
    _PhysicNameToElmtNodesNumList.push_back(make_pair("left",_nNodesPerLineElmt));
    _PhysicNameToElmtNodesNumList.push_back(make_pair("right",_nNodesPerLineElmt));
    _PhysicNameToElmtNodesNumList.push_back(make_pair("bottom",_nNodesPerLineElmt));
    _PhysicNameToElmtNodesNumList.push_back(make_pair("top",_nNodesPerLineElmt));
    _PhysicNameToElmtNodesNumList.push_back(make_pair("alldomain",_nNodesPerBulkElmt));

    _IsMeshCreated=true;
    return _IsMeshCreated;
}