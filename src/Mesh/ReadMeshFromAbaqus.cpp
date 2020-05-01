//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "Mesh/Mesh.h"
#include "Utils/StringUtils.h"

bool Mesh::ReadMeshFromAbaqus(){
    ifstream in;
    in.open(_AbaqusFileName.c_str(),ios::in);
    if(!in.is_open()){
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: can\'t read the .inp file(=%20s)      !!!   ***\n",_GmshFileName.c_str());
        PetscPrintf(PETSC_COMM_WORLD,"***        please make sure file name is correct                !!!   ***\n");
        Msg_AsFem_Exit();
    }


    vector<double> numbers;
    string str,substr;
    vector<int> tempbbulkelmtid;
    _nMaxDim=-1;_nMinDim=4;
    _nPhysicGroups=0;
    _PhysicGroupNameList.clear();
    _PhysicIDToNameList.clear();
    _PhysicNameToIDList.clear();
    _PhysicNameToElmtIndexSet.clear();


    _PhyGroupDimVec.clear();
    _PhyGroupArray.clear();

    _nNodes=GetAbaqusNodesNumFromInp(_AbaqusFileName);
    _nBulkElmts=GetAbaqusElmtsNumFromInp(_AbaqusFileName);
    int nBCElmts=GetAbaqusBCElmtsNumFromInp(_AbaqusFileName);
    _nElmts=_nBulkElmts+nBCElmts;

    _nPhysicGroups=GetNodeSetsNumFromInp(_AbaqusFileName)+1;// plus the bulk
    _PhysicGroupNameList=GetElmtSetsNameFromInp(_AbaqusFileName);
    _PhysicGroupNameList.push_back("alldomain");

    

    
    for(int i=1;i<=(int)_PhysicGroupNameList.size();i++){
        _PhysicIDToNameList.push_back(make_pair(i,_PhysicGroupNameList[i-1]));
        _PhysicNameToDimList.push_back(make_pair(_PhysicGroupNameList[i-1],i));
    }

    _BulkMeshTypeName=GetElmtTypeNameFromInp(_AbaqusFileName);
    int nNodesPerSurfaceElmt=0;
    int nNodesPerLineElmt=0;
    vector<vector<int>> Node2ElmtList;
    while(!in.eof()){
        // now we start to read *.msh file
        getline(in,str);

        if(str.find("*Node")!=string::npos){
            // read the nodes' coordinates
            // node-id, x, y, z
            _NodeCoords.resize(_nNodes*4,0.0);
            Node2ElmtList.resize(_nNodes,vector<int>(0));
            int id=1;
            double x,y,z;

            _Xmax=-1.0e16;_Xmin=1.0e16;
            _Ymax=_Xmax;_Ymin=_Xmin;
            _Zmax=_Xmax;_Zmin=_Xmin;
            for(int i=0;i<_nNodes;i++){
                getline(in,str);
                numbers=SplitStrNum(str);
                id=int(numbers[0]);
                x=numbers[1];y=numbers[2];
                z=0.0;
                if(numbers.size()==4) z=numbers[3];
                _NodeCoords[(id-1)*4+0]=1.0;
                _NodeCoords[(id-1)*4+1]=x;
                _NodeCoords[(id-1)*4+2]=y;
                _NodeCoords[(id-1)*4+3]=z;
                if(x>_Xmax) _Xmax=x;
                if(x<_Xmin) _Xmin=x;
                if(y>_Ymax) _Ymax=y;
                if(y<_Ymin) _Ymin=y;
                if(z>_Zmax) _Zmax=z;
                if(z<_Zmin) _Zmin=z;
            }
        }
        else if(str.find("*Element")!=string::npos){
            // here the element constains all the element(no matter it is line elmt or bulk elmt)
            _ElmtConn.resize(_nElmts,vector<int>(0));
            _ElmtVTKCellType.resize(_nElmts,0);
            _ElmtVolume.resize(_nElmts,0.0);
            _ElmtDimVec.resize(_nElmts,0);
            _ElmtTypeVec.resize(_nElmts,0);
            _ElmtPhyIDVec.resize(_nElmts,0);
            _ElmtGeoIDVec.resize(_nElmts,0);

            tempbbulkelmtid.clear();

            int elmtid,phyid,geoid,elmttype,vtktype;
            int nodes,dim,maxdim,elmtorder,iInd;
            vector<int> tempconn;
            MeshType meshtype,bcmeshtype;
            bool IsElmtCoverThisNode=false;

            _nNodesPerBulkElmt=-1;
            _nNodesPerLineElmt=0;
            _nNodesPerSurfaceElmt=0;
            maxdim=-1;
            _nOrder=1;
            nodes=GetElmtNodesNumFromInpElmtName(_BulkMeshTypeName);
            dim=GetDimFromInpMeshTypeName(_BulkMeshTypeName);
            vtktype=GetVTKCellTypeFormInpMeshTypeName(_BulkMeshTypeName);
            meshtype=GetMeshTypeViaAbaqusMeshName(_BulkMeshTypeName);
            bcmeshtype=GetBCMeshTypeViaAbaqusMeshName(_BulkMeshTypeName);
            elmtorder=GetElmtOrderViaInpElmtTypeName(_BulkMeshTypeName);

            
            if(elmtorder>_nOrder) _nOrder=elmtorder;
            elmtid=0;phyid=1;geoid=1;elmttype=4;

            
            for(int e=1;e<=_nBulkElmts;e++){
                getline(in,str);
                numbers=SplitStrNum(str);
                elmtid=int(numbers[0]);
                

                if(dim==1){
                    _nNodesPerLineElmt=nodes;
                    _nNodesPerBulkElmt=nodes;
                    _LineMeshType=meshtype;
                    _BulkMeshType=meshtype;
                }
                if(dim==2){
                    _nNodesPerSurfaceElmt=nodes;
                    _nNodesPerBulkElmt=nodes;
                    _LineMeshType=bcmeshtype;
                    _SurfaceMeshType=meshtype;
                    _BulkMeshType=meshtype;
                }
                if(dim==3){
                    _nNodesPerBulkElmt=nodes;
                    _SurfaceMeshType=bcmeshtype;
                    _BulkMeshType=meshtype;
                }
                
                if(dim>maxdim) maxdim=dim;
                if(dim>_nMaxDim) _nMaxDim=dim;
                if(dim<_nMinDim) _nMinDim=dim;

                

                _ElmtConn[elmtid-1+nBCElmts].resize(nodes+1,0);
                _ElmtConn[elmtid-1+nBCElmts][0]=nodes;
                for(int j=0;j<nodes;j++){
                    iInd=static_cast<int>(numbers[j+1]);
                    _ElmtConn[elmtid-1+nBCElmts][j+1]=iInd;
                    IsElmtCoverThisNode=true;
                    
                    if(Node2ElmtList[iInd-1].size()<1) IsElmtCoverThisNode=false;
                    for(int i=0;i<Node2ElmtList[iInd-1].size();i++){
                        if(Node2ElmtList[iInd-1][i]==e){
                            IsElmtCoverThisNode=true;
                            break;
                        }
                        else{
                            IsElmtCoverThisNode=false;
                        }
                    }
                    if(!IsElmtCoverThisNode){
                        Node2ElmtList[iInd-1].push_back(e);
                    }
                }
                _ElmtDimVec[elmtid-1+nBCElmts]=dim;
                _ElmtTypeVec[elmtid-1+nBCElmts]=elmttype;
                _ElmtPhyIDVec[elmtid-1+nBCElmts]=phyid;
                _ElmtGeoIDVec[elmtid-1+nBCElmts]=geoid;
                _ElmtVTKCellType[elmtid-1+nBCElmts]=vtktype;
                tempbbulkelmtid.push_back(elmtid+nBCElmts);
                
            }
        }
    }
    in.close();
    
    string phyname;
    vector<int> NodeIndexSets,ElmtIndexSets,tempelmtid;
    
    int count=0;
    int bcnodes=GetAbaqusBCElmtNodesNumFromInp(_BulkMeshTypeName);
    int e,i1,i2;
    _nMinDim=_nMaxDim-1;
    
    for(int i=0;i<_nPhysicGroups-1;i++){
        // now we start to read the bc information
        phyname=_PhysicGroupNameList[i];
        NodeIndexSets=GetNodeIndexVecFromInpNodeSetName(_AbaqusFileName,phyname);
        ElmtIndexSets=GetElmtIndexVecFromInpNodeSetName(_AbaqusFileName,phyname);
        tempelmtid.clear();
        for(e=0;e<(int)ElmtIndexSets.size();e++){
            i1=e*(bcnodes-1);
            i2=i1+bcnodes;
            _ElmtConn[count].resize(bcnodes+1);
            _ElmtConn[count][0]=bcnodes;
            for(i=i1;i<i2;i++){
                _ElmtConn[count][i-i1+1]=NodeIndexSets[i];
            }
            count+=1;
            tempelmtid.push_back(count);
        }
        _PhysicNameToElmtIndexSet.push_back(make_pair(phyname,tempelmtid));
    }

    _PhysicNameToElmtIndexSet.push_back(make_pair("alldomain",tempbbulkelmtid));



    if(count!=nBCElmts){
        return false;
    }
    
    
    return true;
}