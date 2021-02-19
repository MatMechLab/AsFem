//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.02.18
//+++ Purpose: read mesh information from abaqus inp file
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Mesh/AbaqusIO.h"

bool AbaqusIO::ReadMeshFromFile(Mesh &mesh){
    string str,substr;
    vector<double> numbers;
    if(!_HasSetMeshFileName){
        MessagePrinter::PrintErrorTxt("can\'t read mesh, the mesh file name has not been set");
        return false;
    }
    _in.open(_MeshFileName.c_str(),ios::in);
    if(!_in.is_open()){
        MessagePrinter::PrintErrorTxt("can\'t read the .inp file(="+_MeshFileName+"), please make sure the file name is correct or your inp file is in the same folder as your input file");
        MessagePrinter::AsFem_Exit();
    }
    //************************************************************
    //*** init all the mesh array and physical information
    //************************************************************
    mesh.GetBulkMeshNodeCoordsPtr().clear();
    mesh.GetBulkMeshElmtConnPtr().clear();
    mesh.GetBulkMeshElmtVolumePtr().clear();

    mesh.GetBulkMeshElmtVTKCellTypeListPtr().clear();
    mesh.GetBulkMeshElmtPhyIDListPtr().clear();
    mesh.GetBulkMeshElmtDimListPtr().clear();
    mesh.GetBulkMeshElmtMeshTypeListPtr().clear();
    //*** for physical groups
    mesh.GetBulkMeshPhysicalGroupNameListPtr().clear();
    mesh.GetBulkMeshPhysicalGroupIDListPtr().clear();
    mesh.GetBulkMeshPhysicalGroupDimListPtr().clear();
    mesh.GetBulkMeshPhysicalGroupID2NameListPtr().clear();
    mesh.GetBulkMeshPhysicalGroupName2IDListPtr().clear();
    mesh.GetBulkMeshPhysicalGroupName2NodesNumPerElmtListPtr().clear();
    mesh.GetBulkMeshPhysicalName2ElmtIDsListPtr().clear();
    mesh.GetBulkMeshPhysicalName2NodeIDsListPtr().clear();

    //*********************************
    mesh.SetBulkMeshDim(GetElmtDimFromInp());
    mesh.SetBulkMeshMinDim(GetElmtDimFromInp()-1);
    mesh.SetBulkMeshMeshOrder(GetElmtOrderFromInp());
    mesh.SetBulkMeshNodesNum(GetElmtNodesNumFromInp());
    mesh.SetBulkMeshElmtsNum(0);
    mesh.SetBulkMeshBulkElmtsNum(0);
    mesh.SetBulkMeshSurfaceElmtsNum(0);
    mesh.SetBulkMeshLineElmtsNum(0);
    mesh.SetBulkMeshPhysicalGroupNums(0);

    _nMaxDim=-1;_nMinDim=4;
    _nPhysicGroups=0;

    _nMaxDim=GetElmtDimFromInp();
    _nMinDim=_nMaxDim-1;

    int nLowerDimElmts=0;

    map<string,vector<int>> PhysicalName2NodeIDsSet;
    map<string,int> PhysicalName2IDList;

    PhysicalName2NodeIDsSet.clear();
    PhysicalName2IDList.clear();

    _nNodes=0;_nElmts=0;_nBulkElmts=0;
    while(!_in.eof()){
        // now we start to read inp file
        // remember: inp do not have physical group information explicitly!
        //           so we use a default value for all the bulk elements(physical id=0)
        // We! do! not! allow! multiple! parts! in your inp!!!
        getline(_in,str);
        if(str.find("*Node")!=string::npos){
            // start to read the nodes information
            int nodeid;
            double x,y,z;
            _nNodes=GetNodesNumFromInp();
            mesh.GetBulkMeshNodeCoordsPtr().resize(_nNodes*3,0.0);
            mesh.SetBulkMeshNodesNum(_nNodes);
            for(int i=0;i<_nNodes;i++){
                getline(_in,str);
                numbers=StringUtils::SplitStrNum(str);
                if(static_cast<int>(numbers.size())!=_nMaxDim+1){
                    MessagePrinter::PrintErrorTxt("incompatible node coordinates with your mesh dimension");
                    MessagePrinter::AsFem_Exit();
                }
                nodeid=static_cast<int>(numbers[0]);
                x=numbers[1];y=numbers[2];z=0.0;
                if(_nMaxDim==3) z=numbers[3];
                mesh.GetBulkMeshNodeCoordsPtr()[(nodeid-1)*3+1-1]=x;
                mesh.GetBulkMeshNodeCoordsPtr()[(nodeid-1)*3+2-1]=y;
                mesh.GetBulkMeshNodeCoordsPtr()[(nodeid-1)*3+3-1]=z;

                if(x>_Xmax) _Xmax=x;
                if(x<_Xmin) _Xmin=x;
                if(y>_Ymax) _Ymax=y;
                if(y<_Ymin) _Ymin=y;
                if(z>_Zmax) _Zmax=z;
                if(z<_Zmin) _Zmin=z;
            }
        } // end-of-node-block
        else if(str.find("*Element")!=string::npos){
            _nBulkElmts=GetElmtsNumFromInp();// this is only the bulk elements number
            _nNodesPerBulkElmt=GetElmtNodesNumFromInp();
            nLowerDimElmts=GetSurfaceElmtsNumFromInp();
            _nElmts=_nBulkElmts+nLowerDimElmts;
            mesh.SetBulkMeshElmtsNum(_nElmts);
            mesh.SetBulkMeshBulkElmtsNum(_nBulkElmts);
            mesh.SetBulkMeshNodesNumPerBulkElmt(_nNodesPerBulkElmt);
            if(_nMaxDim==2){
                mesh.SetBulkMeshLineElmtsNum(nLowerDimElmts);
                mesh.SetBulkMeshNodesNumPerLineElmt(GetSubElmtNodesNumFromInp());
                mesh.SetBulkMeshLineMeshType(GetSubElmtMeshTypeFromInp());
            }else if(_nMaxDim==3){
                mesh.SetBulkMeshSurfaceElmtsNum(nLowerDimElmts);
                mesh.SetBulkMeshNodesNumPerSurfaceElmt(GetSubElmtNodesNumFromInp());
                mesh.SetBulkMeshNodesNumPerLineElmt(GetSubSubElmtNodesNumFromInp());

                mesh.SetBulkMeshSurfaceMeshType(GetSubElmtMeshTypeFromInp());
                mesh.SetBulkMeshLineMeshType(GetSubSubElmtMeshTypeFromInp());
            }
            mesh.SetBulkMeshMeshType(GetElmtMeshTypeFromInp());
            mesh.SetBulkMeshMeshTypeName(GetElmtMeshTypeNameFromInp());
            // allocate memory for connectivity and other array
            mesh.GetBulkMeshElmtConnPtr().resize(_nElmts,vector<int>(0));
            mesh.GetBulkMeshElmtVTKCellTypeListPtr().resize(_nElmts,0);
            mesh.GetBulkMeshElmtPhyIDListPtr().resize(_nElmts,0);
            mesh.GetBulkMeshElmtMeshTypeListPtr().resize(_nElmts,MeshType::NULLTYPE);
            mesh.GetBulkMeshElmtDimListPtr().resize(_nElmts,0);
            mesh.GetBulkMeshElmtVolumePtr().resize(_nElmts,0.0);

            int phyid,vtktype,elmtid;
            MeshType meshtype;
            phyid=0;// we use 0 for all the bulk mesh !!!
            vtktype=GetElmtVTKCellTypeFromInp();
            meshtype=GetElmtMeshTypeFromInp();
            for(int e=0;e<_nBulkElmts;e++){
                getline(_in,str);
                numbers=StringUtils::SplitStrNum(str);
                if(static_cast<int>(numbers.size())!=1+_nNodesPerBulkElmt){
                    MessagePrinter::PrintErrorTxt("error detected, your nodes number per bulk element is not match with your mesh type, please check your inp file");
                    MessagePrinter::AsFem_Exit();
                }
                elmtid=static_cast<int>(numbers[0]);

                mesh.GetBulkMeshElmtConnPtr()[elmtid-1].resize(_nNodesPerBulkElmt+1,0);
                mesh.GetBulkMeshElmtConnPtr()[elmtid-1][0]=_nNodesPerBulkElmt;
                for(int j=1;j<=_nNodesPerBulkElmt;j++){
                    mesh.GetBulkMeshElmtConnPtr()[elmtid-1][j]=static_cast<int>(numbers[j]);
                }

                // set bulk elements' other properties
                mesh.GetBulkMeshElmtMeshTypeListPtr()[elmtid-1]=meshtype;
                mesh.GetBulkMeshElmtPhyIDListPtr()[elmtid-1]=phyid;
                mesh.GetBulkMeshElmtVTKCellTypeListPtr()[elmtid-1]=vtktype;
                mesh.GetBulkMeshElmtDimListPtr()[elmtid-1]=_nMaxDim;
                mesh.GetBulkMeshElmtVolumePtr()[elmtid-1]=0.0;
            }
        }// end-of-element-block
        else if(str.find("*Surface,")!=string::npos){
            // *Surface, type=ELEMENT, name=Surf-Left
            string SurfaceSetName;
            nLowerDimElmts=GetSurfaceElmtsNumFromInp();
            int surfaceid,i,edge;
            vector<int> ElmtIDList;
            i=str.find_last_of("=");
            SurfaceSetName=str.substr(i+1,string::npos);
            ElmtIDList=GetSurfaceElmtIDViaSurfaceNameFromInp(SurfaceSetName);
            edge=GetSurfaceEdgeIDViaSurfaceNameFromInp(SurfaceSetName);
            if(edge||ElmtIDList.size()||surfaceid){}
            // TODO: here we need some strategy to split the boundary element from the surface set
        }
        else if(str.find("*Nset,")!=string::npos){
            int i,j;
            string nodesetname;
            i=str.find("=");j=str.find_last_of(",");
            nodesetname=str.substr(i+1,j-i);
            PhysicalName2NodeIDsSet[nodesetname].clear();
            if(str.find("generate")!=string::npos){
                getline(_in,str);
                numbers=StringUtils::SplitStrNum(str);
                if(static_cast<int>(numbers.size())!=3){
                    MessagePrinter::PrintErrorTxt("invalid node id index information, you need min-id,max-id,increament for generate node set. Please check your input file");
                    MessagePrinter::AsFem_Exit();
                }
                for(i=static_cast<int>(numbers[0]);i<=static_cast<int>(numbers[1]);i+=static_cast<int>(numbers[2])){
                    PhysicalName2NodeIDsSet[nodesetname].push_back(i);
                }
            }
            else{
                for(i=0;i<numeric_limits<int>::max()-1;i++){
                    getline(_in,str);
                    if(str.find("*")!=string::npos) break;
                    numbers=StringUtils::SplitStrNum(str);
                    for(const auto &val:numbers){
                        PhysicalName2NodeIDsSet[nodesetname].push_back(static_cast<int>(val));
                    }
                }
            }
        }
    }// end-of-file-reading

    _nPhysicGroups=0;
    PhysicalName2IDList.clear();
    for(auto it:PhysicalName2NodeIDsSet){
        _nPhysicGroups+=1;
        PhysicalName2IDList[it.first]=_nPhysicGroups;
    }


    return false;
}