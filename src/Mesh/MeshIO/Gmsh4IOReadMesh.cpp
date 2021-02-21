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
//+++ Date   : 2020.06.26
//+++ Purpose: implement the mesh import for msh(ver=2) file format
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Mesh/Gmsh4IO.h"

bool Gmsh4IO::ReadMeshFromFile(Mesh &mesh){
    if(!_HasSetMeshFileName){
        MessagePrinter::PrintErrorTxt("can\'t read mesh, the mesh file name has not been set");
        return false;
    }
    string str;
    _in.open(_MeshFileName.c_str(),ios::in);
    if(!_in.is_open()){
        str="can\'t read the .msh file(="+_MeshFileName+"), please make sure your mesh file name is correct";
        MessagePrinter::PrintErrorTxt(str);
        MessagePrinter::AsFem_Exit();
    }

    double version;
    int format,size;
    vector<double> numbers;

    //******************************************************
    //*** initialize all the arraies
    //******************************************************
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
    //*** for nodeset physical group information
    mesh.GetBulkMeshNodeSetPhysicalNameListPtr().clear();
    mesh.GetBulkMeshNodeSetPhysicalIDListPtr().clear();
    mesh.GetBulkMeshNodeSetPhysicalID2NameListPtr().clear();
    mesh.GetBulkMeshNodeSetPhysicalName2IDListPtr().clear();
    mesh.GetBulkMeshNodeSetPhysicalName2NodeIDsListPtr().clear();

    mesh.SetBulkMeshDim(0);
    mesh.SetBulkMeshNodesNum(0);
    mesh.SetBulkMeshElmtsNum(0);
    mesh.SetBulkMeshBulkElmtsNum(0);
    mesh.SetBulkMeshSurfaceElmtsNum(0);
    mesh.SetBulkMeshLineElmtsNum(0);
    mesh.SetBulkMeshPhysicalGroupNums(0);
    mesh.SetBulkMeshNodeSetPhysicalGroupNums(0);

    _nMaxDim=-1;_nMinDim=4;
    _nPhysicGroups=0;
    _nNodeSetPhysicalGroups=0;

    vector<pair<int,int>> UniquePhyDim2IDList;
    int MaxPhyIDofPhyGroup=-1,MaxPhyIDofElmt=-1;

    vector<double> NodeCoords;
    vector<int>    NodeIDFlag;//used to remove the discontinue node id
    vector<vector<int>> ElmtConn;
    vector<int> ElmtDimVec,ElmtPhyIDVec,ElmtIDFlag,ElmtVTKCellType;
    vector<MeshType> ElmtTypeVec;
    int nNodes,nElmts;
    int numEntityBlocks,numNodes,minNodeTag,maxNodeTag;
    int numElements,minElementTag,maxElementTag;

    numEntityBlocks=numNodes=minNodeTag=maxNodeTag=0;
    numElements=minElementTag=maxElementTag=0;

    NodeCoords.clear();ElmtConn.clear();
    _PointsEntityPhyIDs.clear();
    _CurvesEntityPhyIDS.clear();
    _SurfaceEntityPhyIDs.clear();
    _VolumesEntityPhyIDs.clear();

    _Xmax=-1.0e16;_Xmin=1.0e16;
    _Ymax=-1.0e16;_Ymin=1.0e16;
    _Ymax=-1.0e16;_Ymin=1.0e16;

    map<int,vector<int>> NodeSetPhyID2NodeIDsList;
    //*****************************************************
    //*** start to read mesh file
    //*****************************************************
    while(!_in.eof()){
        // now we start to read *.msh file
        getline(_in,str);
        if(str.find("$MeshFormat")!=string::npos) {
            _in >> version >> format >> size;
            if ((version != 4.0) && (version != 4.1) && (version != 4.2)) {
                // currently, only the gmsh2 format is supported!!!
                str = "version=" + to_string(version) + " is not supported yet";
                MessagePrinter::PrintErrorTxt(str);
                return false;
            }
        }
        else if(str.find("$PhysicalNames")!=string::npos) {
            _nPhysicGroups = 0;
            mesh.GetBulkMeshPhysicalGroupNameListPtr().clear();
            mesh.GetBulkMeshPhysicalGroupIDListPtr().clear();
            mesh.GetBulkMeshPhysicalGroupDimListPtr().clear();
            mesh.GetBulkMeshPhysicalGroupID2NameListPtr().clear();
            mesh.GetBulkMeshPhysicalGroupName2IDListPtr().clear();
            mesh.GetBulkMeshPhysicalName2ElmtIDsListPtr().clear();
            mesh.GetBulkMeshNodeSetPhysicalName2NodeIDsListPtr().clear();
            mesh.GetBulkMeshPhysicalGroupName2NodesNumPerElmtListPtr().clear();
            mesh.SetBulkMeshPhysicalGroupNums(0);
            //******************************************************
            //*** node-set physical information
            mesh.GetBulkMeshNodeSetPhysicalNameListPtr().clear();
            mesh.GetBulkMeshNodeSetPhysicalIDListPtr().clear();
            mesh.GetBulkMeshNodeSetPhysicalID2NameListPtr().clear();
            mesh.GetBulkMeshNodeSetPhysicalName2IDListPtr().clear();
            mesh.GetBulkMeshNodeSetPhysicalName2NodeIDsListPtr().clear();
            mesh.SetBulkMeshNodeSetPhysicalGroupNums(0);

            int phydim, phyid;
            string phyname;
            _in >> _nPhysicGroups;
            getline(_in, str);//remove \n in this line
            for (int i = 0; i < _nPhysicGroups; i++) {
                getline(_in, str);
                istringstream s_stream(str);
                s_stream >> phydim >> phyid >> phyname;
                // remove the '"' ,keep only the text
                phyname.erase(std::remove(phyname.begin(), phyname.end(), '"'), phyname.end());
                mesh.GetBulkMeshPhysicalGroupNameListPtr().push_back(phyname);
                mesh.GetBulkMeshPhysicalGroupIDListPtr().push_back(phyid);
                mesh.GetBulkMeshPhysicalGroupDimListPtr().push_back(phydim);

                mesh.GetBulkMeshPhysicalGroupID2NameListPtr().push_back(make_pair(phyid, phyname));
                mesh.GetBulkMeshPhysicalGroupName2IDListPtr().push_back(make_pair(phyname, phyid));


                if (phydim > _nMaxDim) _nMaxDim = phydim;
                if (_nMinDim < phydim) _nMinDim = phydim;
                if (phyid > MaxPhyIDofPhyGroup) MaxPhyIDofPhyGroup = phyid;
                if(phydim==0){
                    // for nodal physical information
                    _nNodeSetPhysicalGroups+=1;
                    mesh.GetBulkMeshNodeSetPhysicalNameListPtr().push_back(phyname);
                    mesh.GetBulkMeshNodeSetPhysicalIDListPtr().push_back(phyid);
                    mesh.GetBulkMeshNodeSetPhysicalID2NameListPtr().push_back(make_pair(phyid,phyname));
                    mesh.GetBulkMeshNodeSetPhysicalName2IDListPtr().push_back(make_pair(phyname,phyid));
                }
            }
        } // end-of-physical group information
        else if(str.find("$Entities")!=string::npos){
            // read the entities block
            getline(_in,str);
            numbers=StringUtils::SplitStrNum(str);
            int numPoints,numCurves,numSurfaces,numVolumes;
            if(numbers.size()!=4){
                MessagePrinter::PrintErrorTxt("invalid entity block information in your msh4 file inside the $Entities");
                MessagePrinter::AsFem_Exit();
            }
            numPoints=static_cast<int>(numbers[1-1]);
            numCurves=static_cast<int>(numbers[2-1]);
            numSurfaces=static_cast<int>(numbers[3-1]);
            numVolumes=static_cast<int>(numbers[4-1]);

            _PointsEntityPhyIDs.resize(numPoints);
            _CurvesEntityPhyIDS.resize(numCurves);
            _SurfaceEntityPhyIDs.resize(numSurfaces);
            _VolumesEntityPhyIDs.resize(numVolumes);

            int i;
            int nodeid,curveid,surfaceid,volumeid;

            //********************************************
            //*** read points entities
            for(i=0;i<numPoints;i++){
                getline(_in,str);
                numbers=StringUtils::SplitStrNum(str);
                nodeid=static_cast<int>(numbers[1-1]);
                if(static_cast<int>(numbers[5-1])>0){
                    _PointsEntityPhyIDs[nodeid-1]=static_cast<int>(numbers[6-1]);
                }
            }
            //*** read curvess entities
            for(i=0;i<numCurves;i++){
                getline(_in,str);
                numbers=StringUtils::SplitStrNum(str);
                curveid=static_cast<int>(numbers[1-1]);
                if(static_cast<int>(numbers[8-1])>0){
                    _CurvesEntityPhyIDS[curveid-1]=static_cast<int>(numbers[9-1]);
                }
            }
            //*** read surfaces entities
            for(i=0;i<numSurfaces;i++){
                getline(_in,str);
                numbers=StringUtils::SplitStrNum(str);
                surfaceid=static_cast<int>(numbers[1-1]);
                if(static_cast<int>(numbers[8-1])>0){
                    _SurfaceEntityPhyIDs[surfaceid-1]=static_cast<int>(numbers[9-1]);
                }
            }
            //*** read volumes entities
            for(i=0;i<numVolumes;i++){
                getline(_in,str);
                numbers=StringUtils::SplitStrNum(str);
                volumeid=static_cast<int>(numbers[1-1]);
                if(static_cast<int>(numbers[8-1])>0){
                    _VolumesEntityPhyIDs[volumeid-1]=static_cast<int>(numbers[9-1]);
                }
            }
        }
        else if(str.find("$Nodes")!=string::npos){
            // start to read node coordinates
            getline(_in,str); // read the entities
            numbers=StringUtils::SplitStrNum(str);
            if(numbers.size()!=4){
                MessagePrinter::PrintErrorTxt("invalid entities information in $Nodes block of your msh4 file");
                MessagePrinter::AsFem_Exit();
            }
            numEntityBlocks=static_cast<int>(numbers[1-1]);
            numNodes=static_cast<int>(numbers[2-1]);mesh.SetBulkMeshNodesNum(numNodes);
            minNodeTag=static_cast<int>(numbers[3-1]);
            maxNodeTag=static_cast<int>(numbers[4-1]);
            NodeCoords.resize(3*maxNodeTag,0.0);// here the node id may not be contineous case !!!
            NodeIDFlag.resize(maxNodeTag,0);
            vector<int> nodeid;
            double x,y,z;
            int i,j,numNodesInBlock;
            nNodes=0;
            for(int nBlock=0;nBlock<numEntityBlocks;nBlock++){
                nodeid.clear();
                // entityDim(int) entityTag(int) parametric(int; 0 or 1) numNodesInBlock(size_t)
                getline(_in,str);
                numbers=StringUtils::SplitStrNum(str);
                if(numbers.size()!=4){
                    cout<<"wokr, str="<<str<<endl;
                    MessagePrinter::PrintErrorTxt("invalid node entities in your msh4 file inside the $Nodes block");
                    MessagePrinter::AsFem_Exit();
                }
                numNodesInBlock=static_cast<int>(numbers[4-1]);
                for(i=0;i<numNodesInBlock;i++){
                    getline(_in,str); // read the node id
                    numbers=StringUtils::SplitStrNum(str);
                    if(static_cast<int>(numbers[0])<minNodeTag||static_cast<int>(numbers[0])>maxNodeTag){
                        MessagePrinter::PrintErrorTxt("invalid node Tag in your msh4 file inside the $Nodes");
                        MessagePrinter::AsFem_Exit();
                    }
                    else{
                        nodeid.push_back(static_cast<int>(numbers[0]));// store the node id
                    }
                }
                for(i=0;i<numNodesInBlock;i++){
                    getline(_in,str);
                    numbers=StringUtils::SplitStrNum(str);
                    if(numbers.size()!=3){
                        MessagePrinter::PrintErrorTxt("invalid node coordinates information in your msh4 file inside the $Nodes block");
                        MessagePrinter::AsFem_Exit();
                    }
                    x=numbers[0];y=numbers[1];z=numbers[2];
                    j=nodeid[i];
                    NodeCoords[(j-1)*3+0]=x;
                    NodeCoords[(j-1)*3+1]=y;
                    NodeCoords[(j-1)*3+2]=z;
                    nNodes+=1;
                    NodeIDFlag[j-1]=1;

                    if(x>_Xmax) _Xmax=x;
                    if(x<_Xmin) _Xmin=x;
                    if(y>_Ymax) _Ymax=y;
                    if(y<_Ymin) _Ymin=y;
                    if(z>_Zmax) _Zmax=z;
                    if(z<_Zmin) _Zmin=z;
                }
            }
            if(nNodes!=numNodes){
                MessagePrinter::PrintErrorTxt("something is wrong in your msh4 file inside the $Nodes block, nodes numer is not match with the first line");
                MessagePrinter::AsFem_Exit();
            }
        }// end-of-nodes-read
        else if(str.find("$Elements")!=string::npos){
            // read the element information
            getline(_in,str);// read the element total entities
            numbers=StringUtils::SplitStrNum(str);
            if(numbers.size()!=4){
                //numEntityBlocks(size_t) numElements(size_t) minElementTag(size_t) maxElementTag(size_t)
                MessagePrinter::PrintErrorTxt("invalid element entities in your msh4 file inside the $Elements");
                MessagePrinter::AsFem_Exit();
            }
            numEntityBlocks=static_cast<int>(numbers[1-1]);
            numElements=static_cast<int>(numbers[2-1]);
            minElementTag=static_cast<int>(numbers[3-1]);
            maxElementTag=static_cast<int>(numbers[4-1]);

            nElmts=0;
            _nLineElmts=0;
            _nSurfaceElmts=0;
            _nBulkElmts=0;
            _nNodesPerLineElmt=0;
            _nNodesPerSurfaceElmt=0;
            _nNodesPerBulkElmt=0;

            // here we use the maxElementTag, since the id of elements could be discontineous !!!
            ElmtConn.resize(maxElementTag,vector<int>(0));
            ElmtDimVec.resize(maxElementTag,0);
            ElmtPhyIDVec.resize(maxElementTag,0);
            ElmtTypeVec.resize(maxElementTag,MeshType::NULLTYPE);
            ElmtIDFlag.resize(maxElementTag,0);
            ElmtVTKCellType.resize(maxElementTag,0);

            mesh.SetBulkMeshElmtsNum(numElements);

            int entityDim,entityTag,elementType,numElementsInBlock;
            int i,j,elmtid,nodesnum,phyid,vtktype,order;
            MeshType meshtype;
            bool IsUnique;
            order=0;
            for(int nBlocks=0;nBlocks<numEntityBlocks;nBlocks++){
                getline(_in,str);//get current element entities
                //entityDim(int) entityTag(int) elementType(int; see below) numElementsInBlock(size_t)
                numbers=StringUtils::SplitStrNum(str);
                if(numbers.size()!=4){
                    MessagePrinter::PrintErrorTxt("invalid element entities for current element in your msh4 file inside the $Elements");
                    MessagePrinter::AsFem_Exit();
                }
                entityDim=static_cast<int>(numbers[1-1]);
                entityTag=static_cast<int>(numbers[2-1]);
                elementType=static_cast<int>(numbers[3-1]);
                numElementsInBlock=static_cast<int>(numbers[4-1]);
                phyid=GetPhysicalIDViaEntityTag(entityDim,entityTag);
                vtktype=GetElmtVTKCellTypeFromGmshElmtType(elementType);
                meshtype=GetElmtMeshTypeFromGmshElmtType(elementType);
                order=GetElmtOrderFromGmshElmtType(elementType);

                if(UniquePhyDim2IDList.size()==0){
                    UniquePhyDim2IDList.push_back(make_pair(entityDim,phyid));
                }
                else{
                    IsUnique=true;
                    for(const auto &it:UniquePhyDim2IDList){
                        if(it.first==entityDim && it.second==phyid){
                            IsUnique=false;
                            break;
                        }
                    }
                    if(IsUnique){
                        UniquePhyDim2IDList.push_back(make_pair(entityDim,phyid));
                    }
                }

                if(entityDim==1){
                    _nLineElmts+=1;
                    _nNodesPerLineElmt=GetElmtNodesNumFromGmshElmtType(elementType);
                }
                if(entityDim==2){
                    _nSurfaceElmts+=1;
                    _nNodesPerSurfaceElmt=GetElmtNodesNumFromGmshElmtType(elementType);
                }
                if(entityDim==3){
                    _nNodesPerBulkElmt=GetElmtNodesNumFromGmshElmtType(elementType);
                }

                if(entityDim>_nMaxDim) _nMaxDim=entityDim;
                if(entityDim<_nMinDim) _nMinDim=entityDim;

                if(entityTag>MaxPhyIDofElmt) MaxPhyIDofElmt=entityTag;

                for(i=0;i<numElementsInBlock;i++){
                    getline(_in,str);
                    numbers=StringUtils::SplitStrNum(str);
                    nodesnum=static_cast<int>(numbers.size()-1);
                    elmtid=static_cast<int>(numbers[1-1]);
                    if(elmtid<minElementTag||elmtid>maxElementTag){
                        MessagePrinter::PrintErrorTxt("invalid element Tag in your msh4 file inside the $Elements");
                        MessagePrinter::AsFem_Exit();
                    }
                    ElmtConn[elmtid-1].clear();
                    ElmtConn[elmtid-1].resize(1+nodesnum);
                    ElmtConn[elmtid-1][0]=nodesnum;

                    ElmtDimVec[elmtid-1]=entityDim;
                    ElmtPhyIDVec[elmtid-1]=phyid;
                    ElmtVTKCellType[elmtid-1]=vtktype;
                    ElmtTypeVec[elmtid-1]=meshtype;
                    ElmtIDFlag[elmtid-1]=1;
                    for(j=1;j<=nodesnum;j++){
                        ElmtConn[elmtid-1][j]=static_cast<int>(numbers[j]);
                    }
                    nElmts+=1;
                    if(entityDim==0){
                        // for node set
                        if(nodesnum!=1){
                            MessagePrinter::PrintErrorTxt("invalid node set(element) in your msh2 file, the connectivity array should only contain 1 element(node id)");
                            MessagePrinter::AsFem_Exit();
                        }
                        NodeSetPhyID2NodeIDsList[phyid].push_back(static_cast<int>(numbers[1]));
                    }
                }// end-of-local-element-connectivity
            }// end-of-element-entity-block
            if(nElmts!=numElements){
                MessagePrinter::PrintErrorTxt("something is wrong in your msh4 file inside the $Elements block, elements numer is not match with the first line");
                MessagePrinter::AsFem_Exit();
            }
            _nElmts=nElmts;
            mesh.SetBulkMeshDim(_nMaxDim);
            mesh.SetBulkMeshMaxDim(_nMaxDim);
            mesh.SetBulkMeshMinDim(_nMinDim);
            mesh.SetBulkMeshMeshOrder(order);
            if(_nMaxDim==1){
                mesh.SetBulkMeshBulkElmtsNum(_nLineElmts);
                mesh.SetBulkMeshNodesNumPerBulkElmt(_nNodesPerLineElmt);
            }
            else if(_nMaxDim==2){
                mesh.SetBulkMeshBulkElmtsNum(_nSurfaceElmts);
                mesh.SetBulkMeshNodesNumPerBulkElmt(_nNodesPerSurfaceElmt);
                mesh.SetBulkMeshNodesNumPerLineElmt(_nNodesPerLineElmt);
            }
            else if(_nMaxDim==3){
                mesh.SetBulkMeshBulkElmtsNum(_nBulkElmts);
                mesh.SetBulkMeshNodesNumPerSurfaceElmt(_nNodesPerSurfaceElmt);
                mesh.SetBulkMeshNodesNumPerLineElmt(_nNodesPerLineElmt);
                mesh.SetBulkMeshNodesNumPerBulkElmt(_nNodesPerBulkElmt);
            }

        }// end-of-element-reading

    } // end-of-ifstream-read

    //**********************************************************************************
    //*** now we re-arrange the node id and element id, to make them to be continue
    //**********************************************************************************
    int count=0;
    mesh.GetBulkMeshNodeCoordsPtr().resize(numNodes*3,0.0);
    count=0;
    for(int i=0;i<maxNodeTag;i++){
        if(NodeIDFlag[i]>0){
            count+=1;
            NodeIDFlag[i]=count;// the active node id
            mesh.GetBulkMeshNodeCoordsPtr()[(count-1)*3+1-1]=NodeCoords[i*3+0];
            mesh.GetBulkMeshNodeCoordsPtr()[(count-1)*3+2-1]=NodeCoords[i*3+1];
            mesh.GetBulkMeshNodeCoordsPtr()[(count-1)*3+3-1]=NodeCoords[i*3+2];
        }
    }
    int jj;
    for(auto &it:NodeSetPhyID2NodeIDsList){
        for(int i=0;i<static_cast<int>(it.second.size());i++){
            jj=it.second[i];
            it.second[i]=NodeIDFlag[jj-1];
        }
    }
    // now we rearange the element connectivity
    count=0;
    mesh.SetBulkMeshElmtsNum(_nElmts);
    mesh.GetBulkMeshElmtVolumePtr().resize(_nElmts,0.0);
    mesh.GetBulkMeshElmtDimListPtr().resize(_nElmts,0);
    mesh.GetBulkMeshElmtConnPtr().resize(_nElmts,vector<int>(0));
    mesh.GetBulkMeshElmtVTKCellTypeListPtr().resize(_nElmts,0);
    mesh.GetBulkMeshElmtPhyIDListPtr().resize(_nElmts,0);
    mesh.GetBulkMeshElmtMeshTypeListPtr().resize(_nElmts,MeshType::NULLTYPE);
    for(int e=0;e<maxElementTag;e++){
        if(ElmtIDFlag[e]>0){
            count+=1;
            ElmtIDFlag[e]=count;// the active element id
            mesh.GetBulkMeshElmtDimListPtr()[count-1]=ElmtDimVec[e];
            mesh.GetBulkMeshElmtVTKCellTypeListPtr()[count-1]=ElmtVTKCellType[e];
            mesh.GetBulkMeshElmtPhyIDListPtr()[count-1]=ElmtPhyIDVec[e];
            mesh.GetBulkMeshElmtMeshTypeListPtr()[count-1]=ElmtTypeVec[e];
            mesh.GetBulkMeshElmtConnPtr()[count-1].resize(ElmtConn[e][0]+1);
            mesh.GetBulkMeshElmtConnPtr()[count-1][0]=ElmtConn[e][0];
            for(int i=1;i<=ElmtConn[e][0];i++){
                mesh.GetBulkMeshElmtConnPtr()[count-1][i]=ElmtConn[e][i];
            }
        }
    }

    //********************************************************************
    //*** now we start to lable out all the physical groups
    //********************************************************************
    vector<int> bulkconn;
    vector<vector<int>> phy2conn(UniquePhyDim2IDList.size(),vector<int>(0));

    if(_nPhysicGroups<1){
        // this means, no any physical group information is given in your msh file
        // so, we have to mark all the bulk elements as the "alldomain" physical group
        int maxid=-1;
        int phyid,dim;
        string phyname;
        for(const auto &it:UniquePhyDim2IDList){
            dim=it.first;
            phyid=it.second;
            phyname=to_string(phyid);
            if(phyid>maxid) maxid=phyid;
            if(dim==_nMaxDim){
                mesh.GetBulkMeshPhysicalGroupIDListPtr().push_back(phyid);
                mesh.GetBulkMeshPhysicalGroupDimListPtr().push_back(_nMaxDim);
                mesh.GetBulkMeshPhysicalGroupNameListPtr().push_back(phyname);

                mesh.GetBulkMeshPhysicalGroupID2NameListPtr().push_back(make_pair(phyid,phyname));
                mesh.GetBulkMeshPhysicalGroupName2IDListPtr().push_back(make_pair(phyname,phyid));
            }
        }
        _nBulkElmts=0;
        _nSurfaceElmts=0;
        _nLineElmts=0;


        bulkconn.clear();
        for(int e=0;e<_nElmts;e++){
            if(mesh.GetBulkMeshIthElmtDim(e+1)==_nMaxDim){
                _nBulkElmts+=1;
                bulkconn.push_back(e+1);
                for(int j=0;j<static_cast<int>(phy2conn.size());j++){
                    if(mesh.GetBulkMeshIthElmtPhyID(e+1)==UniquePhyDim2IDList[j].second){
                        phy2conn[j].push_back(e+1);
                        break;
                    }
                }
            }
            else{
                // for lower dimension elements
                if(mesh.GetBulkMeshIthElmtDim(e+1)==2){
                    _nSurfaceElmts+=1;
                }
                else if(mesh.GetBulkMeshIthElmtPhyID(e+1)==1){
                    _nLineElmts+=1;
                }
            }
        }

        mesh.SetBulkMeshBulkElmtsNum(_nBulkElmts);
        mesh.SetBulkMeshSurfaceElmtsNum(_nSurfaceElmts);
        mesh.SetBulkMeshLineElmtsNum(_nLineElmts);

        for(int j=0;j<static_cast<int>(phy2conn.size());j++){
            phyname=mesh.GetBulkMeshIthPhysicalName(j+1);
            mesh.GetBulkMeshPhysicalName2ElmtIDsListPtr().push_back(make_pair(phyname,phy2conn[j]));
        }

        mesh.SetBulkMeshPhysicalGroupNums(1+phy2conn.size());
        mesh.GetBulkMeshPhysicalGroupNameListPtr().push_back("alldomain");
        mesh.GetBulkMeshPhysicalGroupDimListPtr().push_back(_nMaxDim);
        mesh.GetBulkMeshPhysicalGroupIDListPtr().push_back(maxid+1);
        mesh.GetBulkMeshPhysicalGroupID2NameListPtr().push_back(make_pair(maxid+1,"alldomain"));
        mesh.GetBulkMeshPhysicalGroupName2IDListPtr().push_back(make_pair("alldomain",maxid+1));
        mesh.GetBulkMeshPhysicalName2ElmtIDsListPtr().push_back(make_pair("alldomain",bulkconn));

    }
    else if(_nPhysicGroups==1&&UniquePhyDim2IDList.size()==1){
        // for the case where only 1 physical group is defined
        _nBulkElmts=0;
        _nSurfaceElmts=0;
        _nLineElmts=0;
        bulkconn.clear();
        int maxid=-1;
        for(int e=0;e<_nElmts;e++){
            if(mesh.GetBulkMeshIthElmtDim(e+1)==_nMaxDim){
                _nBulkElmts+=1;
                bulkconn.push_back(e+1);
                if(mesh.GetBulkMeshIthElmtPhyID(e+1)>maxid) maxid=mesh.GetBulkMeshIthElmtPhyID(e+1);
            }
        }

        mesh.SetBulkMeshBulkElmtsNum(_nBulkElmts);
        mesh.SetBulkMeshSurfaceElmtsNum(_nSurfaceElmts);
        mesh.SetBulkMeshLineElmtsNum(_nLineElmts);

        mesh.SetBulkMeshPhysicalGroupNums(2);
        mesh.GetBulkMeshPhysicalGroupNameListPtr().push_back("alldomain");
        mesh.GetBulkMeshPhysicalGroupIDListPtr().push_back(maxid+1);
        mesh.GetBulkMeshPhysicalGroupName2IDListPtr().push_back(make_pair("alldomain",maxid+1));
        mesh.GetBulkMeshPhysicalGroupID2NameListPtr().push_back(make_pair(maxid+1,"alldomain"));

        mesh.GetBulkMeshPhysicalName2ElmtIDsListPtr().push_back(make_pair(mesh.GetBulkMeshIthPhysicalName(1),bulkconn));
        mesh.GetBulkMeshPhysicalName2ElmtIDsListPtr().push_back(make_pair("alldomain",bulkconn));

    }
    else{
        string phyname;
        bulkconn.clear();
        phy2conn.resize(_nPhysicGroups,vector<int>(0));

        _nBulkElmts=0;
        _nSurfaceElmts=0;
        _nLineElmts=0;

        for(int e=0;e<_nElmts;e++){
            for(int j=0;j<static_cast<int>(phy2conn.size());j++){
                if(mesh.GetBulkMeshIthElmtPhyID(e+1)==mesh.GetBulkMeshIthPhysicalID(j+1)){
                    phy2conn[j].push_back(e+1);
                    break;
                }
            }
            if(mesh.GetBulkMeshIthElmtDim(e+1)==_nMaxDim){
                _nBulkElmts+=1;
                bulkconn.push_back(e+1);
            }
            else{
                // for lower dimension elements
                if(mesh.GetBulkMeshIthElmtDim(e+1)==2){
                    _nSurfaceElmts+=1;
                }
                else if(mesh.GetBulkMeshIthElmtPhyID(e+1)==1){
                    _nLineElmts+=1;
                }
            }
        }

        mesh.SetBulkMeshBulkElmtsNum(_nBulkElmts);
        mesh.SetBulkMeshSurfaceElmtsNum(_nSurfaceElmts);
        mesh.SetBulkMeshLineElmtsNum(_nLineElmts);

        for(int j=0;j<static_cast<int>(phy2conn.size());j++){
            phyname=mesh.GetBulkMeshIthPhysicalName(j+1);
            mesh.GetBulkMeshPhysicalName2ElmtIDsListPtr().push_back(make_pair(phyname,phy2conn[j]));
        }

        mesh.SetBulkMeshPhysicalGroupNums(1+phy2conn.size());
        mesh.GetBulkMeshPhysicalGroupNameListPtr().push_back("alldomain");
        mesh.GetBulkMeshPhysicalGroupDimListPtr().push_back(_nMaxDim);
        mesh.GetBulkMeshPhysicalGroupIDListPtr().push_back(MaxPhyIDofElmt+1);
        mesh.GetBulkMeshPhysicalGroupID2NameListPtr().push_back(make_pair(MaxPhyIDofElmt+1,"alldomain"));
        mesh.GetBulkMeshPhysicalGroupName2IDListPtr().push_back(make_pair("alldomain",MaxPhyIDofElmt+1));
        mesh.GetBulkMeshPhysicalName2ElmtIDsListPtr().push_back(make_pair("alldomain",bulkconn));

    }

    //*********************************************************
    //*** for node-set physical group information
    //*********************************************************
    int phyid;
    string phyname;
    bool HasNodePhyID;
    for(int i=0;i<_nNodeSetPhysicalGroups;i++){
        phyid=mesh.GetBulkMeshNodeSetPhysicalIDListPtr()[i];
        phyname=mesh.GetBulkMeshNodeSetPhysicalNameListPtr()[i];
        HasNodePhyID=false;
        for(auto it:NodeSetPhyID2NodeIDsList){
            if(phyid==it.first){
                mesh.GetBulkMeshNodeSetPhysicalName2NodeIDsListPtr().push_back(make_pair(phyname,it.second));
                HasNodePhyID=true;
                break;
            }
        }
        if(!HasNodePhyID){
            MessagePrinter::PrintErrorTxt("you defined "+phyname+" in your $Physical block of the msh4 file, however, we can not find the related nodes in $Elements block");
            MessagePrinter::AsFem_Exit();
        }
    }

    return _nBulkElmts>0;
}