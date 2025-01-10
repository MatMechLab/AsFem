//****************************************************************
//* This file is part of the AsFem framework
//* Advanced Simulation kit based on Finite Element Method (AsFem)
//* All rights reserved, Yang Bai/MM-Lab@CopyRight 2020-present
//* https://github.com/MatMechLab/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2024.07.11
//+++ Purpose: Implement the msh file (version-4) import function.
//+++          This mesh file must be the *.msh in version-4.
//+++          For version-2, please use Msh2FileImporter.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FECell/Msh4File2FECellImporter.h"
#include "Utils/StringUtils.h"
#include "MPIUtils/MPIDataBus.h"

int Msh4File2FECellImporter::getMaxMeshDim(const string &filename)const{
    ifstream in;
    string str;
    int maxdim;
    in.open(filename.c_str(),ios::in);
    if(!in.is_open()){
        MessagePrinter::printErrorTxt("can\'t read the .msh file("+filename+"),please make sure file name is correct"
                                      " or you have the access permission");
        return false;
    }
    maxdim=-1;
    while(!in.eof()){
        getline(in,str);
        if(str.find("$Elements")!=string::npos){
            int numEntityBlocks;
            int entityDim,numElementsInBlock;
            vector<double> numbers;

            getline(in,str);
            numbers=StringUtils::splitStrNum(str);

            if(numbers.size()!=4){
                //numEntityBlocks(size_t) numElements(size_t) minElementTag(size_t) maxElementTag(size_t)
                MessagePrinter::printErrorTxt("Invalid element entities in your msh4 file inside the $Elements");
                MessagePrinter::exitAsFem();
            }
            numEntityBlocks=static_cast<int>(numbers[1-1]);

            for(int nblock=0;nblock<numEntityBlocks;nblock++){
                getline(in,str);
                numbers=StringUtils::splitStrNum(str);
                if(numbers.size()!=4){
                    MessagePrinter::printErrorTxt("Invalid element entities for current element in your msh4 file inside the $Elements");
                    MessagePrinter::exitAsFem();
                }
                entityDim=static_cast<int>(numbers[1-1]);
                if(entityDim>maxdim) maxdim=entityDim;
                numElementsInBlock=static_cast<int>(numbers[4-1]);
                for(int i=0;i<numElementsInBlock;i++){
                    getline(in,str);
                }
            }
            break;
        }
    }
    in.close();
    return maxdim;
}

int Msh4File2FECellImporter::getPhysicalIDViaEntityTag(const int &entityDim,const int &entityTag)const{
    if(entityDim==0){
        return m_PointsEntityPhyIDs[entityTag-1];
    }
    else if(entityDim==1){
        return m_CurvesEntityPhyIDS[entityTag-1];
    }
    else if(entityDim==2){
        return m_SurfaceEntityPhyIDs[entityTag-1];
    }
    else if(entityDim==3){
        return m_VolumesEntityPhyIDs[entityTag-1];
    }
    else{
        MessagePrinter::printErrorTxt("entityDim="+to_string(entityDim)+" is invalid for msh4 file, please check your mesh");
        MessagePrinter::exitAsFem();
    }
    return -1;
}

bool Msh4File2FECellImporter::importMeshFile(const string &filename,FECellData &t_celldata){
    ifstream in;
    in.open(filename.c_str(),ios::in);
    if(!in.is_open()){
        MessagePrinter::printErrorTxt("can\'t read the .msh file("+filename+"),please make sure file name is correct"
                                      " or you have the access permission");
        return false;
    }
    

    int rank,size;

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    
    t_celldata.MaxDim=getMaxMeshDim(filename);
    t_celldata.MinDim=0;
    t_celldata.ActiveDofsNum=0;
    t_celldata.TotalDofsNum=0;
    t_celldata.MaxDofsPerNode=0;

    if(rank!=0) in.close();// close the ifstream on other ranks

    if(rank==0){
        int mshMaxDim;

        int mshPhyGroupNum;
        int maxPhyGroupID;
        int maxPhyGroupDim,minPhyGroupDim;
        int mshBulkPhyGroupNum;

        vector<int> mshPhyGroupIDVec,mshPhyGroupDimVec;
        vector<string> mshPhyGroupNameVec;

        vector<int> mshBulkElmtUniquePhyIDVec;
        vector<int> mshBulkElmtPhyIDs;

        int mshNodalPhyGroupNum;
        vector<int> mshNodalPhyGroupIDVec;
        vector<string> mshNodalPhyGroupNameVec;

        map<int,vector<int>> NodeSetPhyID2IndexMap;
        map<int,string> mshPhyID2NameMap;

        vector<int> ElmtPhyIDVec;
        vector<int> ElmtDimVec;
        vector<vector<int>> ElmtConn;
        vector<int> ElmtIDFlag;
        vector<SingleMeshCell> TempCellVec,BulkCellVec;
        int PointElmtsNum;
        string str;
        double version;
        int format;
        vector<double> numbers;
        
        int numNodes=0;
        int minNodeTag=0;
        int maxNodeTag=0;
        vector<double> NodeCoords;// here the node id may not be contineous case !!!
        vector<int> NodeIDFlag;
        
        int numEntityBlocks=0;
        int numElements=0;
        int minElementTag=0;
        int maxElementTag=0;
        int numElementsInBlock=0;
        
        mshMaxDim=getMaxMeshDim(filename);
        
        t_celldata.MaxDim=mshMaxDim;
        t_celldata.MinDim=10;
        
        mshBulkPhyGroupNum=0;
        mshPhyGroupNum=0;
        mshPhyGroupDimVec.clear();
        mshPhyGroupNameVec.clear();
        mshPhyGroupIDVec.clear();
        
        maxPhyGroupID=-1;
        maxPhyGroupDim=-1;
        minPhyGroupDim=100;
        
        mshNodalPhyGroupNum=0;
        
        t_celldata.LineElmtsNum=0;
        t_celldata.SurfElmtsNum=0;
        t_celldata.BulkElmtsNum=0;
        t_celldata.ElmtsNum=0;
        while(!in.eof()){
            getline(in,str);
            if(str.find("$MeshFormat")!=string::npos){
                // read the version, format, size
                int mysize;
                in>>version>>format>>mysize;
                if(version<4.0 || version>4.2){
                    MessagePrinter::printErrorTxt("version="+to_string(version)+" is not supported for msh4 file importer, "
                                                  "please check your mesh file");
                    return false;
                }
            }
            else if(str.find("$PhysicalNames")!=string::npos){
                int phydim,phyid;
                string phyname;
                in>>mshPhyGroupNum;
                getline(in,str);//remove \n in this line
                for(int i=0;i<mshPhyGroupNum;i++){
                    getline(in,str);
                    istringstream s_stream(str);
                    s_stream>>phydim>>phyid>>phyname;
                    // remove the '"' ,keep only the text
                    phyname.erase(remove(phyname.begin(),phyname.end(),'"'),phyname.end());
                    if(phydim>maxPhyGroupDim) maxPhyGroupDim=phydim;
                    if(minPhyGroupDim<phydim) minPhyGroupDim=phydim;
                    if(phyid>maxPhyGroupID) maxPhyGroupID=phyid;
                    
                    if(phydim==0){
                        // for nodal physical information
                        mshNodalPhyGroupNum+=1;
                        mshNodalPhyGroupNameVec.push_back(phyname);
                        mshNodalPhyGroupIDVec.push_back(phyid);
                    }
                    else{
                        mshPhyGroupDimVec.push_back(phydim);
                        mshPhyGroupIDVec.push_back(phyid);
                        mshPhyGroupNameVec.push_back(phyname);
                    }
                    if(phydim==mshMaxDim){
                        mshBulkPhyGroupNum+=1;
                    }
                    mshPhyID2NameMap[phyid]=phyname;
                }
                getline(in,str);
            }// end-of-physical-group-reading
            else if(str.find("$Entities")!=string::npos){
                getline(in,str);
                numbers=StringUtils::splitStrNum(str);
                int numPoints,numCurves,numSurfaces,numVolumes;
                if(numbers.size()!=4){
                    MessagePrinter::printErrorTxt("Invalid entity block information in your msh4 file inside the $Entities");
                    MessagePrinter::exitAsFem();
                }
                
                numPoints=static_cast<int>(numbers[1-1]);
                numCurves=static_cast<int>(numbers[2-1]);
                numSurfaces=static_cast<int>(numbers[3-1]);
                numVolumes=static_cast<int>(numbers[4-1]);
                
                // factor=20 to avoid the unordered entitie index(especially the nodal entities)
                m_PointsEntityPhyIDs.resize(numPoints*100,0);
                m_CurvesEntityPhyIDS.resize(numCurves*100,0);
                m_SurfaceEntityPhyIDs.resize(numSurfaces*100,0);
                m_VolumesEntityPhyIDs.resize(numVolumes*100,0);
                
                int i,j;
                int nodeid,curveid,surfaceid,volumeid;
                //*** read points entities
                for(i=0;i<numPoints;i++){
                    getline(in,str);
                    numbers=StringUtils::splitStrNum(str);
                    nodeid=static_cast<int>(numbers[1-1]);
                    if(static_cast<int>(numbers[5-1])>0){
                        m_PointsEntityPhyIDs[nodeid-1]=static_cast<int>(numbers[6-1]);
                    }
                }
                //*** read curvess entities
                for(i=0;i<numCurves;i++){
                    getline(in,str);
                    numbers=StringUtils::splitStrNum(str);
                    curveid=static_cast<int>(numbers[1-1]);
                    j=static_cast<int>(numbers[8-1]);
                    if(j>0){
                        m_CurvesEntityPhyIDS[curveid-1]=static_cast<int>(numbers[8+j-1]);
                    }
                }
                //*** read surfaces entities
                for(i=0;i<numSurfaces;i++){
                    getline(in,str);
                    numbers=StringUtils::splitStrNum(str);
                    surfaceid=static_cast<int>(numbers[1-1]);
                    j=static_cast<int>(numbers[8-1]);
                    if(j>0){
                        m_SurfaceEntityPhyIDs[surfaceid-1]=static_cast<int>(numbers[8+j-1]);
                    }
                }
                //*** read volumes entities
                for(i=0;i<numVolumes;i++){
                    getline(in,str);
                    numbers=StringUtils::splitStrNum(str);
                    volumeid=static_cast<int>(numbers[1-1]);
                    j=static_cast<int>(numbers[8-1]);
                    if(j>0){
                        m_VolumesEntityPhyIDs[volumeid-1]=static_cast<int>(numbers[8+j-1]);
                    }
                }
            }
            else if(str.find("$Nodes")!=string::npos){
                // read the nodes' coordinates
                // node-id, x, y, z
                vector<int> nodeid;
                
                t_celldata.NodesNum=0;
                getline(in,str); // read the entities
                numbers=StringUtils::splitStrNum(str);
                
                numEntityBlocks=static_cast<int>(numbers[1-1]);
                numNodes=static_cast<int>(numbers[2-1]);
                minNodeTag=static_cast<int>(numbers[3-1]);
                maxNodeTag=static_cast<int>(numbers[4-1]);
                
                NodeCoords.resize(3*maxNodeTag,0.0);// here the node id may not be contineous case !!!
                NodeIDFlag.resize(maxNodeTag,0);

                t_celldata.NodesNum=numNodes;
                
                double x,y,z;
                int i,j,numNodesInBlock;

                t_celldata.Xmin=t_celldata.Ymin=t_celldata.Zmin= 1.0e16;
                t_celldata.Xmax=t_celldata.Ymax=t_celldata.Zmax=-1.0e16;
                
                int nNodes=0;
                for(int nBlock=0;nBlock<numEntityBlocks;nBlock++){
                    nodeid.clear();
                    getline(in,str);
                    numbers=StringUtils::splitStrNum(str);
                    if(numbers.size()!=4){
                        MessagePrinter::printErrorTxt("Invalid node entities in your msh4 file inside the $Nodes block");
                        MessagePrinter::exitAsFem();
                    }
                    numNodesInBlock=static_cast<int>(numbers[4-1]);
                    for(i=0;i<numNodesInBlock;i++){
                        getline(in,str); // read the node id
                        if(str.size()>1){
                            if(str.at(str.size()-1)<'0' || str.at(str.size()-1)>'9'){
                                // Remove the last invalid char, this may come from the msh file generated in different platform.
                                // In some cases, the last char in this line is not '\n' or ' ', but other unrecognized character.
                                str.pop_back();
                            }
                        }
                        numbers=StringUtils::splitStrNum(str);
                        if(numbers.size()<1){
                            MessagePrinter::printErrorTxt("Can't find a valid node Tag in your msh4 file inside the $Nodes");
                            MessagePrinter::exitAsFem();
                        }
                        if(static_cast<int>(numbers[0])<minNodeTag||static_cast<int>(numbers[0])>maxNodeTag){
                            MessagePrinter::printErrorTxt("Invalid node Tag in your msh4 file inside the $Nodes");
                            MessagePrinter::exitAsFem();
                        }
                        else{
                            nodeid.push_back(static_cast<int>(numbers[0]));// store the node id
                        }
                    }
                    for(i=0;i<numNodesInBlock;i++){
                        getline(in,str);
                        numbers=StringUtils::splitStrNum(str);
                        if(numbers.size()!=3){
                            MessagePrinter::printErrorTxt("Invalid node coordinates information in your msh4 file inside the $Nodes block");
                            MessagePrinter::exitAsFem();
                        }
                        x=numbers[0];y=numbers[1];z=numbers[2];
                        j=nodeid[i];
                        NodeCoords[(j-1)*3+0]=x;
                        NodeCoords[(j-1)*3+1]=y;
                        NodeCoords[(j-1)*3+2]=z;
                        nNodes+=1;
                        NodeIDFlag[j-1]=1;

                        if(x>t_celldata.Xmax) t_celldata.Xmax=x;
                        if(x<t_celldata.Xmin) t_celldata.Xmin=x;
                        if(y>t_celldata.Ymax) t_celldata.Ymax=y;
                        if(y<t_celldata.Ymin) t_celldata.Ymin=y;
                        if(z>t_celldata.Zmax) t_celldata.Zmax=z;
                        if(z<t_celldata.Zmin) t_celldata.Zmin=z;
                    } // end-of-nodes-block-reading
                } // end-of-block-reading
                
                if(nNodes!=numNodes){
                    MessagePrinter::printErrorTxt("Something is wrong in your msh4 file inside the $Nodes block, nodes numer is not match with the first line");
                    MessagePrinter::exitAsFem();
                }
            }// end-of-node-coordinates-reading
            else if(str.find("$Elements")!=string::npos){
                vector<int> tempconn;
                int elmtid,nodeid,phyid,entityTag,elmttype,vtktype;
                int nodes,entityDim,elmtorder;
                string meshtypename;
                MeshType meshtype;

                SingleMeshCell TempCell;
                
                numEntityBlocks=0;
                numElements=0;
                minElementTag=0;
                maxElementTag=0;
                numElementsInBlock=0;

                t_celldata.LineElmtMeshType=MeshType::EDGE2;
                t_celldata.SurfElmtMeshType=MeshType::TRI3;

                PointElmtsNum=0;
                t_celldata.LineElmtsNum=0;
                t_celldata.SurfElmtsNum=0;
                t_celldata.BulkElmtsNum=0;
                t_celldata.ElmtsNum=0;
                
                getline(in,str);// read the element total entities
                numbers=StringUtils::splitStrNum(str);
                if(numbers.size()!=4){
                    //numEntityBlocks(size_t) numElements(size_t) minElementTag(size_t) maxElementTag(size_t)
                    MessagePrinter::printErrorTxt("Invalid element entities in your msh4 file inside the $Elements");
                    MessagePrinter::exitAsFem();
                }
                
                numEntityBlocks=static_cast<int>(numbers[1-1]);
                numElements=static_cast<int>(numbers[2-1]);
                minElementTag=static_cast<int>(numbers[3-1]);
                maxElementTag=static_cast<int>(numbers[4-1]);

                t_celldata.ElmtsNum=numElements;
                ElmtPhyIDVec.resize(maxElementTag,0);
                ElmtDimVec.resize(maxElementTag,0);
                ElmtConn.resize(maxElementTag);
                ElmtIDFlag.resize(maxElementTag,0);
                
                mshBulkElmtUniquePhyIDVec.clear();

                t_celldata.MinDim=10;

                TempCellVec.clear();
                BulkCellVec.clear();
                
                for(int block=0;block<numEntityBlocks;block++){
                    getline(in,str);//get current element entities
                    //entityDim(int) entityTag(int) elementType(int; see below) numElementsInBlock(size_t)
                    numbers=StringUtils::splitStrNum(str);
                    if(numbers.size()!=4){
                        MessagePrinter::printErrorTxt("Invalid element entities for current element in your msh4 file inside the $Elements");
                        MessagePrinter::exitAsFem();
                    }
                    
                    entityDim=static_cast<int>(numbers[1-1]);
                    entityTag=static_cast<int>(numbers[2-1]);
                    elmttype=static_cast<int>(numbers[3-1]);
                    numElementsInBlock=static_cast<int>(numbers[4-1]);
                    phyid=getPhysicalIDViaEntityTag(entityDim,entityTag);
                    vtktype=MshFileUtils::getElmtVTKCellTypeFromElmtType(elmttype);
                    meshtype=MshFileUtils::getElmtMeshTypeFromElmtType(elmttype);
                    elmtorder=MshFileUtils::getElmtOrderFromElmtType(elmttype);
                    nodes=MshFileUtils::getElmtNodesNumFromElmtType(elmttype);
                    meshtypename=MshFileUtils::getElmtMeshTypeNameFromElmtType(elmttype);
                    
                    for(int i=0;i<numElementsInBlock;i++){
                        getline(in,str);
                        numbers=StringUtils::splitStrNum(str);
                        nodes=static_cast<int>(numbers.size()-1);
                        elmtid=static_cast<int>(numbers[1-1]);
                        if(elmtid<minElementTag||elmtid>maxElementTag){
                            MessagePrinter::printErrorTxt("Invalid element Tag in your msh4 file inside the $Elements");
                            MessagePrinter::exitAsFem();
                        }
                        
                        ElmtIDFlag[elmtid-1]=1;
                        ElmtDimVec[elmtid-1]=entityDim;
                        ElmtPhyIDVec[elmtid-1]=phyid;
                        tempconn.resize(nodes);
                        for(int j=1;j<=nodes;j++){
                            tempconn[j-1]=static_cast<int>(numbers[j]);
                        }
                        
                        MshFileUtils::reorderNodesIndex(elmttype,tempconn);
                        ElmtConn[elmtid-1]=tempconn;

                        if(entityDim<t_celldata.MinDim) t_celldata.MinDim=entityDim;
                        
                        if(entityDim==0 && entityDim<mshMaxDim){
                            PointElmtsNum+=1;

                            TempCell.Dim=0;
                            TempCell.NodesNumPerElmt=nodes;
                            TempCell.VTKCellType=vtktype;

                            TempCell.ElmtConn.clear();
                            TempCell.ElmtNodeCoords0.resize(nodes);
                            for(int k=0;k<nodes;k++){
                                nodeid=tempconn[k];
                                TempCell.ElmtConn.push_back(nodeid);
                                TempCell.ElmtNodeCoords0(k+1,1)=NodeCoords[(nodeid-1)*3+1-1];
                                TempCell.ElmtNodeCoords0(k+1,2)=NodeCoords[(nodeid-1)*3+2-1];
                                TempCell.ElmtNodeCoords0(k+1,3)=NodeCoords[(nodeid-1)*3+3-1];
                            }
                            TempCell.ElmtNodeCoords=TempCell.ElmtNodeCoords0;

                            TempCellVec.push_back(TempCell);
                        }
                        if(entityDim==1 && entityDim<mshMaxDim){
                            t_celldata.LineElmtsNum+=1;
                            t_celldata.NodesNumPerLineElmt=nodes;
                            t_celldata.LineElmtMeshType=meshtype;

                            TempCell.Dim=1;
                            TempCell.NodesNumPerElmt=nodes;
                            TempCell.VTKCellType=vtktype;

                            TempCell.ElmtConn.clear();
                            TempCell.ElmtNodeCoords0.resize(nodes);
                            for(int k=0;k<nodes;k++){
                                nodeid=tempconn[k];
                                TempCell.ElmtConn.push_back(nodeid);
                                TempCell.ElmtNodeCoords0(k+1,1)=NodeCoords[(nodeid-1)*3+1-1];
                                TempCell.ElmtNodeCoords0(k+1,2)=NodeCoords[(nodeid-1)*3+2-1];
                                TempCell.ElmtNodeCoords0(k+1,3)=NodeCoords[(nodeid-1)*3+3-1];
                            }
                            TempCell.ElmtNodeCoords=TempCell.ElmtNodeCoords0;

                            TempCellVec.push_back(TempCell);
                        }
                        
                        if(entityDim==2 && entityDim<mshMaxDim){
                            t_celldata.SurfElmtsNum+=1;
                            t_celldata.NodesNumPerSurfElmt=nodes;
                            t_celldata.SurfElmtMeshType=meshtype;
                            //
                            TempCell.Dim=2;
                            TempCell.NodesNumPerElmt=nodes;
                            TempCell.VTKCellType=vtktype;

                            TempCell.ElmtConn.clear();
                            TempCell.ElmtNodeCoords0.resize(nodes);
                            for(int k=0;k<nodes;k++){
                                nodeid=tempconn[k];
                                TempCell.ElmtConn.push_back(nodeid);
                                TempCell.ElmtNodeCoords0(k+1,1)=NodeCoords[(nodeid-1)*3+1-1];
                                TempCell.ElmtNodeCoords0(k+1,2)=NodeCoords[(nodeid-1)*3+2-1];
                                TempCell.ElmtNodeCoords0(k+1,3)=NodeCoords[(nodeid-1)*3+3-1];
                            }
                            TempCell.ElmtNodeCoords=TempCell.ElmtNodeCoords0;

                            TempCellVec.push_back(TempCell);
                        }
                        if(entityDim==mshMaxDim){
                            t_celldata.BulkElmtsNum+=1;
                            t_celldata.MeshOrder=elmtorder;
                            t_celldata.BulkElmtMeshType=meshtype;
                            t_celldata.BulkElmtVTKCellType=vtktype;
                            t_celldata.BulkMeshTypeName=meshtypename;
                            t_celldata.BulkMeshTypeName=meshtypename;
                            t_celldata.NodesNumPerBulkElmt=nodes;
                            
                            mshBulkElmtUniquePhyIDVec.push_back(phyid);

                            TempCell.Dim=mshMaxDim;
                            TempCell.NodesNumPerElmt=nodes;
                            TempCell.VTKCellType=vtktype;
                            TempCell.CellMeshType=meshtype;

                            TempCell.ElmtConn.clear();
                            TempCell.ElmtNodeCoords0.resize(nodes);
                            for(int k=0;k<nodes;k++){
                                nodeid=tempconn[k];
                                TempCell.ElmtConn.push_back(nodeid);
                                TempCell.ElmtNodeCoords0(k+1,1)=NodeCoords[(nodeid-1)*3+1-1];
                                TempCell.ElmtNodeCoords0(k+1,2)=NodeCoords[(nodeid-1)*3+2-1];
                                TempCell.ElmtNodeCoords0(k+1,3)=NodeCoords[(nodeid-1)*3+3-1];
                            }
                            TempCell.ElmtNodeCoords=TempCell.ElmtNodeCoords0;

                            TempCellVec.push_back(TempCell);
                            BulkCellVec.push_back(TempCell);
                        }
                    } // end-of-element-in-block-reading
                } // end-of-element-block-loop
                
                // before we jump out, we should check the consistency between different elements
                if(PointElmtsNum
                  +t_celldata.LineElmtsNum
                  +t_celldata.SurfElmtsNum
                  +t_celldata.BulkElmtsNum!=numElements){
                    MessagePrinter::printErrorTxt("The elements number dosen\'t match with the total one, please check your msh(2) file");
                    return false;
                }
            } // end-of-element-reading
        }// end-of-msh-file-reading
        in.close();

        //**********************************************************************************
        //*** now we re-arrange the node id and element id, to make them to be continue
        //**********************************************************************************
        int count=0;
        t_celldata.NodeCoords_Global.resize(numNodes*3,0.0);
        count=0;
        for(int i=0;i<maxNodeTag;i++){
            if(NodeIDFlag[i]>0){
                count+=1;
                NodeIDFlag[i]=count;// the active node id
                t_celldata.NodeCoords_Global[(count-1)*3+1-1]=NodeCoords[i*3+0];
                t_celldata.NodeCoords_Global[(count-1)*3+2-1]=NodeCoords[i*3+1];
                t_celldata.NodeCoords_Global[(count-1)*3+3-1]=NodeCoords[i*3+2];
            }
        }
        if(count!=numNodes){
            MessagePrinter::printErrorTxt("Nodes number dose not match with your total nodes in msh file");
            MessagePrinter::exitAsFem();
        }
        NodeCoords.clear();
        
        int nodeid=0,j,maxnodeid;
        count=0;maxnodeid=-1;
        for(int e=0;e<maxElementTag;e++){
            if(ElmtIDFlag[e]>0){
                count+=1;
                ElmtIDFlag[e]=count;// the active element id
                for(int i=1;i<=static_cast<int>(ElmtConn[e].size());i++){
                    j=ElmtConn[e][i-1];
                    if(NodeIDFlag[j-1]){
                        nodeid=NodeIDFlag[j-1];
                        ElmtConn[e][i-1]=nodeid;
                        if(nodeid>maxnodeid) maxnodeid=nodeid;
                    }
                    else{
                        MessagePrinter::printErrorTxt("Invalid node id in "+to_string(e+1)+"-th element "+to_string(i)+"-th node, please check your msh4 file");
                        MessagePrinter::exitAsFem();
                    }
                }
            }
        }
        if(maxnodeid<=numNodes){
            t_celldata.NodesNum=maxnodeid;
        }
        else if(maxnodeid>numNodes){
            MessagePrinter::printErrorTxt("The maximum node id used by your msh4 element is larger than your total nodes number,"
                                          " please check your msh4 file");
                                          MessagePrinter::exitAsFem();
        }
        
        vector<int> ElmtDimVecCopy,ElmtPhyIDVecCopy;
        vector<vector<int>> ElmtConnCopy;
        ElmtDimVecCopy.resize(numElements,0);
        ElmtPhyIDVecCopy.resize(numElements,0);
        ElmtConnCopy.resize(numElements);
        int elmtid;
        for(int e=0;e<maxElementTag;e++){
            if(ElmtIDFlag[e]){
                // if current element is active
                elmtid=ElmtIDFlag[e];
                ElmtDimVecCopy[elmtid-1]=ElmtDimVec[e];
                ElmtPhyIDVecCopy[elmtid-1]=ElmtPhyIDVec[e];
                ElmtConnCopy[elmtid-1]=ElmtConn[e];
            }
        }
        ElmtDimVec=ElmtDimVecCopy;
        ElmtPhyIDVec=ElmtPhyIDVecCopy;
        ElmtConn=ElmtConnCopy;

        //save mesh conn into global fe cell vector
        t_celldata.MeshCell_Total.clear();
        SingleMeshCell TempCell;
        int elmts,phyid;
        string phyname;
        elmts=0;
        for(int e=0;e<numElements;e++){
            if(ElmtDimVec[e]==mshMaxDim){
                TempCell.Dim=mshMaxDim;
                TempCell.NodesNumPerElmt=BulkCellVec[elmts].NodesNumPerElmt;
                TempCell.VTKCellType=BulkCellVec[elmts].VTKCellType;
                TempCell.CellMeshType=BulkCellVec[elmts].CellMeshType;

                TempCell.ElmtConn=ElmtConn[e];
                TempCell.ElmtNodeCoords0.resize(BulkCellVec[elmts].NodesNumPerElmt);
                for(int k=0;k<BulkCellVec[elmts].NodesNumPerElmt;k++){
                    nodeid=TempCell.ElmtConn[k];
                    TempCell.ElmtNodeCoords0(k+1,1)=t_celldata.NodeCoords_Global[(nodeid-1)*3+1-1];
                    TempCell.ElmtNodeCoords0(k+1,2)=t_celldata.NodeCoords_Global[(nodeid-1)*3+2-1];
                    TempCell.ElmtNodeCoords0(k+1,3)=t_celldata.NodeCoords_Global[(nodeid-1)*3+3-1];
                }
                TempCell.ElmtNodeCoords=TempCell.ElmtNodeCoords0;
                t_celldata.MeshCell_Total.push_back(TempCell);
                phyid=ElmtPhyIDVec[e];phyname=mshPhyID2NameMap[phyid];
                elmts+=1;
            }
        }
        if(elmts!=t_celldata.BulkElmtsNum){
            MessagePrinter::printErrorTxt("The bulk mesh num dosent equal to the elements num of the fecell class");
            MessagePrinter::exitAsFem();
        }
        
        ElmtDimVecCopy.clear();
        ElmtPhyIDVecCopy.clear();
        ElmtConnCopy.clear();
        
        mshPhyGroupNum=static_cast<int>(mshPhyGroupDimVec.size());// remove the nodal phy info
        if(mshBulkPhyGroupNum==0){
            // if no any physical volume (for bulk elmt) are given, we manually add them into the physical group
            sort(mshBulkElmtUniquePhyIDVec.begin(),mshBulkElmtUniquePhyIDVec.end());
            mshBulkElmtUniquePhyIDVec.erase(unique(mshBulkElmtUniquePhyIDVec.begin(),mshBulkElmtUniquePhyIDVec.end()),mshBulkElmtUniquePhyIDVec.end());
            if(mshPhyGroupNum==0){
                // if no any phy group is defined for the lower dimension mesh
                if(mshBulkElmtUniquePhyIDVec.size()<1){
                    MessagePrinter::printErrorTxt("Invalid msh2 mesh file, no any physical id is assigned to the volume mesh, please check your mesh file");
                    return false;
                }
                // now we prepare and create the 'new' physical group info
                t_celldata.PhyGroupNum_Global=static_cast<int>(mshBulkElmtUniquePhyIDVec.size())+1;
                
                t_celldata.PhyDimVector_Global.clear();
                t_celldata.PhyIDVector_Global.clear();
                t_celldata.PhyNameVector_Global.clear();

                t_celldata.PhyName2IDMap_Global.clear();
                t_celldata.PhyID2NameMap_Global.clear();
                t_celldata.PhyGroupElmtsNumVector_Global.clear();
                t_celldata.PhyName2MeshCellVectorMap_Global.clear();
                
                int maxphyid,dim;
                maxphyid=-1;
                for(int i=0;i<t_celldata.PhyGroupNum_Global-1;i++){
                    phyid=mshBulkElmtUniquePhyIDVec[i];
                    if(phyid>maxphyid) maxphyid=phyid;
                    phyname=to_string(phyid);

                    t_celldata.PhyDimVector_Global.push_back(mshMaxDim);
                    t_celldata.PhyIDVector_Global.push_back(phyid);
                    t_celldata.PhyNameVector_Global.push_back(phyname);

                    t_celldata.PhyName2IDMap_Global[phyname]=phyid;
                    t_celldata.PhyID2NameMap_Global[phyid]=phyname;
                }
                // for all the domain
                t_celldata.PhyDimVector_Global.push_back(mshMaxDim);
                t_celldata.PhyIDVector_Global.push_back(maxphyid+1);
                t_celldata.PhyNameVector_Global.push_back("alldomain");

                t_celldata.PhyName2IDMap_Global["alldomain"]=maxphyid+1;
                t_celldata.PhyID2NameMap_Global[maxphyid+1]="alldomain";

                // now we need to loop all the elements and clasify all the element sets
                t_celldata.PhyName2MeshCellVectorMap_Global.clear();
                t_celldata.PhyID2MeshCellVectorMap_Global.clear();
                //
                t_celldata.PhyID2BulkFECellIDMap_Global.clear();
                t_celldata.PhyName2BulkFECellIDMap_Global.clear();
                elmts=0;
                for(int e=0;e<t_celldata.ElmtsNum;e++){
                    phyid=ElmtPhyIDVec[e];
                    dim=ElmtDimVec[e];
                    if(dim==mshMaxDim){
                        // only account for the volume mesh
                        for(int i=0;i<t_celldata.PhyGroupNum_Global-1;i++){
                            if(phyid==t_celldata.PhyIDVector_Global[i]){
                                t_celldata.PhyGroupElmtsNumVector_Global[i]+=1;
                                phyname=t_celldata.PhyNameVector_Global[i];

                                t_celldata.PhyName2MeshCellVectorMap_Global[phyname].push_back(TempCellVec[e]);
                                t_celldata.PhyID2MeshCellVectorMap_Global[phyid].push_back(TempCellVec[e]);

                            }
                        }
                        t_celldata.PhyName2MeshCellVectorMap_Global["alldomain"].push_back(TempCellVec[e]);
                        t_celldata.PhyID2MeshCellVectorMap_Global[maxphyid+1].push_back(TempCellVec[e]);

                        elmts+=1;
                        phyname=t_celldata.PhyID2NameMap_Global[phyid];
                        t_celldata.PhyID2BulkFECellIDMap_Global[phyid].push_back(elmts);
                        t_celldata.PhyName2BulkFECellIDMap_Global[phyname].push_back(elmts);
                        //
                        t_celldata.PhyID2BulkFECellIDMap_Global[maxphyid+1].push_back(elmts);
                        t_celldata.PhyName2BulkFECellIDMap_Global["alldomain"].push_back(elmts);
                    }
                }
                if(elmts!=t_celldata.BulkElmtsNum){
                    MessagePrinter::printErrorTxt("The bulk mesh num dosent equal to the elements num of the fecell class");
                    MessagePrinter::exitAsFem();
                }
                t_celldata.PhyGroupElmtsNumVector_Global[t_celldata.PhyGroupNum_Global-1]=t_celldata.BulkElmtsNum;
            }
            else{
                // if we have the lower dimension physical group but zero volume mesh physical info, we should 
                // add them into the group 
                mshBulkPhyGroupNum=static_cast<int>(mshBulkElmtUniquePhyIDVec.size());
                t_celldata.PhyGroupNum_Global=mshPhyGroupNum+mshBulkPhyGroupNum+1;

                t_celldata.PhyDimVector_Global.clear();
                t_celldata.PhyIDVector_Global.clear();
                t_celldata.PhyNameVector_Global.clear();

                t_celldata.PhyName2IDMap_Global.clear();
                t_celldata.PhyID2NameMap_Global.clear();

                t_celldata.PhyGroupElmtsNumVector_Global.clear();

                t_celldata.PhyName2MeshCellVectorMap_Global.clear();
                t_celldata.PhyID2MeshCellVectorMap_Global.clear();

                t_celldata.PhyID2BulkFECellIDMap_Global.clear();
                t_celldata.PhyName2BulkFECellIDMap_Global.clear();

                // for the predefined phy group
                for(int i=0;i<mshPhyGroupNum;i++){
                    t_celldata.PhyDimVector_Global.push_back(mshPhyGroupDimVec[i]);
                    t_celldata.PhyIDVector_Global.push_back(mshPhyGroupIDVec[i]);
                    t_celldata.PhyNameVector_Global.push_back(mshPhyGroupNameVec[i]);

                    t_celldata.PhyName2IDMap_Global[mshPhyGroupNameVec[i]]=mshPhyGroupIDVec[i];
                    t_celldata.PhyID2NameMap_Global[mshPhyGroupIDVec[i]]=mshPhyGroupNameVec[i];
                    
                    if(mshPhyGroupDimVec[i]==t_celldata.MaxDim){
                        MessagePrinter::printErrorTxt("Invalid info, your phy group dim="+to_string(t_celldata.MaxDim)+"(name="+mshPhyGroupNameVec[i]+"), however, it is not a volume mesh, please check your mesh file");
                        return false;
                    }
                }
                // now we add the volume physical info
                int maxphyid,dim;
                maxphyid=-1;
                for(int i=0;i<mshBulkPhyGroupNum;i++){
                    phyid=mshBulkElmtUniquePhyIDVec[i];
                    phyname=to_string(phyid);
                    if(phyid>maxphyid) maxphyid=phyid;

                    t_celldata.PhyDimVector_Global.push_back(mshMaxDim);
                    t_celldata.PhyIDVector_Global.push_back(phyid);
                    t_celldata.PhyNameVector_Global.push_back(phyname);

                    t_celldata.PhyID2NameMap_Global[phyid]=phyname;
                    t_celldata.PhyName2IDMap_Global[phyname]=phyid;
                }
                // for all domain
                phyid=maxphyid+1;
                phyname="alldomain";

                t_celldata.PhyDimVector_Global.push_back(mshMaxDim);
                t_celldata.PhyIDVector_Global.push_back(phyid);
                t_celldata.PhyNameVector_Global.push_back(phyname);

                t_celldata.PhyID2NameMap_Global[phyid]=phyname;
                t_celldata.PhyName2IDMap_Global[phyname]=phyid;

                t_celldata.PhyGroupElmtsNumVector_Global.resize(t_celldata.PhyGroupNum_Global,0);

                elmts=0;
                for(int e=0;e<t_celldata.ElmtsNum;e++){
                    phyid=ElmtPhyIDVec[e];
                    dim=ElmtDimVec[e];
                    for(int i=0;i<t_celldata.PhyGroupNum_Global-1;i++){
                        if(phyid==t_celldata.PhyIDVector_Global[i]){
                            phyname=t_celldata.PhyNameVector_Global[i];
                            t_celldata.PhyGroupElmtsNumVector_Global[i]+=1;

                            t_celldata.PhyName2MeshCellVectorMap_Global[phyname].push_back(TempCellVec[e]);
                            t_celldata.PhyID2MeshCellVectorMap_Global[phyid].push_back(TempCellVec[e]);
                            if(dim==mshMaxDim){
                                // for "alldomain"
                                t_celldata.PhyName2MeshCellVectorMap_Global["alldomain"].push_back(TempCellVec[e]);
                                t_celldata.PhyID2MeshCellVectorMap_Global[maxphyid+1].push_back(TempCellVec[e]);
                                //
                                elmts+=1;
                                t_celldata.PhyID2BulkFECellIDMap_Global[phyid].push_back(elmts);
                                t_celldata.PhyName2BulkFECellIDMap_Global[phyname].push_back(elmts);
                                //
                                t_celldata.PhyID2BulkFECellIDMap_Global[maxphyid+1].push_back(elmts);
                                t_celldata.PhyName2BulkFECellIDMap_Global["alldomain"].push_back(elmts);
                            }
                        }
                    }// end-of-phygroup-loop
                }// end-of-element-loop
                if(elmts!=t_celldata.BulkElmtsNum){
                    MessagePrinter::printErrorTxt("The bulk mesh num dosent equal to the elements num read from gmsh2 file, please check it");
                    MessagePrinter::exitAsFem();
                }
                // add "alldomain" info
                t_celldata.PhyGroupElmtsNumVector_Global[t_celldata.PhyGroupNum_Global-1]=t_celldata.BulkElmtsNum;
            }
        } // end-of-zero-volume-phy-info-case
        else{
            // the volume physical info is given, then we directly add them into meshdata
            t_celldata.PhyGroupNum_Global=mshPhyGroupNum+1;

            t_celldata.PhyDimVector_Global.clear();
            t_celldata.PhyIDVector_Global.clear();
            t_celldata.PhyNameVector_Global.clear();

            t_celldata.PhyName2IDMap_Global.clear();
            t_celldata.PhyID2NameMap_Global.clear();

            t_celldata.PhyName2MeshCellVectorMap_Global.clear();
            t_celldata.PhyID2MeshCellVectorMap_Global.clear();

            t_celldata.PhyGroupElmtsNumVector_Global.clear();

            t_celldata.PhyID2BulkFECellIDMap_Global.clear();
            t_celldata.PhyName2BulkFECellIDMap_Global.clear();
            
            int dim,maxphyid;
            maxphyid=-1;
            for(int i=0;i<mshPhyGroupNum;i++){
                dim=mshPhyGroupDimVec[i];
                phyid=mshPhyGroupIDVec[i];
                phyname=mshPhyGroupNameVec[i];
                
                if(phyid>maxphyid) maxphyid=phyid;

                t_celldata.PhyDimVector_Global.push_back(dim);
                t_celldata.PhyIDVector_Global.push_back(phyid);
                t_celldata.PhyNameVector_Global.push_back(phyname);

                t_celldata.PhyID2NameMap_Global[phyid]=phyname;
                t_celldata.PhyName2IDMap_Global[phyname]=phyid;
            }
            t_celldata.PhyGroupElmtsNumVector_Global.resize(t_celldata.PhyGroupNum_Global,0);

            elmts=0;
            for(int e=0;e<t_celldata.ElmtsNum;e++){
                phyid=ElmtPhyIDVec[e];
                dim=ElmtDimVec[e];
                for(int i=0;i<mshPhyGroupNum;i++){
                    if(phyid==t_celldata.PhyIDVector_Global[i]){
                        phyname=t_celldata.PhyNameVector_Global[i];
                        t_celldata.PhyGroupElmtsNumVector_Global[i]+=1;

                        t_celldata.PhyName2MeshCellVectorMap_Global[phyname].push_back(TempCellVec[e]);
                        t_celldata.PhyID2MeshCellVectorMap_Global[phyid].push_back(TempCellVec[e]);
                        //
                        if (dim==mshMaxDim) {
                            elmts+=1;
                            t_celldata.PhyID2BulkFECellIDMap_Global[phyid].push_back(elmts);
                            t_celldata.PhyName2BulkFECellIDMap_Global[phyname].push_back(elmts);
                            t_celldata.PhyID2BulkFECellIDMap_Global[maxphyid+1].push_back(elmts);
                            t_celldata.PhyName2BulkFECellIDMap_Global["alldomain"].push_back(elmts);
                        }
                    }
                }
            }
            if (elmts!=t_celldata.BulkElmtsNum) {
                MessagePrinter::printErrorTxt("Your bulk elements number dosent equal to the one read from gmsh2 file, please check your gmsh2 file");
                MessagePrinter::exitAsFem();
            }
            // now we add "alldomain" info
            dim=mshMaxDim;
            phyname="alldomain";
            phyid=maxphyid+1;

            t_celldata.PhyDimVector_Global.push_back(dim);
            t_celldata.PhyIDVector_Global.push_back(phyid);
            t_celldata.PhyNameVector_Global.push_back(phyname);

            t_celldata.PhyName2IDMap_Global[phyname]=phyid;
            t_celldata.PhyID2NameMap_Global[phyid]=phyname;

            t_celldata.PhyID2MeshCellVectorMap_Global[phyid]=BulkCellVec;
            t_celldata.PhyName2MeshCellVectorMap_Global[phyname]=BulkCellVec;

            t_celldata.PhyGroupElmtsNumVector_Global[t_celldata.PhyGroupNum_Global-1]=t_celldata.BulkElmtsNum;

        }// end-of-the volume physical info is given-if

        // for nodal physical group
        if(mshNodalPhyGroupNum){
            t_celldata.NodalPhyGroupNum_Global=mshNodalPhyGroupNum;

            t_celldata.NodalPhyIDVector_Global.resize(t_celldata.NodalPhyGroupNum_Global,0);
            t_celldata.NodalPhyNameVector_Global.resize(t_celldata.NodalPhyGroupNum_Global);

            t_celldata.NodalPhyGroupNodesNumVector_Global.resize(t_celldata.NodalPhyGroupNum_Global,0);


            int dim;
            for(int i=0;i<mshNodalPhyGroupNum;i++){
                phyid=mshNodalPhyGroupIDVec[i];
                phyname=mshNodalPhyGroupNameVec[i];

                t_celldata.NodalPhyIDVector_Global.push_back(phyid);
                t_celldata.NodalPhyNameVector_Global.push_back(phyname);

                t_celldata.NodalPhyID2NameMap_Global[phyid]=phyname;
                t_celldata.NodalPhyName2IDMap_Global[phyname]=phyid;

            }
            for(int e=0;e<numElements;e++){
                dim=ElmtDimVec[e];
                phyid=ElmtPhyIDVec[e];
                if(dim==0){
                    for(int i=0;i<mshNodalPhyGroupNum;i++){
                        if(phyid==mshNodalPhyGroupIDVec[i]){
                            phyname=mshNodalPhyGroupNameVec[i];

                            t_celldata.NodalPhyName2NodeIDVecMap_Global[phyname].push_back(TempCellVec[e].ElmtConn[0]);
                            t_celldata.NodalPhyGroupNodesNumVector_Global[i]+=1;
                        }
                    }
                }
            }
        }
        else{

            t_celldata.NodalPhyGroupNum_Global=0;

            t_celldata.NodalPhyNameVector_Global.clear();
            t_celldata.NodalPhyIDVector_Global.clear();
            t_celldata.NodalPhyGroupNodesNumVector_Global.clear();

            t_celldata.NodalPhyName2IDMap_Global.clear();
            t_celldata.NodalPhyID2NameMap_Global.clear();

            t_celldata.NodalPhyName2NodeIDVecMap_Global.clear();
        }

        t_celldata.BulkCellPartionInfo_Global.resize(t_celldata.BulkElmtsNum,0);

        // send basic physical group info to other ranks
        int tag;
        for(int cpuid=1;cpuid<size;cpuid++){
            tag=cpuid*1000+1;
            MPIDataBus::sendIntegerToOthers(t_celldata.PhyGroupNum_Global,tag,cpuid);
            tag=cpuid*1000+2;
            MPIDataBus::sentIntegerVectorToOthers(t_celldata.PhyDimVector_Global,tag,cpuid);
            tag=cpuid*1000+3;
            MPIDataBus::sentStringVectorToOthers(t_celldata.PhyNameVector_Global,tag,cpuid);
            //
            tag=cpuid*1000+4;
            MPIDataBus::sendMeshTypeToOthers(t_celldata.BulkElmtMeshType,tag,cpuid);
            tag=cpuid*1000+5;
            MPIDataBus::sendMeshTypeToOthers(t_celldata.SurfElmtMeshType,tag,cpuid);
            tag=cpuid*1000+6;
            MPIDataBus::sendMeshTypeToOthers(t_celldata.LineElmtMeshType,tag,cpuid);
        }

    } // end-of-rank-0
    else{
        // setup the basic physical group info
        int tag;
        tag=rank*1000+1;
        MPIDataBus::receiveIntegerFromMaster(t_celldata.PhyGroupNum_Global,tag);
        tag=rank*1000+2;
        MPIDataBus::receiveIntegerVectorFromMaster(t_celldata.PhyDimVector_Global,tag);
        tag=rank*1000+3;
        MPIDataBus::receiveStringVectorFromMaster(t_celldata.PhyNameVector_Global,tag);
        //
        tag=rank*1000+4;
        MPIDataBus::receiveMeshTypeFromMater(t_celldata.BulkElmtMeshType,tag);
        //
        tag=rank*1000+5;
        MPIDataBus::receiveMeshTypeFromMater(t_celldata.SurfElmtMeshType,tag);
        //
        tag=rank*1000+6;
        MPIDataBus::receiveMeshTypeFromMater(t_celldata.LineElmtMeshType,tag);
    }

    /**
     * Share some common variables among ranks
     */
    MPI_Bcast(&t_celldata.MeshOrder,1,MPI_INT,0,MPI_COMM_WORLD);

    MPI_Bcast(&t_celldata.BulkElmtsNum,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&t_celldata.LineElmtsNum,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&t_celldata.SurfElmtsNum,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&t_celldata.ElmtsNum,1,MPI_INT,0,MPI_COMM_WORLD);

    MPI_Bcast(&t_celldata.NodesNum,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&t_celldata.NodesNumPerBulkElmt,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&t_celldata.NodesNumPerSurfElmt,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&t_celldata.NodesNumPerLineElmt,1,MPI_INT,0,MPI_COMM_WORLD);

    MPI_Bcast(&t_celldata.BulkElmtVTKCellType,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&t_celldata.SurfElmtVTKCellType,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&t_celldata.LineElmtVTKCellType,1,MPI_INT,0,MPI_COMM_WORLD);

    return true;
}