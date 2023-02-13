//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2022.08.28
//+++ Purpose: Implement the msh file (version-4) import function.
//+++          This mesh file must be the *.msh in version-4.
//+++          For version-2, please use Msh2FileImporter.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Mesh/Msh4FileImporter.h"
#include "Utils/StringUtils.h"

int Msh4FileImporter::getMaxMeshDim(const string &filename)const{
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

int Msh4FileImporter::getPhysicalIDViaEntityTag(const int &entityDim,const int &entityTag)const{
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

bool Msh4FileImporter::importMeshFile(const string &filename,MeshData &meshdata){
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

    vector<int> ElmtPhyIDVec;
    vector<int> ElmtDimVec;
    vector<vector<int>> ElmtConn;
    vector<int> ElmtIDFlag;

    ifstream in;
    in.open(filename.c_str(),ios::in);
    if(!in.is_open()){
        MessagePrinter::printErrorTxt("can\'t read the .msh file("+filename+"),please make sure file name is correct"
                                      " or you have the access permission");
        return false;
    }
    string str;
    double version;
    int format,size;
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

    meshdata.m_maxdim=mshMaxDim;
    meshdata.m_mindim=10;

    mshBulkPhyGroupNum=0;
    mshPhyGroupNum=0;
    mshPhyGroupDimVec.clear();
    mshPhyGroupNameVec.clear();
    mshPhyGroupIDVec.clear();

    maxPhyGroupID=-1;
    maxPhyGroupDim=-1;
    minPhyGroupDim=100;

    mshNodalPhyGroupNum=0;

    meshdata.m_pointelmts=0;
    meshdata.m_lineelmts=0;
    meshdata.m_surfaceelmts=0;
    meshdata.m_bulkelmts=0;

    while(!in.eof()){
        getline(in,str);

        if(str.find("$MeshFormat")!=string::npos){
            // read the version, format, size
            in>>version>>format>>size;
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

            meshdata.m_nodes=0;
            getline(in,str); // read the entities
            numbers=StringUtils::splitStrNum(str);

            numEntityBlocks=static_cast<int>(numbers[1-1]);
            numNodes=static_cast<int>(numbers[2-1]);
            minNodeTag=static_cast<int>(numbers[3-1]);
            maxNodeTag=static_cast<int>(numbers[4-1]);

            NodeCoords.resize(3*maxNodeTag,0.0);// here the node id may not be contineous case !!!
            NodeIDFlag.resize(maxNodeTag,0);

            meshdata.m_nodes=numNodes;

            double x,y,z;
            int i,j,numNodesInBlock;
            
            meshdata.m_xmin=meshdata.m_ymin=meshdata.m_zmin= 1.0e16;
            meshdata.m_xmax=meshdata.m_ymax=meshdata.m_zmax=-1.0e16;

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

                    if(x>meshdata.m_xmax) meshdata.m_xmax=x;
                    if(x<meshdata.m_xmin) meshdata.m_xmin=x;
                    if(y>meshdata.m_ymax) meshdata.m_ymax=y;
                    if(y<meshdata.m_ymin) meshdata.m_ymin=y;
                    if(z>meshdata.m_zmax) meshdata.m_zmax=z;
                    if(z<meshdata.m_zmin) meshdata.m_zmin=z;
                } // end-of-nodes-block-reading
            } // end-of-block-reading

            if(nNodes!=numNodes){
                MessagePrinter::printErrorTxt("Something is wrong in your msh4 file inside the $Nodes block, nodes numer is not match with the first line");
                MessagePrinter::exitAsFem();
            }
        }// end-of-node-coordinates-reading
        else if(str.find("$Elements")!=string::npos){
            vector<int> tempconn;
            int elmtid,phyid,entityTag,elmttype,vtktype;
            int nodes,entityDim,elmtorder;
            string meshtypename;
            MeshType meshtype;

            numEntityBlocks=0;
            numElements=0;
            minElementTag=0;
            maxElementTag=0;
            numElementsInBlock=0;

            meshdata.m_pointelmt_connectivity.clear();
            meshdata.m_lineelmt_connectivity.clear();
            meshdata.m_surfaceelmt_connectivity.clear();
            meshdata.m_bulkelmt_connectivity.clear();

            meshdata.m_pointelmts=0;
            meshdata.m_lineelmts=0;
            meshdata.m_surfaceelmts=0;
            meshdata.m_bulkelmts=0;
            meshdata.m_elements=0;

            meshdata.m_lineelmt_type=MeshType::EDGE2;
            meshdata.m_surfaceelmt_type=MeshType::TRI3;

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

            meshdata.m_elements=numElements;

            ElmtPhyIDVec.resize(maxElementTag,0);
            ElmtDimVec.resize(maxElementTag,0);
            ElmtConn.resize(maxElementTag);
            ElmtIDFlag.resize(maxElementTag,0);
            mshBulkElmtUniquePhyIDVec.clear();

            meshdata.m_mindim=10;

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
                    if(entityDim<meshdata.m_mindim) meshdata.m_mindim=entityDim;

                    if(entityDim==0 && entityDim<mshMaxDim){
                        meshdata.m_pointelmts+=1;
                        meshdata.m_pointelmt_connectivity.push_back(tempconn);
                        meshdata.m_pointelmt_volume.push_back(0.0);
                    }
                    if(entityDim==1 && entityDim<mshMaxDim){
                        meshdata.m_lineelmts+=1;
                        meshdata.m_lineelmt_connectivity.push_back(tempconn);
                        meshdata.m_lineelmt_type=meshtype;
                        meshdata.m_lineelmt_volume.push_back(0.0);
                        meshdata.m_nodesperlineelmt=nodes;
                    }
                    if(entityDim==2 && entityDim<mshMaxDim){
                        meshdata.m_surfaceelmts+=1;
                        meshdata.m_surfaceelmt_connectivity.push_back(tempconn);
                        meshdata.m_surfaceelmt_type=meshtype;
                        meshdata.m_surfaceelmt_volume.push_back(0.0);
                        meshdata.m_nodespersurfaceelmt=nodes;
                    }
                    if(entityDim==mshMaxDim){
                        meshdata.m_bulkelmts+=1;
                        mshBulkElmtUniquePhyIDVec.push_back(phyid);
                        meshdata.m_bulkelmt_connectivity.push_back(tempconn);
                        meshdata.m_bulkelmt_type=meshtype;
                        meshdata.m_bulkelmt_typename=meshtypename;
                        meshdata.m_bulkelmt_volume.push_back(0.0);
                        meshdata.m_bulkelmt_vtktype=vtktype;
                        meshdata.m_order=elmtorder;
                        meshdata.m_nodesperbulkelmt=nodes;
                    }
                } // end-of-element-in-block-reading
            } // end-of-element-block-loop

            // before we jump out, we check the consistency between different elements
            if(meshdata.m_pointelmts
              +meshdata.m_lineelmts
              +meshdata.m_surfaceelmts
              +meshdata.m_bulkelmts!=meshdata.m_elements){
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
    meshdata.m_nodecoords0.resize(numNodes*3,0.0);
    meshdata.m_nodecoords.resize(numNodes*3,0.0);
    count=0;
    for(int i=0;i<maxNodeTag;i++){
        if(NodeIDFlag[i]>0){
            count+=1;
            NodeIDFlag[i]=count;// the active node id
            meshdata.m_nodecoords0[(count-1)*3+1-1]=NodeCoords[i*3+0];
            meshdata.m_nodecoords0[(count-1)*3+2-1]=NodeCoords[i*3+1];
            meshdata.m_nodecoords0[(count-1)*3+3-1]=NodeCoords[i*3+2];

            meshdata.m_nodecoords[(count-1)*3+1-1]=NodeCoords[i*3+0];
            meshdata.m_nodecoords[(count-1)*3+2-1]=NodeCoords[i*3+1];
            meshdata.m_nodecoords[(count-1)*3+3-1]=NodeCoords[i*3+2];
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
    if(maxnodeid<numNodes){
        meshdata.m_nodes=maxnodeid;
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
            meshdata.m_phygroups=static_cast<int>(mshBulkElmtUniquePhyIDVec.size())+1;
            meshdata.m_phygroup_dimvec.clear();
            meshdata.m_phygroup_phyidvec.clear();
            meshdata.m_phygroup_phynamevec.clear();

            meshdata.m_phygroup_name2dimvec.clear();
            meshdata.m_phygroup_name2phyidvec.clear();
            meshdata.m_phygroup_phyid2namevec.clear();

            meshdata.m_phygroup_elmtnumvec.clear();

            meshdata.m_phygroup_name2bulkelmtidvec.clear();
            meshdata.m_phygroup_name2elmtconnvec.clear();

            meshdata.m_phygroup_nodesnumperelmtvec.clear();

            int phyid,maxphyid,dim;
            string phyname;
            maxphyid=-1;
            for(int i=0;i<meshdata.m_phygroups-1;i++){
                phyid=mshBulkElmtUniquePhyIDVec[i];
                if(phyid>maxphyid) maxphyid=phyid;
                phyname=to_string(phyid);
                meshdata.m_phygroup_dimvec.push_back(mshMaxDim);
                meshdata.m_phygroup_phyidvec.push_back(phyid);
                meshdata.m_phygroup_phynamevec.push_back(phyname);

                meshdata.m_phygroup_name2dimvec.push_back(make_pair(phyname,mshMaxDim));
                meshdata.m_phygroup_name2phyidvec.push_back(make_pair(phyname,phyid));
                meshdata.m_phygroup_phyid2namevec.push_back(make_pair(phyid,phyname));

                meshdata.m_phygroup_nodesnumperelmtvec.push_back(meshdata.m_nodesperbulkelmt);
            }
            // for all the domain
            meshdata.m_phygroup_dimvec.push_back(mshMaxDim);
            meshdata.m_phygroup_phyidvec.push_back(maxphyid+1);
            meshdata.m_phygroup_phynamevec.push_back("alldomain");

            meshdata.m_phygroup_name2dimvec.push_back(make_pair("alldomain",mshMaxDim));
            meshdata.m_phygroup_name2phyidvec.push_back(make_pair("alldomain",maxphyid+1));
            meshdata.m_phygroup_phyid2namevec.push_back(make_pair(maxphyid+1,"alldomain"));

            meshdata.m_phygroup_nodesnumperelmtvec.push_back(meshdata.m_nodesperbulkelmt);// for alldomain

            meshdata.m_phygroup_name2bulkelmtidvec.resize(meshdata.m_phygroups);
            meshdata.m_phygroup_name2elmtconnvec.resize(meshdata.m_phygroups);
            meshdata.m_phygroup_elmtnumvec.resize(meshdata.m_phygroups);

            // now we need to loop all the elements and clasify all the element sets
            // meshdata.m_phygroup_elmtnumvec.clear();
            // meshdata.m_phygroup_name2bulkelmtidvec.clear();
            // meshdata.m_phygroup_name2elmtconnvec.clear();
            int subelmts;
            subelmts=meshdata.m_pointelmts+meshdata.m_lineelmts+meshdata.m_surfaceelmts;
            for(int e=subelmts;e<meshdata.m_elements;e++){
                phyid=ElmtPhyIDVec[e];
                dim=ElmtDimVec[e];
                if(dim==mshMaxDim){
                    // only account for the volume mesh
                    for(int i=0;i<meshdata.m_phygroups-1;i++){
                        if(phyid==meshdata.m_phygroup_phyidvec[i]){
                            meshdata.m_phygroup_elmtnumvec[i]+=1;
                            phyname=meshdata.m_phygroup_phynamevec[i];

                            meshdata.m_phygroup_name2bulkelmtidvec[i].first=phyname;
                            meshdata.m_phygroup_name2bulkelmtidvec[i].second.push_back(e+1-subelmts);

                            meshdata.m_phygroup_name2elmtconnvec[i].first=phyname;
                            meshdata.m_phygroup_name2elmtconnvec[i].second.push_back(meshdata.m_bulkelmt_connectivity[e]);
                        }
                    }
                    meshdata.m_phygroup_name2bulkelmtidvec[meshdata.m_phygroups-1].first="alldomain";
                    meshdata.m_phygroup_name2bulkelmtidvec[meshdata.m_phygroups-1].second.push_back(e+1-subelmts);
                }
            }
            meshdata.m_phygroup_elmtnumvec[meshdata.m_phygroups-1]=meshdata.m_bulkelmts;
            meshdata.m_phygroup_name2elmtconnvec[meshdata.m_phygroups-1].first="alldomain";
            meshdata.m_phygroup_name2elmtconnvec[meshdata.m_phygroups-1].second=meshdata.m_bulkelmt_connectivity;
        }
        else{
            // if we have the lower dimension physical group but zero volume mesh physical info, we should 
            // add them into the group 
            mshBulkPhyGroupNum=static_cast<int>(mshBulkElmtUniquePhyIDVec.size());
            meshdata.m_phygroups=mshPhyGroupNum+mshBulkPhyGroupNum+1;

            meshdata.m_phygroup_dimvec.clear();
            meshdata.m_phygroup_phyidvec.clear();
            meshdata.m_phygroup_phynamevec.clear();

            meshdata.m_phygroup_name2dimvec.clear();
            meshdata.m_phygroup_name2phyidvec.clear();
            meshdata.m_phygroup_phyid2namevec.clear();

            meshdata.m_phygroup_elmtnumvec.clear();

            meshdata.m_phygroup_name2bulkelmtidvec.clear();
            meshdata.m_phygroup_name2elmtconnvec.clear();

            meshdata.m_phygroup_nodesnumperelmtvec.clear();

            // for the predefined phy group
            for(int i=0;i<mshPhyGroupNum;i++){
                meshdata.m_phygroup_dimvec.push_back(mshPhyGroupDimVec[i]);
                meshdata.m_phygroup_phyidvec.push_back(mshPhyGroupIDVec[i]);
                meshdata.m_phygroup_phynamevec.push_back(mshPhyGroupNameVec[i]);

                meshdata.m_phygroup_name2dimvec.push_back(make_pair(mshPhyGroupNameVec[i],mshPhyGroupDimVec[i]));
                meshdata.m_phygroup_name2phyidvec.push_back(make_pair(mshPhyGroupNameVec[i],mshPhyGroupIDVec[i]));
                meshdata.m_phygroup_phyid2namevec.push_back(make_pair(mshPhyGroupIDVec[i],mshPhyGroupNameVec[i]));
                if(mshPhyGroupDimVec[i]==1){
                    meshdata.m_phygroup_nodesnumperelmtvec.push_back(meshdata.m_nodesperlineelmt);
                }
                else if(mshPhyGroupDimVec[i]==2){
                    meshdata.m_phygroup_nodesnumperelmtvec.push_back(meshdata.m_nodespersurfaceelmt);
                }
                else if(mshPhyGroupDimVec[i]==3){
                    MessagePrinter::printErrorTxt("Invalid info, your phy group dim=3(name="+mshPhyGroupNameVec[i]+"), however, it is not a volume mesh, please check your mesh file");
                    return false;
                }
            }
            // now we add the volume physical info
            int phyid,maxphyid,dim,subelmts;
            string phyname;
            maxphyid=-1;
            for(int i=0;i<mshBulkPhyGroupNum;i++){
                phyid=mshBulkElmtUniquePhyIDVec[i];
                phyname=to_string(phyid);
                if(phyid>maxphyid) maxphyid=phyid;
                
                meshdata.m_phygroup_dimvec.push_back(mshMaxDim);
                meshdata.m_phygroup_phyidvec.push_back(phyid);
                meshdata.m_phygroup_phynamevec.push_back(phyname);

                meshdata.m_phygroup_name2dimvec.push_back(make_pair(phyname,mshMaxDim));
                meshdata.m_phygroup_name2phyidvec.push_back(make_pair(phyname,phyid));
                meshdata.m_phygroup_phyid2namevec.push_back(make_pair(phyid,phyname));
                meshdata.m_phygroup_nodesnumperelmtvec.push_back(meshdata.m_nodesperbulkelmt);
                
            }
            // for all domain
            phyid=maxphyid+1;
            phyname="alldomain";
            meshdata.m_phygroup_dimvec.push_back(mshMaxDim);
            meshdata.m_phygroup_phyidvec.push_back(phyid);
            meshdata.m_phygroup_phynamevec.push_back(phyname);

            meshdata.m_phygroup_name2dimvec.push_back(make_pair(phyname,mshMaxDim));
            meshdata.m_phygroup_name2phyidvec.push_back(make_pair(phyname,phyid));
            meshdata.m_phygroup_phyid2namevec.push_back(make_pair(phyid,phyname));
            meshdata.m_phygroup_nodesnumperelmtvec.push_back(meshdata.m_nodesperbulkelmt);

            meshdata.m_phygroup_elmtnumvec.resize(meshdata.m_phygroups,0);
            meshdata.m_phygroup_name2bulkelmtidvec.resize(mshBulkPhyGroupNum+1);
            meshdata.m_phygroup_name2elmtconnvec.resize(meshdata.m_phygroups);

            subelmts=meshdata.m_pointelmts+meshdata.m_lineelmts+meshdata.m_surfaceelmts;
            for(int e=0;e<meshdata.m_elements;e++){
                phyid=ElmtPhyIDVec[e];
                dim=ElmtDimVec[e];
                for(int i=0;i<meshdata.m_phygroups-1;i++){
                    if(phyid==meshdata.m_phygroup_phyidvec[i]){
                        phyname=meshdata.m_phygroup_phynamevec[i];
                        meshdata.m_phygroup_elmtnumvec[i]+=1;

                        meshdata.m_phygroup_name2elmtconnvec[i].first=phyname;
                        meshdata.m_phygroup_name2elmtconnvec[i].second.push_back(ElmtConn[e]);
                        if(dim==mshMaxDim){
                            // for volume mesh
                            meshdata.m_phygroup_name2bulkelmtidvec[i-mshPhyGroupNum].first=phyname;
                            meshdata.m_phygroup_name2bulkelmtidvec[i-mshPhyGroupNum].second.push_back(e-subelmts+1);
                        }
                    }
                }// end-of-phygroup-loop
                if(dim==mshMaxDim){
                    // for "alldomain"
                    meshdata.m_phygroup_name2bulkelmtidvec[mshBulkPhyGroupNum].first=phyname;
                    meshdata.m_phygroup_name2bulkelmtidvec[mshBulkPhyGroupNum].second.push_back(e-subelmts+1);

                    meshdata.m_phygroup_name2elmtconnvec[meshdata.m_phygroups-1].first="alldomain";
                    meshdata.m_phygroup_name2elmtconnvec[meshdata.m_phygroups-1].second.push_back(ElmtConn[e]);
                }
            }// end-of-element-loop
            // add "alldomain" info
            phyname="alldomain";
            meshdata.m_phygroup_elmtnumvec[meshdata.m_phygroups-1]=meshdata.m_bulkelmts;
        }
    } // end-of-zero-volume-phy-info-case
    else{
        // the volume physical info is given, then we directly add them into meshdata
        meshdata.m_phygroups=mshPhyGroupNum+1;

        meshdata.m_phygroup_dimvec.clear();
        meshdata.m_phygroup_phyidvec.clear();
        meshdata.m_phygroup_phynamevec.clear();

        meshdata.m_phygroup_name2dimvec.clear();
        meshdata.m_phygroup_name2phyidvec.clear();
        meshdata.m_phygroup_phyid2namevec.clear();

        meshdata.m_phygroup_elmtnumvec.clear();

        meshdata.m_phygroup_name2bulkelmtidvec.clear();
        meshdata.m_phygroup_name2elmtconnvec.clear();

        meshdata.m_phygroup_nodesnumperelmtvec.clear();

        int dim,phyid,maxphyid,subelmts;
        string phyname;
        maxphyid=-1;
        subelmts=meshdata.m_pointelmts+meshdata.m_lineelmts+meshdata.m_surfaceelmts;
        for(int i=0;i<mshPhyGroupNum;i++){
            dim=mshPhyGroupDimVec[i];
            phyid=mshPhyGroupIDVec[i];
            phyname=mshPhyGroupNameVec[i];

            if(phyid>maxphyid) maxphyid=phyid;

            meshdata.m_phygroup_dimvec.push_back(dim);
            meshdata.m_phygroup_phyidvec.push_back(phyid);
            meshdata.m_phygroup_phynamevec.push_back(phyname);

            meshdata.m_phygroup_name2dimvec.push_back(make_pair(phyname,dim));
            meshdata.m_phygroup_name2phyidvec.push_back(make_pair(phyname,phyid));
            meshdata.m_phygroup_phyid2namevec.push_back(make_pair(phyid,phyname));
            if(dim==1 && dim<mshMaxDim){
                meshdata.m_phygroup_nodesnumperelmtvec.push_back(meshdata.m_nodesperlineelmt);
            }
            else if(dim==2 && dim<mshMaxDim){
                meshdata.m_phygroup_nodesnumperelmtvec.push_back(meshdata.m_nodespersurfaceelmt);
            }
            if(dim==mshMaxDim){
                meshdata.m_phygroup_nodesnumperelmtvec.push_back(meshdata.m_nodesperbulkelmt);
            }
        }
        meshdata.m_phygroup_elmtnumvec.resize(meshdata.m_phygroups,0);
        meshdata.m_phygroup_name2bulkelmtidvec.resize(mshBulkPhyGroupNum+1);
        meshdata.m_phygroup_name2elmtconnvec.resize(meshdata.m_phygroups);
        for(int e=0;e<meshdata.m_elements;e++){
            phyid=ElmtPhyIDVec[e];
            dim=ElmtDimVec[e];
            for(int i=0;i<mshPhyGroupNum;i++){
                if(phyid==meshdata.m_phygroup_phyidvec[i]){
                    phyname=meshdata.m_phygroup_phynamevec[i];
                    meshdata.m_phygroup_elmtnumvec[i]+=1;

                    meshdata.m_phygroup_name2elmtconnvec[i].first=phyname;
                    meshdata.m_phygroup_name2elmtconnvec[i].second.push_back(ElmtConn[e]);

                    if(dim==mshMaxDim){
                        meshdata.m_phygroup_name2bulkelmtidvec[i-(mshPhyGroupNum-mshBulkPhyGroupNum)].first=phyname;
                        meshdata.m_phygroup_name2bulkelmtidvec[i-(mshPhyGroupNum-mshBulkPhyGroupNum)].second.push_back(e+1-subelmts);
                    }
                }
            }
            if(dim==mshMaxDim){
                // for "alldomain"
                meshdata.m_phygroup_name2bulkelmtidvec[mshBulkPhyGroupNum].first="alldomain";
                meshdata.m_phygroup_name2bulkelmtidvec[mshBulkPhyGroupNum].second.push_back(e+1-subelmts);
            }
        }
        // now we add "alldomain" info
        dim=mshMaxDim;
        phyname="alldomain";
        phyid=maxphyid+1;
        meshdata.m_phygroup_dimvec.push_back(dim);
        meshdata.m_phygroup_phyidvec.push_back(phyid);
        meshdata.m_phygroup_phynamevec.push_back(phyname);

        meshdata.m_phygroup_name2dimvec.push_back(make_pair(phyname,dim));
        meshdata.m_phygroup_name2phyidvec.push_back(make_pair(phyname,phyid));
        meshdata.m_phygroup_phyid2namevec.push_back(make_pair(phyid,phyname));

        meshdata.m_phygroup_nodesnumperelmtvec.push_back(meshdata.m_nodesperbulkelmt);

        meshdata.m_phygroup_elmtnumvec[meshdata.m_phygroups-1]=meshdata.m_bulkelmts;

        meshdata.m_phygroup_name2elmtconnvec[meshdata.m_phygroups-1].first=phyname;
        meshdata.m_phygroup_name2elmtconnvec[meshdata.m_phygroups-1].second=meshdata.m_bulkelmt_connectivity;

    }

    // for nodal physical group
    if(mshNodalPhyGroupNum){
        meshdata.m_nodal_phygroups=mshNodalPhyGroupNum;

        meshdata.m_nodephygroup_name2nodeidvec.resize(meshdata.m_nodal_phygroups);

        string phyname;
        int phyid,dim;
        for(int i=0;i<mshNodalPhyGroupNum;i++){
            phyid=mshNodalPhyGroupIDVec[i];
            phyname=mshNodalPhyGroupNameVec[i];
            meshdata.m_nodephygroup_name2phyidvec.push_back(make_pair(phyname,phyid));
            meshdata.m_nodephygroup_phyid2namevec.push_back(make_pair(phyid,phyname));
            meshdata.m_nodephygroup_phynamevec.push_back(phyname);
            meshdata.m_nodephygroup_phyidvec.push_back(phyid);
        }
        for(int e=0;e<numElements;e++){
            dim=ElmtDimVec[e];
            phyid=ElmtPhyIDVec[e];
            if(dim==0){
                for(int i=0;i<mshNodalPhyGroupNum;i++){
                    if(phyid==mshNodalPhyGroupIDVec[i]){
                        phyname=mshNodalPhyGroupNameVec[i];
                        meshdata.m_nodephygroup_name2nodeidvec[i].first=phyname;
                        meshdata.m_nodephygroup_name2nodeidvec[i].second.push_back(ElmtConn[e][0]);
                    }
                }
            }
        }
    }
    else{
        meshdata.m_nodal_phygroups=0;
        meshdata.m_nodephygroup_name2nodeidvec.clear();
        meshdata.m_nodephygroup_name2phyidvec.clear();
        meshdata.m_nodephygroup_phyid2namevec.clear();
        meshdata.m_nodephygroup_phynamevec.clear();
        meshdata.m_nodephygroup_phyidvec.clear();
    }

    return true;
}