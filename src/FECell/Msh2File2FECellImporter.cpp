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
//+++ Date   : 2024.02.03
//+++ Purpose: Implement the msh file (version-2) import function.
//+++          This mesh file must be the *.msh in version-2.
//+++          For version-4, please use Msh4FileImporter.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FECell/Msh2File2FECellImporter.h"
#include "MPIUtils/MPIDataBus.h"

int Msh2File2FECellImporter::getMaxMeshDim(const string &filename)const{
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
            int nelmts;
            int elmtid,phyid,geoid,ntags,elmttype,dim,nodes;
            in>>nelmts;
            maxdim=-1;
            for(int e=0;e<nelmts;e++){
                in>>elmtid>>elmttype>>ntags>>phyid>>geoid;
                dim=MshFileUtils::getElmtDimFromElmtType(elmttype);
                nodes=MshFileUtils::getElmtNodesNumFromElmtType(elmttype);
                for(int i=0;i<nodes;i++) in>>geoid;
                if(dim>maxdim) maxdim=dim;
            }
            break;
        }
    }
    in.close();
    return maxdim;
}

bool Msh2File2FECellImporter::importMeshFile(const string &filename,FECellData &t_celldata){

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

    if(rank!=0) in.close();// close the ifstream on other ranks

    if(rank==0){
        int mshMaxDim,mshMinDim;
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

        map<int,string> mshPhyID2NameMap;
        map<string,int> mshPhyName2IDMap;
        map<int,int> mshPhyID2DimMap;
        map<string,int> mshPhyName2DimMap;

        map<string,vector<SingleMeshCell>> mshPhyName2CellVecMap;
        map<int,vector<SingleMeshCell>>    mshPhyID2CellVecMap;
        
        string str;
        double version;
        int format;
        
        mshMaxDim=getMaxMeshDim(filename);
        t_celldata.MaxDim=mshMaxDim;
        
        mshBulkPhyGroupNum=0;
        mshPhyGroupNum=0;
        mshPhyGroupDimVec.clear();
        mshPhyGroupNameVec.clear();
        mshPhyGroupIDVec.clear();

        mshPhyID2NameMap.clear();
        mshPhyID2DimMap.clear();
        mshPhyID2CellVecMap.clear();

        mshPhyName2IDMap.clear();
        mshPhyName2DimMap.clear();
        mshPhyName2CellVecMap.clear();
        
        maxPhyGroupID=-1;
        maxPhyGroupDim=-1;
        minPhyGroupDim=100;
        
        mshNodalPhyGroupNum=0;
        
        t_celldata.LineElmtsNum=0;
        t_celldata.SurfElmtsNum=0;
        t_celldata.BulkElmtsNum=0;
        // only read mesh from master rank
        while(!in.eof()){
            getline(in,str);

            if(str.find("$MeshFormat")!=string::npos){
                // read the version, format, size
                int mshsize;
                in>>version>>format>>mshsize;
                if(version<2.0 || version>2.2){
                    MessagePrinter::printErrorTxt("version="+to_string(version)+" is not supported for msh2 file importer, "
                                                  "please check your mesh file");
                    return false;
                }
            }
            else if(str.find("$PhysicalNames")!=string::npos){
                int phydim,phyid;
                string phyname;
                mshBulkPhyGroupNum=0;
                in>>mshPhyGroupNum;
                getline(in,str);//remove \n in this line
                for(int i=0;i<mshPhyGroupNum;i++){
                    getline(in,str);
                    istringstream s_stream(str);
                    s_stream>>phydim>>phyid>>phyname;
                    
                    // remove the '"' ,keep only the text
                    phyname.erase(remove(phyname.begin(),phyname.end(),'"'),phyname.end());

                    mshPhyID2NameMap[phyid]=phyname;
                    mshPhyID2DimMap[phyid]=phydim;
                    
                    mshPhyName2IDMap[phyname]=phyid;
                    mshPhyName2DimMap[phyname]=phyid;
                    
                    if(phydim>maxPhyGroupDim) maxPhyGroupDim=phydim;
                    if(minPhyGroupDim<phydim) minPhyGroupDim=phydim;
                    if(phyid>maxPhyGroupID)   maxPhyGroupID=phyid;
                    
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
            else if(str.find("$Nodes")!=string::npos){
                // read the nodes' coordinates
                // node-id, x, y, z
                t_celldata.NodesNum=0;
                in>>t_celldata.NodesNum;
                int nodeid;
                double x,y,z;
                t_celldata.Xmin=t_celldata.Ymin=t_celldata.Zmin= 1.0e16;
                t_celldata.Xmax=t_celldata.Ymax=t_celldata.Zmax=-1.0e16;

                t_celldata.NodeCoords_Global.resize(t_celldata.NodesNum*3,0.0);
                for(int i=0;i<t_celldata.NodesNum;i++){
                    in>>nodeid>>x>>y>>z;
                    t_celldata.NodeCoords_Global[(nodeid-1)*3+1-1]=x;
                    t_celldata.NodeCoords_Global[(nodeid-1)*3+2-1]=y;
                    t_celldata.NodeCoords_Global[(nodeid-1)*3+3-1]=z;
                    
                    if(x>t_celldata.Xmax) t_celldata.Xmax=x;
                    if(x<t_celldata.Xmin) t_celldata.Xmin=x;

                    if(y>t_celldata.Ymax) t_celldata.Ymax=y;
                    if(y<t_celldata.Ymin) t_celldata.Ymin=y;
                    
                    if(z>t_celldata.Zmax) t_celldata.Zmax=z;
                    if(z<t_celldata.Zmin) t_celldata.Zmin=z;
                }
                getline(in,str);
            }// end-of-node-coordinates-reading
            else if(str.find("$Elements")!=string::npos){
                t_celldata.ElmtsNum=0;
                in>>t_celldata.ElmtsNum;// total elements number
                vector<int> tempconn;
                int elmtid,phyid,geoid,ntags,elmttype,vtktype;
                int nodes,dim,elmtorder,nodeid;
                int PointsNum;
                string meshtypename,phyname;
                // MeshType meshtype;

                vector<SingleMeshCell> LineCell;
                vector<SingleMeshCell> SurfCell;
                vector<SingleMeshCell> BulkCell;

                vector<int> LineCellPhyID,SurfCellPhyID,BulkCellPhyID;

                SingleMeshCell TempCell;

                t_celldata.LineElmtsNum=0;
                t_celldata.SurfElmtsNum=0;
                t_celldata.BulkElmtsNum=0;

                t_celldata.LineElmtMeshType=MeshType::EDGE2;
                t_celldata.SurfElmtMeshType=MeshType::TRI3;
                
                ElmtPhyIDVec.resize(t_celldata.ElmtsNum,0);
                ElmtDimVec.resize(t_celldata.ElmtsNum,0);
                ElmtConn.resize(t_celldata.ElmtsNum);

                LineCell.clear();LineCellPhyID.clear();
                SurfCell.clear();SurfCellPhyID.clear();
                BulkCell.clear();BulkCellPhyID.clear();

                t_celldata.MeshCell_Local.clear();
                t_celldata.MeshCell_Total.clear();//stores the bulk mesh cell

                mshPhyID2CellVecMap.clear();
                mshPhyName2CellVecMap.clear();

                nodes=dim=vtktype=elmtorder=1;
                meshtypename.clear();
                
                PointsNum=0;
                for(int e=0;e<t_celldata.ElmtsNum;e++){
                    in>>elmtid>>elmttype>>ntags>>phyid>>geoid;
                    
                    nodes=MshFileUtils::getElmtNodesNumFromElmtType(elmttype);
                    dim=MshFileUtils::getElmtDimFromElmtType(elmttype);
                    vtktype=MshFileUtils::getElmtVTKCellTypeFromElmtType(elmttype);
                    // meshtype=MshFileUtils::getElmtMeshTypeFromElmtType(elmttype);
                    meshtypename=MshFileUtils::getElmtMeshTypeNameFromElmtType(elmttype);
                    elmtorder=MshFileUtils::getElmtOrderFromElmtType(elmttype);
                    
                    // read the connectivity info
                    tempconn.resize(nodes,0);
                    for(int j=0;j<nodes;j++) in>>tempconn[j];
                    
                    MshFileUtils::reorderNodesIndex(elmttype,tempconn);
                    if(dim<mshMinDim) mshMinDim=dim;
                    
                    ElmtDimVec[elmtid-1]=dim;
                    ElmtPhyIDVec[elmtid-1]=phyid;
                    ElmtConn[elmtid-1]=tempconn;
                    
                    if(dim==0){
                        // for node case
                        if(tempconn.size()!=1){
                            MessagePrinter::printErrorTxt("Invalid node set (element) in your msh(2) file, the nodal connectivity array should contain only 1 node");
                            MessagePrinter::exitAsFem();
                        }
                        NodeSetPhyID2IndexMap[phyid].push_back(tempconn[0]);
                        PointsNum+=1;
                    }
                    
                    if(dim==1 && dim<mshMaxDim){
                        TempCell.Dim=1;
                        TempCell.NodesNumPerElmt=nodes;
                        TempCell.VTKCellType=vtktype;
                        TempCell.ElmtConn.clear();
                        TempCell.ElmtNodeCoords0.resize(nodes);
                        for(int i=0;i<nodes;i++){
                            nodeid=tempconn[i];
                            TempCell.ElmtConn.push_back(nodeid);
                            TempCell.ElmtNodeCoords0(i+1,1)=t_celldata.NodeCoords_Global[(nodeid-1)*3+1-1];
                            TempCell.ElmtNodeCoords0(i+1,2)=t_celldata.NodeCoords_Global[(nodeid-1)*3+2-1];
                            TempCell.ElmtNodeCoords0(i+1,3)=t_celldata.NodeCoords_Global[(nodeid-1)*3+3-1];
                        }
                        TempCell.ElmtNodeCoords=TempCell.ElmtNodeCoords0;
                        LineCell.push_back(TempCell);
                        LineCellPhyID.push_back(phyid);
                        t_celldata.LineElmtsNum+=1;
                    }
                    if(dim==2 && dim<mshMaxDim){
                        TempCell.Dim=2;
                        TempCell.NodesNumPerElmt=nodes;
                        TempCell.VTKCellType=vtktype;
                        TempCell.ElmtConn.clear();
                        TempCell.ElmtNodeCoords0.resize(nodes);
                        for(int i=0;i<nodes;i++){
                            nodeid=tempconn[i];
                            TempCell.ElmtConn.push_back(nodeid);
                            TempCell.ElmtNodeCoords0(i+1,1)=t_celldata.NodeCoords_Global[(nodeid-1)*3+1-1];
                            TempCell.ElmtNodeCoords0(i+1,2)=t_celldata.NodeCoords_Global[(nodeid-1)*3+2-1];
                            TempCell.ElmtNodeCoords0(i+1,3)=t_celldata.NodeCoords_Global[(nodeid-1)*3+3-1];
                        }
                        TempCell.ElmtNodeCoords=TempCell.ElmtNodeCoords0;
                        SurfCell.push_back(TempCell);
                        SurfCellPhyID.push_back(phyid);
                        t_celldata.SurfElmtsNum+=1;
                    }
                    if(dim==mshMaxDim){
                        mshBulkElmtUniquePhyIDVec.push_back(phyid);
                        TempCell.Dim=dim;
                        t_celldata.MeshOrder=elmtorder;
                        TempCell.NodesNumPerElmt=nodes;
                        TempCell.VTKCellType=vtktype;
                        TempCell.ElmtConn.clear();
                        TempCell.ElmtNodeCoords.resize(nodes);
                        for(int i=0;i<nodes;i++){
                            nodeid=tempconn[i];
                            TempCell.ElmtConn.push_back(nodeid);
                            TempCell.ElmtNodeCoords0(i+1,1)=t_celldata.NodeCoords_Global[(nodeid-1)*3+1-1];
                            TempCell.ElmtNodeCoords0(i+1,2)=t_celldata.NodeCoords_Global[(nodeid-1)*3+2-1];
                            TempCell.ElmtNodeCoords0(i+1,3)=t_celldata.NodeCoords_Global[(nodeid-1)*3+3-1];
                        }
                        TempCell.ElmtNodeCoords=TempCell.ElmtNodeCoords0;
                        BulkCell.push_back(TempCell);
                        BulkCellPhyID.push_back(phyid);

                        t_celldata.MeshCell_Total.push_back(TempCell);
                        t_celldata.BulkElmtsNum+=1;
                    }

                    mshPhyID2CellVecMap[phyid].push_back(TempCell);
                    if(mshPhyID2NameMap.count(phyid)){
                        phyname=mshPhyID2NameMap[phyid];
                        mshPhyName2CellVecMap[phyname].push_back(TempCell);
                    }
                }// end-of-element-loop-for-element-reading
                
                // before we jump out, we should check the consistency between different cells
                if(PointsNum
                  +t_celldata.LineElmtsNum
                  +t_celldata.SurfElmtsNum
                  +t_celldata.BulkElmtsNum!=t_celldata.ElmtsNum){
                    MessagePrinter::printErrorTxt("The elements number dosen\'t match with the total one, please check your msh(2) file");
                    return false;
                }
            } // end-of-element-reading
        }// end-of-in-reading
        in.close();

        // setup the phyname<--->cell mapping
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

                t_celldata.PhyGroupElmtsNumVector_Global.clear();
                t_celldata.PhyID2NameMap_Global.clear();
                t_celldata.PhyName2IDMap_Global.clear();

                t_celldata.PhyName2MeshCellVectorMap_Global.clear();
                t_celldata.PhyID2MeshCellVectorMap_Global.clear();
                
                int phyid,maxphyid;
                string phyname;
                maxphyid=-1;
                t_celldata.PhyGroupElmtsNumVector_Global.resize(t_celldata.PhyGroupNum_Global,0);
                for(int i=0;i<t_celldata.PhyGroupNum_Global-1;i++){
                    phyid=mshBulkElmtUniquePhyIDVec[i];
                    if(phyid>maxphyid) maxphyid=phyid;
                    phyname=to_string(phyid);
                    
                    t_celldata.PhyDimVector_Global.push_back(mshMaxDim);
                    t_celldata.PhyIDVector_Global.push_back(phyid);
                    t_celldata.PhyNameVector_Global.push_back(phyname);
                    
                    t_celldata.PhyID2NameMap_Global[phyid]=phyname;
                    t_celldata.PhyName2IDMap_Global[phyname]=phyid;

                    if(mshPhyID2CellVecMap.count(phyid)){
                        t_celldata.PhyGroupElmtsNumVector_Global[i]+=1;
                        t_celldata.PhyName2MeshCellVectorMap_Global[phyname]=mshPhyID2CellVecMap[phyid];
                        t_celldata.PhyID2MeshCellVectorMap_Global[phyid]=mshPhyID2CellVecMap[phyid];
                    }
                }
                // for "alldomain" cell set
                t_celldata.PhyGroupElmtsNumVector_Global[t_celldata.PhyGroupNum_Global-1]=t_celldata.BulkElmtsNum;
                phyid=maxphyid;
                phyname="alldomain";
                t_celldata.PhyDimVector_Global.push_back(mshMaxDim);
                t_celldata.PhyIDVector_Global.push_back(phyid);
                t_celldata.PhyNameVector_Global.push_back(phyname);
                
                t_celldata.PhyID2NameMap_Global[phyid]=phyname;
                t_celldata.PhyName2IDMap_Global[phyname]=phyid;

                t_celldata.PhyName2MeshCellVectorMap_Global[phyname]=t_celldata.MeshCell_Total;
                t_celldata.PhyID2MeshCellVectorMap_Global[phyid]=t_celldata.MeshCell_Total;
            }
            else{
                // if we have the lower dimension physical group but zero volume mesh physical info, we should 
                // add them into the group 
                mshBulkPhyGroupNum=static_cast<int>(mshBulkElmtUniquePhyIDVec.size());
                
                t_celldata.PhyGroupNum_Global=mshPhyGroupNum+mshBulkPhyGroupNum+1;
                t_celldata.PhyDimVector_Global.clear();
                t_celldata.PhyIDVector_Global.clear();
                t_celldata.PhyNameVector_Global.clear();
                
                t_celldata.PhyID2NameMap_Global.clear();
                t_celldata.PhyName2IDMap_Global.clear();
                
                t_celldata.PhyName2MeshCellVectorMap_Global.clear();
                t_celldata.PhyID2MeshCellVectorMap_Global.clear();

                t_celldata.PhyGroupElmtsNumVector_Global.resize(t_celldata.PhyGroupNum_Global,0);

                int phyid,maxphyid,dim;
                string phyname;
                // for the predefined sub-dimensional phy group
                for(int i=0;i<mshPhyGroupNum;i++){
                    dim=mshPhyGroupDimVec[i];
                    phyid=mshPhyGroupIDVec[i];
                    phyname=mshPhyGroupNameVec[i];
                    t_celldata.PhyDimVector_Global.push_back(dim);
                    t_celldata.PhyIDVector_Global.push_back(phyid);
                    t_celldata.PhyNameVector_Global.push_back(phyname);
                    
                    t_celldata.PhyID2NameMap_Global[phyid]=phyname;
                    t_celldata.PhyName2IDMap_Global[phyname]=phyid;

                    if(mshPhyID2CellVecMap.count(phyid)){
                        t_celldata.PhyGroupElmtsNumVector_Global[i]+=static_cast<int>(mshPhyID2CellVecMap[phyid].size());
                        t_celldata.PhyID2MeshCellVectorMap_Global[phyid]=mshPhyID2CellVecMap[phyid];
                        t_celldata.PhyName2MeshCellVectorMap_Global[phyname]=mshPhyID2CellVecMap[phyid];
                    }
                    
                    if(dim>=maxPhyGroupDim){
                        MessagePrinter::printErrorTxt("Invalid info, your phy group dim="+to_string(dim)+"(name="+phyname+"), however, it is not a bulk mesh cell, please check your mesh file");
                        return false;
                    }
                }
                // now we add the volume physical info
                int iStart;
                maxphyid=-1;
                iStart=mshPhyGroupNum;
                for(int i=0;i<mshBulkPhyGroupNum;i++){
                    phyid=mshBulkElmtUniquePhyIDVec[i];
                    phyname=to_string(phyid);
                    if(phyid>maxphyid) maxphyid=phyid;
                    if(phyid>maxPhyGroupID){
                        t_celldata.PhyDimVector_Global.push_back(mshMaxDim);
                        t_celldata.PhyIDVector_Global.push_back(phyid);
                        t_celldata.PhyNameVector_Global.push_back(phyname);

                        t_celldata.PhyID2NameMap_Global[phyid]=phyname;
                        t_celldata.PhyName2IDMap_Global[phyname]=phyid;

                        if(mshPhyID2CellVecMap.count(phyid)){
                            t_celldata.PhyGroupElmtsNumVector_Global[i+iStart]+=static_cast<int>(mshPhyID2CellVecMap[phyid].size());
                            t_celldata.PhyID2MeshCellVectorMap_Global[phyid]=mshPhyID2CellVecMap[phyid];
                            t_celldata.PhyName2MeshCellVectorMap_Global[phyname]=mshPhyID2CellVecMap[phyid];
                        }

                    }
                }
                // for all domain
                phyid=maxphyid+1;
                phyname="alldomain";
                
                t_celldata.PhyDimVector_Global.push_back(mshMaxDim);
                t_celldata.PhyIDVector_Global.push_back(phyid);
                t_celldata.PhyNameVector_Global.push_back(phyname);
                
                t_celldata.PhyID2NameMap_Global[phyid]=phyname;
                t_celldata.PhyName2IDMap_Global[phyname]=phyid;

                t_celldata.PhyGroupElmtsNumVector_Global[t_celldata.PhyGroupNum_Global-1]=t_celldata.BulkElmtsNum;
                t_celldata.PhyID2MeshCellVectorMap_Global[phyid]=t_celldata.MeshCell_Total;
                t_celldata.PhyName2MeshCellVectorMap_Global[phyname]=t_celldata.MeshCell_Total;
            }// end-of-zero-volume-phy-but-nonzero-subdim-phygroup
        }
        else{
            // the volume physical info is given, then we directly add them into meshdata
            t_celldata.PhyGroupNum_Global=mshPhyGroupNum+1;

            t_celldata.PhyDimVector_Global.clear();
            t_celldata.PhyIDVector_Global.clear();
            t_celldata.PhyNameVector_Global.clear();
                
            t_celldata.PhyID2NameMap_Global.clear();
            t_celldata.PhyName2IDMap_Global.clear();
                
            t_celldata.PhyName2MeshCellVectorMap_Global.clear();
            t_celldata.PhyID2MeshCellVectorMap_Global.clear();

            t_celldata.PhyGroupElmtsNumVector_Global.resize(t_celldata.PhyGroupNum_Global,0);

            int phyid,maxphyid,dim;
            string phyname;

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

                if(mshPhyID2CellVecMap.count(phyid)){
                    t_celldata.PhyGroupElmtsNumVector_Global[i]+=static_cast<int>(mshPhyID2CellVecMap[phyid].size());
                    t_celldata.PhyID2MeshCellVectorMap_Global[phyid]=mshPhyID2CellVecMap[phyid];
                    t_celldata.PhyName2MeshCellVectorMap_Global[phyname]=mshPhyID2CellVecMap[phyid];
                }
            }
            // for all domain
            phyid=maxphyid+1;
            phyname="alldomain";
                
            t_celldata.PhyDimVector_Global.push_back(mshMaxDim);
            t_celldata.PhyIDVector_Global.push_back(phyid);
            t_celldata.PhyNameVector_Global.push_back(phyname);
                
            t_celldata.PhyID2NameMap_Global[phyid]=phyname;
            t_celldata.PhyName2IDMap_Global[phyname]=phyid;

            t_celldata.PhyGroupElmtsNumVector_Global[t_celldata.PhyGroupNum_Global-1]=t_celldata.BulkElmtsNum;
            t_celldata.PhyID2MeshCellVectorMap_Global[phyid]=t_celldata.MeshCell_Total;
            t_celldata.PhyName2MeshCellVectorMap_Global[phyname]=t_celldata.MeshCell_Total;
        } // end-of-nonzero-volume-phy-info-case

        /**
         * distribute the total physical group info to each local rank
        */
        int iStart,iEnd,ranksize;
        vector<SingleMeshCell> LocalCellVector;
        vector<int> nodeids;
        int phyid,dim,k;
        string phyname;
        for(int cpuid=0;cpuid<size;cpuid++){
            // send out the phynum, phyname, phyid info (except the "alldomain" group !!!)
            MPIDataBus::sendIntegerToOthers(t_celldata.PhyGroupNum_Global,cpuid*10000+0,cpuid);
            for(k=0;k<t_celldata.PhyGroupNum_Global-1;k++){
                phyid=t_celldata.PhyIDVector_Global[k];
                dim=t_celldata.PhyDimVector_Global[k];
                phyname=t_celldata.PhyNameVector_Global[k];

                MPIDataBus::sendIntegerToOthers( phyid,cpuid*10000+k*1000+1,cpuid);
                MPIDataBus::sendIntegerToOthers(   dim,cpuid*10000+k*1000+2,cpuid);
                MPIDataBus::sendStringToOthers(phyname,cpuid*10000+k*1000+3,cpuid);

                ranksize=static_cast<int>(t_celldata.PhyID2MeshCellVectorMap_Global[phyid].size())/size;
                iStart=cpuid*ranksize;
                iEnd=(cpuid+1)*ranksize;
                if(cpuid==size-1) iEnd=static_cast<int>(t_celldata.PhyID2MeshCellVectorMap_Global[phyid].size());
                
                LocalCellVector.clear();
                for(int e=iStart;e<iEnd;e++){
                    LocalCellVector.push_back(t_celldata.MeshCell_Total[e]);
                }
                if(cpuid==0){
                    t_celldata.PhyID2MeshCellVectorMap_Local[phyid]=LocalCellVector;
                    t_celldata.PhyName2MeshCellVectorMap_Local[phyname]=LocalCellVector;// each local rank share the same phy name as the master rank!!!
                }
                else{
                    MPIDataBus::sendPhyID2MeshCellMapToOthers(    phyid,LocalCellVector,10000*cpuid+k*1000+4   ,cpuid);
                    MPIDataBus::sendPhyName2MeshCellMapToOthers(phyname,LocalCellVector,10000*cpuid+k*1000+4+20,cpuid);
                }
            }
            // send out the "alldomain" mesh cell
            k=t_celldata.PhyGroupNum_Global-1;
            ranksize=t_celldata.BulkElmtsNum/size;
            iStart=cpuid*ranksize;
            iEnd=(cpuid+1)*ranksize;
            if(cpuid==size-1) iEnd=t_celldata.BulkElmtsNum;

            LocalCellVector.clear();
            for(int e=iStart;e<iEnd;e++){
                LocalCellVector.push_back(t_celldata.MeshCell_Total[e]);
            }
            if(cpuid==0){
                t_celldata.MeshCell_Local=LocalCellVector;
                t_celldata.PhyID2MeshCellVectorMap_Local[0]=LocalCellVector;
                t_celldata.PhyName2MeshCellVectorMap_Local["alldomain"]=t_celldata.MeshCell_Local;// each local rank share the same phy name as the master rank!!!
            }
            else{
                MPIDataBus::sendMeshCellToOthers(                       LocalCellVector,10000*cpuid+k*1000+4+40,cpuid);
                MPIDataBus::sendPhyID2MeshCellMapToOthers(            0,LocalCellVector,10000*cpuid+k*1000+4+60,cpuid);
                MPIDataBus::sendPhyName2MeshCellMapToOthers("alldomain",LocalCellVector,10000*cpuid+k*1000+4+80,cpuid);
            }
        } 

    }// end-of-master-rank-process
    else{
        // each local rank receive info from master rank
        int phyid,k,dim,maxdim;
        string phyname;
        MPIDataBus::receiveIntegerFromMaster(k,rank*10000+0);
        t_celldata.PhyGroupNum_Global=k;
        t_celldata.PhyDimVector_Global.resize(k,0);
        t_celldata.PhyIDVector_Global.resize(k,0);
        t_celldata.PhyNameVector_Global.resize(k);
        maxdim=-1;
        for(k=0;k<t_celldata.PhyGroupNum_Global-1;k++){
            MPIDataBus::receiveIntegerFromMaster( phyid,rank*10000+k*1000+1);
            MPIDataBus::receiveIntegerFromMaster(   dim,rank*10000+k*1000+2);
            MPIDataBus::receiveStringFromMaster(phyname,rank*10000+k*1000+3);

            t_celldata.PhyIDVector_Global[k]=phyid;
            t_celldata.PhyDimVector_Global[k]=dim;
            t_celldata.PhyNameVector_Global[k]=phyname;

            if(dim>maxdim) maxdim=dim;

            MPIDataBus::receivePhyID2MeshCellMapFromMaster(t_celldata.PhyID2MeshCellVectorMap_Local,10000*rank+k*1000+4);
            MPIDataBus::receivePhyName2MeshCellMapFromMaster(t_celldata.PhyName2MeshCellVectorMap_Local,10000*rank+k*1000+4+20);
        }
        k=t_celldata.PhyGroupNum_Global-1;
        t_celldata.PhyIDVector_Global[k]=0;
        t_celldata.PhyDimVector_Global[k]=maxdim;
        t_celldata.PhyNameVector_Global[k]="alldomain";

        MPIDataBus::receiveMeshCellFromMaster(t_celldata.MeshCell_Local,10000*rank+k*1000+4+40);
        MPIDataBus::receivePhyID2MeshCellMapFromMaster(t_celldata.PhyID2MeshCellVectorMap_Local,10000*rank+k*1000+4+60);
        MPIDataBus::receivePhyName2MeshCellMapFromMaster(t_celldata.PhyName2MeshCellVectorMap_Local,10000*rank+k*1000+4+80);
    }

    return true;
}