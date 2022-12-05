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
//+++ Date   : 2022.08.26
//+++ Purpose: Implement the msh file (version-2) import function.
//+++          This mesh file must be the *.msh in version-2.
//+++          For version-4, please use Msh4FileImporter.
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Mesh/Msh2FileImporter.h"

int Msh2FileImporter::getMaxMeshDim(const string &filename)const{
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

bool Msh2FileImporter::importMeshFile(const string &filename,MeshData &meshdata){
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

    mshMaxDim=getMaxMeshDim(filename);

    meshdata.m_maxdim=mshMaxDim;

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
            if(version<2.0 || version>2.2){
                MessagePrinter::printErrorTxt("version="+to_string(version)+" is not supported for msh2 file importer, "
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
        else if(str.find("$Nodes")!=string::npos){
            // read the nodes' coordinates
            // node-id, x, y, z
            meshdata.m_nodes=0;
            in>>meshdata.m_nodes;
            int nodeid;
            double x,y,z;
            meshdata.m_xmin=meshdata.m_ymin=meshdata.m_zmin= 1.0e16;
            meshdata.m_xmax=meshdata.m_ymax=meshdata.m_zmax=-1.0e16;
            meshdata.m_nodecoords0.resize(meshdata.m_nodes*3,0.0);
            meshdata.m_nodecoords.resize(meshdata.m_nodes*3,0.0);
            for(int i=0;i<meshdata.m_nodes;i++){
                in>>nodeid>>x>>y>>z;
                meshdata.m_nodecoords0[(nodeid-1)*3+1-1]=x;
                meshdata.m_nodecoords0[(nodeid-1)*3+2-1]=y;
                meshdata.m_nodecoords0[(nodeid-1)*3+3-1]=z;

                meshdata.m_nodecoords[(nodeid-1)*3+1-1]=x;
                meshdata.m_nodecoords[(nodeid-1)*3+2-1]=y;
                meshdata.m_nodecoords[(nodeid-1)*3+3-1]=z;

                if(x>meshdata.m_xmax) meshdata.m_xmax=x;
                if(x<meshdata.m_xmin) meshdata.m_xmin=x;
                if(y>meshdata.m_ymax) meshdata.m_ymax=y;
                if(y<meshdata.m_ymin) meshdata.m_ymin=y;
                if(z>meshdata.m_zmax) meshdata.m_zmax=z;
                if(z<meshdata.m_zmin) meshdata.m_zmin=z;
            }
            getline(in,str);
        }// end-of-node-coordinates-reading
        else if(str.find("$Elements")!=string::npos){
            meshdata.m_elements=0;
            in>>meshdata.m_elements;// total elements number
            vector<int> tempconn;
            int elmtid,phyid,geoid,ntags,elmttype,vtktype;
            int nodes,dim,elmtorder;
            string meshtypename;
            MeshType meshtype;

            meshdata.m_pointelmt_connectivity.clear();
            meshdata.m_lineelmt_connectivity.clear();
            meshdata.m_surfaceelmt_connectivity.clear();
            meshdata.m_bulkelmt_connectivity.clear();

            meshdata.m_pointelmts=0;
            meshdata.m_lineelmts=0;
            meshdata.m_surfaceelmts=0;
            meshdata.m_bulkelmts=0;

            meshdata.m_lineelmt_type=MeshType::EDGE2;
            meshdata.m_surfaceelmt_type=MeshType::TRI3;

            ElmtPhyIDVec.resize(meshdata.m_elements,0);
            ElmtDimVec.resize(meshdata.m_elements,0);
            ElmtConn.resize(meshdata.m_elements);

            for(int e=0;e<meshdata.m_elements;e++){
                in>>elmtid>>elmttype>>ntags>>phyid>>geoid;

                nodes=MshFileUtils::getElmtNodesNumFromElmtType(elmttype);
                dim=MshFileUtils::getElmtDimFromElmtType(elmttype);
                vtktype=MshFileUtils::getElmtVTKCellTypeFromElmtType(elmttype);
                meshtype=MshFileUtils::getElmtMeshTypeFromElmtType(elmttype);
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
                    meshdata.m_pointelmts+=1;
                    meshdata.m_pointelmt_volume.push_back(0.0);
                    meshdata.m_pointelmt_connectivity.push_back(tempconn);
                }
                
                if(dim==1 && dim<mshMaxDim){
                    meshdata.m_lineelmts+=1;
                    meshdata.m_lineelmt_connectivity.push_back(tempconn);
                    meshdata.m_lineelmt_type=meshtype;
                    meshdata.m_lineelmt_volume.push_back(0.0);
                    meshdata.m_nodesperlineelmt=nodes;
                }
                if(dim==2 && dim<mshMaxDim){
                    meshdata.m_surfaceelmts+=1;
                    meshdata.m_surfaceelmt_connectivity.push_back(tempconn);
                    meshdata.m_surfaceelmt_type=meshtype;
                    meshdata.m_surfaceelmt_volume.push_back(0.0);
                    meshdata.m_nodespersurfaceelmt=nodes;
                }
                if(dim==mshMaxDim){
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
            }// end-of-element-loop-for-element-reading

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
                    meshdata.m_phygroup_name2bulkelmtidvec[meshdata.m_phygroups-1].second.push_back(e+1);
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
                if(phyid>maxPhyGroupID){
                    meshdata.m_phygroup_dimvec.push_back(mshMaxDim);
                    meshdata.m_phygroup_phyidvec.push_back(phyid);
                    meshdata.m_phygroup_phynamevec.push_back(phyname);

                    meshdata.m_phygroup_name2dimvec.push_back(make_pair(phyname,mshMaxDim));
                    meshdata.m_phygroup_name2phyidvec.push_back(make_pair(phyname,phyid));
                    meshdata.m_phygroup_phyid2namevec.push_back(make_pair(phyid,phyname));
                    meshdata.m_phygroup_nodesnumperelmtvec.push_back(meshdata.m_nodesperbulkelmt);
                }
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
                        dim=meshdata.m_phygroup_dimvec[i];
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
        for(int e=0;e<meshdata.m_elements;e++){
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