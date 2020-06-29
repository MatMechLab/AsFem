//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.06.26
//+++ Purpose: implement the mesh import for msh(ver=2) file format
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "Mesh/Gmsh2IO.h"

bool Gmsh2IO::ReadMeshFromFile(Mesh &mesh){
    char longbuff[105],buff[55];
    string str,substr;
    if(!_HasSetMeshFileName){
        MessagePrinter::PrintErrorTxt("can\'t read mesh, the mesh file name has not been set");
        return false;
    }
    _in.open(_MeshFileName.c_str(),ios::in);
    if(!_in.is_open()){
        snprintf(longbuff,105,"can\'t read the .msh file(=%25s)      ,please make sure file name is correct",_MeshFileName.c_str());
        str=longbuff;
        MessagePrinter::PrintErrorTxt(str);
        MessagePrinter::AsFem_Exit();
    }
    double version;
    int format,size;

    //*** initialize
    mesh.GetNodeCoordsPtr().clear();
    mesh.GetElmtConnPtr().clear();
    mesh.GetElmtVolumePtr().clear();

    mesh.GetElmtVTKCellTypeListPtr().clear();
    mesh.GetElmtPhyIDListPtr().clear();
    mesh.GetElmtDimListPtr().clear();
    mesh.GetElmtMeshTypeListPtr().clear();
    //*** for physical groups
    mesh.GetPhysicalGroupNameListPtr().clear();
    mesh.GetPhysicalGroupIDListPtr().clear();
    mesh.GetPhysicalGroupDimListPtr().clear();
    mesh.GetPhysicalGroupID2NameListPtr().clear();
    mesh.GetPhysicalGroupName2IDListPtr().clear();
    mesh.GetPhysicalGroupName2NodesNumPerElmtListPtr().clear();
    mesh.GetPhysicalName2ElmtIDsListPtr().clear();
    mesh.GetPhysicalName2NodeIDsListPtr().clear();

    mesh.SetDim(0);
    mesh.SetNodesNum(0);
    mesh.SetElmtsNum(0);
    mesh.SetBulkElmtsNum(0);
    mesh.SetSurfaceElmtsNum(0);
    mesh.SetLineElmtsNum(0);
    mesh.SetPhysicalGroupNums(0);

    _nMaxDim=-1;_nMinDim=4;
    _nPhysicGroups=0;

    vector<pair<int,int>> UniquePhyDim2IDList;
    int MaxPhyIDofPhyGroup=-1,MaxPhyIDofElmt=-1;

    while(!_in.eof()){
        // now we start to read *.msh file
        getline(_in,str);
        if(str.find("$MeshFormat")!=string::npos){
            _in>>version>>format>>size;
            if((version!=2.0)&&(version!=2.1)&&(version!=2.2)){
                // currently, only the gmsh2 format is supported!!!
                snprintf(buff,50,"version=%12.3f is not supported yet",version);
                str=buff;
                MessagePrinter::PrintErrorTxt(str);
                return false;
            }
        }// end-of-versioin-read
        else if(str.find("$PhysicalNames")!=string::npos){
            _nPhysicGroups=0;
            mesh.GetPhysicalGroupNameListPtr().clear();
            mesh.GetPhysicalGroupIDListPtr().clear();
            mesh.GetPhysicalGroupDimListPtr().clear();
            mesh.GetPhysicalGroupID2NameListPtr().clear();
            mesh.GetPhysicalGroupName2IDListPtr().clear();
            mesh.GetPhysicalName2ElmtIDsListPtr().clear();
            mesh.GetPhysicalName2NodeIDsListPtr().clear();
            mesh.GetPhysicalGroupName2NodesNumPerElmtListPtr().clear();
            mesh.SetPhysicalGroupNums(0);

            int phydim,phyid;
            string phyname;
            _in>>_nPhysicGroups;
            getline(_in,str);//remove \n in this line
            for(int i=0;i<_nPhysicGroups;i++){
                getline(_in,str);
                istringstream s_stream(str);
                s_stream>>phydim>>phyid>>phyname;
                // remove the '"' ,keep only the text
                phyname.erase(std::remove(phyname.begin(),phyname.end(),'"'),phyname.end());
                mesh.GetPhysicalGroupNameListPtr().push_back(phyname);
                mesh.GetPhysicalGroupIDListPtr().push_back(phyid);
                mesh.GetPhysicalGroupDimListPtr().push_back(phydim);

                mesh.GetPhysicalGroupID2NameListPtr().push_back(make_pair(phyid,phyname));
                mesh.GetPhysicalGroupName2IDListPtr().push_back(make_pair(phyname,phyid));

                if(phydim>_nMaxDim) _nMaxDim=phydim;
                if(_nMinDim<phydim) _nMinDim=phydim;
                if(phyid>MaxPhyIDofPhyGroup) MaxPhyIDofPhyGroup=phyid;
            }
        }//end-of-physical-group-information
        else if(str.find("$Nodes")!=string::npos){
            // read the nodes' coordinates
            // node-id, x, y, z
            _nNodes=0;
            _in>>_nNodes;
            mesh.GetNodeCoordsPtr().resize(_nNodes*3,0.0);
            int id;
            double x,y,z;
            _Xmax=-1.0e16;_Xmin=1.0e16;
            _Ymax=_Xmax;_Ymin=_Xmin;
            _Zmax=_Xmax;_Zmin=_Xmin;
            for(int i=0;i<_nNodes;i++){
                _in>>id>>x>>y>>z;
                mesh.GetNodeCoordsPtr()[(id-1)*3+1-1]=x;
                mesh.GetNodeCoordsPtr()[(id-1)*3+2-1]=y;
                mesh.GetNodeCoordsPtr()[(id-1)*3+3-1]=z;

                if(x>_Xmax) _Xmax=x;
                if(x<_Xmin) _Xmin=x;
                if(y>_Ymax) _Ymax=y;
                if(y<_Ymin) _Ymin=y;
                if(z>_Zmax) _Zmax=z;
                if(z<_Zmin) _Zmin=z;
            }
            getline(_in,str);
        }//end-of-node-reading
        else if(str.find("$Elements")!=string::npos){
            // here the element constains all the element(no matter it is line elmt or bulk elmt)
            _nElmts=0;
            _in>>_nElmts;
            mesh.SetElmtsNum(_nElmts);
            mesh.GetElmtConnPtr().resize(_nElmts,vector<int>(0));
            mesh.GetElmtVTKCellTypeListPtr().resize(_nElmts,0);
            mesh.GetElmtVTKCellTypeListPtr().resize(_nElmts,0.0);
            mesh.GetElmtPhyIDListPtr().resize(_nElmts,0);
            mesh.GetElmtMeshTypeListPtr().resize(_nElmts,MeshType::NULLTYPE);

            int elmtid,phyid,geoid,ntags,elmttype,vtktype;
            int nodes,dim,maxdim,elmtorder;
            vector<int> tempconn;
            MeshType meshtype,bcmeshtype;

            _nNodesPerBulkElmt=-1;
            _nNodesPerLineElmt=0;
            _nNodesPerSurfaceElmt=0;
            maxdim=-1;
            _nOrder=1;

            bool IsUnique;

        

            for(int e=0;e<_nElmts;e++){
                _in>>elmtid>>elmttype>>ntags>>phyid>>geoid;

                if(phyid==0) phyid=geoid;

                nodes=GetElmtNodesNumFromGmshElmtType(elmttype);
                dim=GetElmtDimFromGmshElmtType(elmttype);
                vtktype=GetElmtVTKCellTypeFromGmshElmtType(elmttype);
                meshtype=GetElmtMeshTypeFromGmshElmtType(elmttype);
                bcmeshtype=GetSubElmtMeshTypeFromGmshElmtType(elmttype);
                elmtorder=GetElmtOrderFromGmshElmtType(elmttype);
                if(elmtorder>_nOrder) _nOrder=elmtorder;

                if(dim==1){
                    _nNodesPerLineElmt=nodes;
                    _nNodesPerBulkElmt=nodes;
                    mesh.SetMeshType(meshtype);
                    mesh.SetLineMeshType(meshtype);
                }
                if(dim==2){
                    _nNodesPerSurfaceElmt=nodes;
                    _nNodesPerBulkElmt=nodes;
                    mesh.SetMeshType(meshtype);
                    mesh.SetLineMeshType(bcmeshtype);
                }
                if(dim==3){
                    _nNodesPerBulkElmt=nodes;
                    mesh.SetMeshType(meshtype);
                    mesh.SetSurfaceMeshType(bcmeshtype);
                }
                
                if(dim>maxdim) maxdim=dim;
                if(dim>_nMaxDim) _nMaxDim=dim;
                if(dim<_nMinDim) _nMinDim=dim;
                if(phyid>MaxPhyIDofElmt) MaxPhyIDofElmt=phyid;

                mesh.GetElmtConnPtr()[elmtid-1].resize(nodes+1,0);
                mesh.GetElmtConnPtr()[elmtid-1][0]=nodes;

                tempconn.resize(nodes,0);
                for(int j=0;j<nodes;j++){
                    _in>>tempconn[j];
                }

                for(int j=0;j<nodes;j++){
                    mesh.GetElmtConnPtr()[elmtid-1][j+1]=tempconn[j];
                }

                mesh.GetElmtDimListPtr()[elmtid-1]=dim;
                mesh.GetElmtMeshTypeListPtr()[elmtid-1]=meshtype;
                mesh.GetElmtPhyIDListPtr()[elmtid-1]=phyid;
                mesh.GetElmtVTKCellTypeListPtr()[elmtid-1]=vtktype;

                if(UniquePhyDim2IDList.size()==0){
                    UniquePhyDim2IDList.push_back(make_pair(dim,phyid));
                }
                else{
                    IsUnique=true;
                    for(const auto &it:UniquePhyDim2IDList){
                        if(it.first==dim && it.second==phyid){
                            IsUnique=false;
                            break;
                        }
                    }
                    if(IsUnique){
                        UniquePhyDim2IDList.push_back(make_pair(dim,phyid));
                    }
                }

            }
            mesh.SetMeshOrder(_nOrder);
            mesh.SetMinDim(_nMinDim);
            mesh.SetMaxDim(_nMaxDim);
            mesh.SetNodesNumPerBulkElmt(_nNodesPerBulkElmt);
            mesh.SetNodesNumPerSurfaceElmt(_nNodesPerSurfaceElmt);
            mesh.SetNodesNumPerLineElmt(_nNodesPerLineElmt);
            
        }//end-of-read-element-info
        
    }// end-of-while(!_in.eof())


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
                mesh.GetPhysicalGroupIDListPtr().push_back(phyid);
                mesh.GetPhysicalGroupDimListPtr().push_back(_nMaxDim);
                mesh.GetPhysicalGroupNameListPtr().push_back(phyname);

                mesh.GetPhysicalGroupID2NameListPtr().push_back(make_pair(phyid,phyname));
                mesh.GetPhysicalGroupName2IDListPtr().push_back(make_pair(phyname,phyid));
            }
        }
        _nBulkElmts=0;
        _nSurfaceElmts=0;
        _nLineElmts=0;
        

        bulkconn.clear();
        for(int e=0;e<_nElmts;e++){
            if(mesh.GetIthElmtDim(e+1)==_nMaxDim){
                _nBulkElmts+=1;
                bulkconn.push_back(e+1);
                for(int j=0;j<static_cast<int>(phy2conn.size());j++){
                    if(mesh.GetIthElmtPhyID(e+1)==UniquePhyDim2IDList[j].second){
                        phy2conn[j].push_back(e+1);
                        break;
                    }
                }
            }
            else{
                // for lower dimension elements
                if(mesh.GetIthElmtDim(e+1)==2){
                    _nSurfaceElmts+=1;
                }
                else if(mesh.GetIthElmtPhyID(e+1)==1){
                    _nLineElmts+=1;
                }
            }
        }

        mesh.SetBulkElmtsNum(_nBulkElmts);
        mesh.SetSurfaceElmtsNum(_nSurfaceElmts);
        mesh.SetLineElmtsNum(_nLineElmts);

        for(int j=0;j<static_cast<int>(phy2conn.size());j++){
            phyname=mesh.GetIthPhysicalName(j+1);
            mesh.GetPhysicalName2ElmtIDsListPtr().push_back(make_pair(phyname,phy2conn[j]));
        }

        mesh.SetPhysicalGroupNums(1+phy2conn.size());
        mesh.GetPhysicalGroupNameListPtr().push_back("alldomain");
        mesh.GetPhysicalGroupDimListPtr().push_back(_nMaxDim);
        mesh.GetPhysicalGroupIDListPtr().push_back(maxid+1);
        mesh.GetPhysicalGroupID2NameListPtr().push_back(make_pair(maxid+1,"alldomain"));
        mesh.GetPhysicalGroupName2IDListPtr().push_back(make_pair("alldomain",maxid+1));
        mesh.GetPhysicalName2ElmtIDsListPtr().push_back(make_pair("alldomain",bulkconn));

        return _nBulkElmts>0;

    }
    else if(_nPhysicGroups==1&&UniquePhyDim2IDList.size()==1){
        // for the case where only 1 physical group is defined
        _nBulkElmts=0;
        _nSurfaceElmts=0;
        _nLineElmts=0;
        bulkconn.clear();
        int maxid=-1;
        for(int e=0;e<_nElmts;e++){
            _nBulkElmts+=1;
            bulkconn.push_back(e+1);
            if(mesh.GetIthElmtPhyID(e+1)>maxid) maxid=mesh.GetIthElmtPhyID(e+1);
        }

        mesh.SetBulkElmtsNum(_nBulkElmts);
        mesh.SetSurfaceElmtsNum(_nSurfaceElmts);
        mesh.SetLineElmtsNum(_nLineElmts);

        mesh.SetPhysicalGroupNums(2);
        mesh.GetPhysicalGroupNameListPtr().push_back("alldomain");
        mesh.GetPhysicalGroupIDListPtr().push_back(maxid+1);
        mesh.GetPhysicalGroupName2IDListPtr().push_back(make_pair("alldomain",maxid+1));
        mesh.GetPhysicalGroupID2NameListPtr().push_back(make_pair(maxid+1,"alldomain"));

        mesh.GetPhysicalName2ElmtIDsListPtr().push_back(make_pair(mesh.GetIthPhysicalName(1),bulkconn));
        mesh.GetPhysicalName2ElmtIDsListPtr().push_back(make_pair("alldomain",bulkconn));

        return _nBulkElmts>0;

    }
    else{
        int maxid=-1;
        string phyname;
        bulkconn.clear();
        phy2conn.resize(_nPhysicGroups,vector<int>(0));
        for(int e=0;e<_nElmts;e++){
            if(mesh.GetIthElmtDim(e+1)==_nMaxDim){
                _nBulkElmts+=1;
                bulkconn.push_back(e+1);
                for(int j=0;j<static_cast<int>(phy2conn.size());j++){
                    if(mesh.GetIthElmtPhyID(e+1)==mesh.GetIthPhysicalID(j+1)){
                        phy2conn[j].push_back(e+1);
                        break;
                    }
                }
            }
            else{
                // for lower dimension elements
                if(mesh.GetIthElmtDim(e+1)==2){
                    _nSurfaceElmts+=1;
                }
                else if(mesh.GetIthElmtPhyID(e+1)==1){
                    _nLineElmts+=1;
                }
            }
        }

        mesh.SetBulkElmtsNum(_nBulkElmts);
        mesh.SetSurfaceElmtsNum(_nSurfaceElmts);
        mesh.SetLineElmtsNum(_nLineElmts);

        for(int j=0;j<static_cast<int>(phy2conn.size());j++){
            phyname=mesh.GetIthPhysicalName(j+1);
            mesh.GetPhysicalName2ElmtIDsListPtr().push_back(make_pair(phyname,phy2conn[j]));
        }

        mesh.SetPhysicalGroupNums(1+phy2conn.size());
        mesh.GetPhysicalGroupNameListPtr().push_back("alldomain");
        mesh.GetPhysicalGroupDimListPtr().push_back(_nMaxDim);
        mesh.GetPhysicalGroupIDListPtr().push_back(maxid+1);
        mesh.GetPhysicalGroupID2NameListPtr().push_back(make_pair(maxid+1,"alldomain"));
        mesh.GetPhysicalGroupName2IDListPtr().push_back(make_pair("alldomain",maxid+1));
        mesh.GetPhysicalName2ElmtIDsListPtr().push_back(make_pair("alldomain",bulkconn));

        return _nBulkElmts>0;
    }

    return false;
}