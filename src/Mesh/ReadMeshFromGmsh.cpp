//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "Mesh/Mesh.h"

bool Mesh::ReadMeshFromGmsh(){
    ifstream in;
    in.open(_GmshFileName,ios::in);
    if(!in.is_open()){
        PetscPrintf(PETSC_COMM_WORLD,"*** Error: can\'t read the .msh file(=%20s)      !!!   ***\n",_GmshFileName.c_str());
        PetscPrintf(PETSC_COMM_WORLD,"***        please make sure file name is correct                !!!   ***\n");
        Msg_AsFem_Exit();
    }


    double version;
    int format,size;
    string str,substr;
    _nMaxDim=-1;_nMinDim=4;
    _nPhysicGroups=0;
    _PhysicGroupNameList.clear();
    _PhysicIDToNameList.clear();
    _PhysicNameToIDList.clear();
    _PhysicNameToElmtIndexSet.clear();

    _nPhysicGroups=0;
    _PhysicGroupNameList.clear();
    _PhysicIDToNameList.clear();
    _PhysicNameToIDList.clear();

    _PhyGroupDimVec.clear();
    _PhyGroupArray.clear();


    while(!in.eof()){
        // now we start to read *.msh file
        getline(in,str);

        if(str.find("$MeshFormat")!=string::npos){
            in>>version>>format>>size;
            if((version!=2.0)&&(version!=2.1)&&(version!=2.2)){
                // currently, only the gmsh2 format is supported!!!
                PetscPrintf(PETSC_COMM_WORLD,"*** Error: version=%12.3f is not supported yet        !!!   ***\n");
                return false;
            }
        }
        else if(str.find("$PhysicalNames")!=string::npos){
            // start to read the physical group information
            // remember!!! it is also possible that
            // msh file dosent has any physical information
            // read PhysicalName group information
            _nPhysicGroups=0;
            _PhysicGroupNameList.clear();
            _PhysicIDToNameList.clear();
            _PhysicNameToIDList.clear();

            _PhyGroupDimVec.clear();
            _PhyGroupArray.clear();

            int phydim,phyid;
            string phyname;
            in>>_nPhysicGroups;
            getline(in,str);//remove \n in this line
            for(int i=0;i<_nPhysicGroups;i++){
                getline(in,str);
                istringstream s_stream(str);
                s_stream>>phydim>>phyid>>phyname;
                // remove the '"' ,keep only the text
                phyname.erase(std::remove(phyname.begin(),phyname.end(),'"'),phyname.end());
                _PhyGroupDimVec.push_back(phydim);
                _PhyGroupArray.push_back(make_pair(phyid,phyname));
                if(phydim>_nMaxDim) _nMaxDim=phydim;
                if(_nMinDim<phydim) _nMinDim=phydim;
            }
        }
        else if(str.find("$Nodes")!=string::npos){
            // read the nodes' coordinates
            // node-id, x, y, z
            _nNodes=0;
            in>>_nNodes;
            _NodeCoords.resize(_nNodes*4,0.0);
            int id;
            double x,y,z;

            _Xmax=-1.0e16;_Xmin=1.0e16;
            _Ymax=_Xmax;_Ymin=_Xmin;
            _Zmax=_Xmax;_Zmin=_Xmin;
            for(long int i=0;i<_nNodes;i++){
                in>>id>>x>>y>>z;
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
            getline(in,str);
        }
        else if(str.find("$Elements")!=string::npos){
            // here the element constains all the element(no matter it is line elmt or bulk elmt)
            in>>_nElmts;// number of total elements
            _ElmtConn.resize(_nElmts,vector<int>(0));
            _ElmtVTKCellType.resize(_nElmts,0);
            _ElmtVolume.resize(_nElmts,0.0);
            _ElmtDimVec.resize(_nElmts,0);
            _ElmtTypeVec.resize(_nElmts,0);
            _ElmtPhyIDVec.resize(_nElmts,0);
            _ElmtGeoIDVec.resize(_nElmts,0);

            _MeshUniGeoID.clear();
            _MeshUniPhyID.clear();
            _MeshUniPhyDim.clear();

            int elmtid,phyid,geoid,ntags,elmttype,vtktype;
            int nodes,dim,maxdim,elmtorder;
            vector<int> tempconn;
            MeshType meshtype,bcmeshtype;

            _nNodesPerBulkElmt=-1;
            _nNodesPerLineElmt=0;
            _nNodesPerSurfaceElmt=0;
            maxdim=-1;
            _nOrder=1;
            for(int e=0;e<_nElmts;++e){
                in>>elmtid>>elmttype>>ntags>>phyid>>geoid;

                nodes=GetNodesNumViaGmshElmtType(elmttype);
                dim=GetElmtDimViaGmshElmtType(elmttype);
                vtktype=GetElmtVTKCellTypeViaGmshElmtType(elmttype);
                meshtype=GetElmtTypeViaGmshElmtType(elmttype);
                bcmeshtype=GetBCElmtTypeViaGmshElmtType(elmttype);
                elmtorder=GetElmtOrderViaGmshElmtType(elmttype);
                if(elmtorder>_nOrder) _nOrder=elmtorder;


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


                _ElmtConn[elmtid-1].resize(nodes+1,0);
                _ElmtConn[elmtid-1][0]=nodes;
                tempconn.resize(nodes,0);
                for(int j=0;j<nodes;j++){
                    in>>tempconn[j];
                }
                // modify the order
                ModifyElmtConnViaGmshElmtType(elmttype,tempconn);
                for(int j=0;j<nodes;j++){
                    _ElmtConn[elmtid-1][j+1]=tempconn[j];
                }
                _ElmtDimVec[elmtid-1]=dim;
                _ElmtTypeVec[elmtid-1]=elmttype;
                _ElmtPhyIDVec[elmtid-1]=phyid;
                _ElmtGeoIDVec[elmtid-1]=geoid;
                _ElmtVTKCellType[elmtid-1]=vtktype;

                // set the unique id for geometry and physical
                if(_MeshUniPhyID.size()==0){
                    _MeshUniGeoID.push_back(geoid);
                    _MeshUniPhyID.push_back(phyid);
                    _MeshUniPhyDim.push_back(dim);
                }
                else{
                    bool IsUni=true;
                    for(unsigned int i=0;i<_MeshUniPhyID.size();i++){
                        if(phyid==_MeshUniPhyID[i]){
                            IsUni=false;
                            break;
                        }
                    }
                    if(IsUni){
                        _MeshUniGeoID.push_back(geoid);
                        _MeshUniPhyID.push_back(phyid);
                        _MeshUniPhyDim.push_back(dim);
                    }
                }
            }
        }
    }

    // now we check wether all the elements' phy id is already list in PhyGroup info
    // if not, we will add the new phy info to the list
    if(_nPhysicGroups<1){
        //if phy group is empty; use UniPhyID as the phy group information
        //thus the physical name will be given as "number", the number to be the physical name
        // _BulkPhyGroupNameList.clear();
        // _BounPhyGroupNameList.clear();
        int maxid=-1;
        for(int i=0;i<int(_MeshUniPhyID.size());++i){
            // only the max dim related id will be put into PhyGroup
            // for lower dimension? No, since you don't define the physical information for them
            // they will be ignored, if and only if you remember to define them in Gmsh
            if(_MeshUniPhyDim[i]==_nMaxDim){
                _PhyGroupDimVec.push_back(_MeshUniPhyDim[i]);
                _PhyGroupArray.push_back(make_pair(_MeshUniPhyID[i],to_string(_MeshUniPhyID[i])));
                if(_MeshUniPhyID[i]>maxid) maxid=_MeshUniPhyID[i];
                
            }
        }
        _nBulkElmts=0;
        vector<int> tempconn;
        vector<vector<int>> phyconn(_PhyGroupArray.size(),vector<int>(0));
        tempconn.clear();
        for(int e=0;e<_nElmts;++e){
            if(_ElmtDimVec[e]==_nMaxDim){
                _nBulkElmts+=1;
                tempconn.push_back(e+1);// store all the elmts for "alldomain"
                for(int j=0;j<int(phyconn.size());++j){
                    if(_ElmtPhyIDVec[e]==_MeshUniPhyID[j]){
                        phyconn[j].push_back(e+1);
                        break;
                    }
                }
            }
            else{
                continue;
            }
        }
        //************************
        _PhysicGroupNameList.clear();
        _PhysicIDToNameList.clear();
        _PhysicNameToIDList.clear();
        _PhysicNameToElmtIndexSet.clear();

        _PhysicNameToDimList.clear();
        _PhysicNameToElmtNodesNumList.clear();
        for(int j=0;j<int(phyconn.size());++j){
            _PhysicGroupNameList.push_back(_PhyGroupArray[j].second);
            _PhysicIDToNameList.push_back(make_pair(_PhyGroupArray[j].first,_PhyGroupArray[j].second));
            _PhysicNameToIDList.push_back(make_pair(_PhyGroupArray[j].second,_PhyGroupArray[j].first));
            _PhysicNameToElmtIndexSet.push_back(make_pair(_PhyGroupArray[j].second,phyconn[j]));

            _PhysicNameToDimList.push_back(make_pair(_PhyGroupArray[j].second,_PhyGroupDimVec[j]));
            _PhysicNameToElmtNodesNumList.push_back(make_pair(_PhyGroupArray[j].second,_nNodesPerBulkElmt));
        }
        _PhysicGroupNameList.push_back("alldomain");
        _PhysicIDToNameList.push_back(make_pair(maxid+1,"alldomain"));
        _PhysicNameToIDList.push_back(make_pair("alldomain",maxid+1));
        _PhysicNameToElmtIndexSet.push_back(make_pair("alldomain",tempconn));

        //****************************************
        _PhysicNameToDimList.push_back(make_pair("alldomain",_nMaxDim));
        _PhysicNameToElmtNodesNumList.push_back(make_pair("alldomain",_nNodesPerBulkElmt));

        return (_nBulkElmts>0);
    }
    else if(_nPhysicGroups==1&&_MeshUniPhyID.size()==1){
        // only the bulk elmts is defined (in this case, on BC elmts is defined!!!)
        // then all the elmts should be stored into "alldomain"
        string blockname=_PhyGroupArray[0].second;
        _PhysicGroupNameList.push_back(blockname);
        _PhysicIDToNameList.push_back(make_pair(_PhyGroupArray[0].first,blockname));
        _PhysicNameToIDList.push_back(make_pair(blockname,_PhyGroupArray[0].first));

        // _BulkPhyGroupNameList.clear();
        // _BounPhyGroupNameList.clear();

        // _BulkPhyGroupNameList.push_back(blockname);

        _nBulkElmts=0;
        vector<int> tempconn;
        tempconn.clear();
        for(int e=0;e<_nElmts;++e){
            _nBulkElmts+=1;
            tempconn.push_back(e+1);// store all the elmts for "alldomain"
        }
        _PhysicGroupNameList.push_back("alldomain");
        _PhysicIDToNameList.push_back(make_pair(_PhyGroupArray[0].first,"alldomain"));
        _PhysicNameToIDList.push_back(make_pair("alldomain",_PhyGroupArray[0].first));

        _PhysicNameToElmtIndexSet.push_back(make_pair(blockname,tempconn));
        _PhysicNameToElmtIndexSet.push_back(make_pair("alldomain",tempconn));

        // _BulkPhyGroupNameList.push_back("alldomain");
        _PhysicNameToDimList.push_back(make_pair("alldomain",_nMaxDim));
        _PhysicNameToElmtNodesNumList.push_back(make_pair("alldomain",_nNodesPerBulkElmt));

        return _nBulkElmts>0;
    }
    else{
        // for the case, which defined both the bulk phy and bc phy
        int maxid=-1;
        for(int i=0;i<int(_MeshUniPhyID.size());++i){
            bool IsInList=false;
            if(_MeshUniPhyID[i]>maxid) maxid=_MeshUniPhyID[i];
            for(int j=0;j<int(_PhyGroupArray.size());++j){
                if(_PhyGroupArray[j].first>maxid) maxid=_PhyGroupArray[j].first;

                if(_MeshUniPhyID[i]==_PhyGroupArray[j].first){
                    // if the phy id is not in the list of PhyInfo array
                    // add it as the new phy
                    IsInList=true;
                    break;
                }
            }
            if(!IsInList){
                _PhyGroupDimVec.push_back(_MeshUniPhyDim[i]);
                _PhyGroupArray.push_back(make_pair(_MeshUniPhyID[i],"block"+to_string(_MeshUniPhyID[i])));
            }
        }
        _nBulkElmts=0;
        vector<int> tempconn;
        vector<vector<int>> phyconn(_PhyGroupArray.size(),vector<int>(0));
        tempconn.clear();
        for(int e=0;e<_nElmts;++e){
            if(_ElmtDimVec[e]==_nMaxDim){
                _nBulkElmts+=1;
                tempconn.push_back(e+1);// store all the elmts for "alldomain"
            }
            
            for(int j=0;j<int(phyconn.size());++j){
                if(_ElmtPhyIDVec[e]==_PhyGroupArray[j].first){
                    phyconn[j].push_back(e+1);
                    break;
                }
            }
        }
        //************************
        _PhysicGroupNameList.clear();
        _PhysicIDToNameList.clear();
        _PhysicNameToIDList.clear();
        _PhysicNameToElmtIndexSet.clear();
        // _BounPhyGroupNameList.clear();
        // _BulkPhyGroupNameList.clear();
        _PhysicNameToDimList.clear();
        _PhysicNameToElmtNodesNumList.clear();

        for(int j=0;j<int(phyconn.size());++j){
            _PhysicGroupNameList.push_back(_PhyGroupArray[j].second);
            _PhysicIDToNameList.push_back(make_pair(_PhyGroupArray[j].first,_PhyGroupArray[j].second));
            _PhysicNameToIDList.push_back(make_pair(_PhyGroupArray[j].second,_PhyGroupArray[j].first));
            _PhysicNameToElmtIndexSet.push_back(make_pair(_PhyGroupArray[j].second,phyconn[j]));
            // if(_PhyGroupDimVec[j]==_nMaxDim){
            //     _BulkPhyGroupNameList.push_back(_PhyGroupArray[j].second);
            // }
            
            _PhysicNameToDimList.push_back(make_pair(_PhyGroupArray[j].second,_PhyGroupDimVec[j]));
            if(_nMaxDim==1){
                if(_PhyGroupDimVec[j]==0){
                    _PhysicNameToElmtNodesNumList.push_back(make_pair(_PhyGroupArray[j].second,1));
                }
                else{
                    _PhysicNameToElmtNodesNumList.push_back(make_pair(_PhyGroupArray[j].second,_nNodesPerBulkElmt));
                }
            }
            else if(_nMaxDim==2){
                if(_PhyGroupDimVec[j]==0){
                    _PhysicNameToElmtNodesNumList.push_back(make_pair(_PhyGroupArray[j].second,1));
                }
                else if(_PhyGroupDimVec[j]==1){
                    _PhysicNameToElmtNodesNumList.push_back(make_pair(_PhyGroupArray[j].second,_nNodesPerLineElmt));
                }
                else{
                    _PhysicNameToElmtNodesNumList.push_back(make_pair(_PhyGroupArray[j].second,_nNodesPerBulkElmt));
                }
            }
            else if(_nMaxDim==3){
                if(_PhyGroupDimVec[j]==0){
                    _PhysicNameToElmtNodesNumList.push_back(make_pair(_PhyGroupArray[j].second,1));
                }
                else if(_PhyGroupDimVec[j]==1){
                    _PhysicNameToElmtNodesNumList.push_back(make_pair(_PhyGroupArray[j].second,_nNodesPerLineElmt));
                }
                else if(_PhyGroupDimVec[j]==2){
                    _PhysicNameToElmtNodesNumList.push_back(make_pair(_PhyGroupArray[j].second,_nNodesPerSurfaceElmt));
                }
                else{
                    _PhysicNameToElmtNodesNumList.push_back(make_pair(_PhyGroupArray[j].second,_nNodesPerBulkElmt));
                }
            }
        }
        _PhysicGroupNameList.push_back("alldomain");
        _PhysicIDToNameList.push_back(make_pair(maxid+1,"alldomain"));
        _PhysicNameToIDList.push_back(make_pair("alldomain",maxid+1));
        _PhysicNameToElmtIndexSet.push_back(make_pair("alldomain",tempconn));

        // _BulkPhyGroupNameList.push_back("alldomain");

        _PhysicNameToDimList.push_back(make_pair("alldomain",_nMaxDim));
        _PhysicNameToElmtNodesNumList.push_back(make_pair("alldomain",_nNodesPerBulkElmt));

        return (_nBulkElmts>0);
    }
    return false;
}