#include "Mesh/Mesh.h"

bool Mesh::ReadMeshFromGmsh(){
    ifstream in;
    in.open(_GmshFileName,ios::in);
    if(!in.is_open()){
        cout<<"*** Error: can\'t read .msh file(="<<_GmshFileName<<")!!!"<<endl;
        cout<<"***        please make sure file name is correct !!!  ***"<<endl;
        Msg_AsFem_Exit();
    }


    double version;
    int format,size;
    string str,substr;
    _nDim=-1;_nDimMin=4;
    _nPhyGroups=0;
    _PhyGroupNameList.clear();
    _PhyIDToNameList.clear();
    _PhyNameToIDList.clear();
    _PhyNameToElmtIndexList.clear();

    while(!in.eof()){
        // now we start to read *.msh file
        getline(in,str);

        if(str.find("$MeshFormat")!=string::npos){
            in>>version>>format>>size;
            if((version!=2.0)&&(version!=2.1)&&(version!=2.2)){
                // currently, only the gmsh2 format is supported!!!
                cout<<"*** Error: version="<<version<<" is not supported !!!"<<endl;
                return false;
            }
        }
        else if(str.find("$PhysicalNames")!=string::npos){
            // start to read the physical group information
            // remember!!! it is also possible that
            // msh file dosent has any physical information
            // read PhysicalName group information
            _nPhyGroups=0;
            _PhyGroupNameList.clear();
            _PhyIDToNameList.clear();
            _PhyNameToIDList.clear();

            _PhyGroupDimVec.clear();
            _PhyGroupArray.clear();

            int phydim,phyid;
            string phyname;
            in>>_nPhyGroups;
            getline(in,str);//remove \n in this line
            for(int i=0;i<_nPhyGroups;i++){
                getline(in,str);
                istringstream s_stream(str);
                s_stream>>phydim>>phyid>>phyname;
                // remove the '"' ,keep only the text
                phyname.erase(std::remove(phyname.begin(),phyname.end(),'"'),phyname.end());
                _PhyGroupDimVec.push_back(phydim);
                _PhyGroupArray.push_back(make_pair(phyid,phyname));
                if(phydim>_nDim) _nDim=phydim;
                if(_nDimMin<phydim) _nDimMin=phydim;
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
            _CellVTKType.resize(_nElmts,0);
            _ElmtDimVec.resize(_nElmts,0);
            _ElmtTypeVec.resize(_nElmts,0);
            _ElmtPhyIDVec.resize(_nElmts,0);
            _ElmtGeoIDVec.resize(_nElmts,0);

            _MeshUniGeoID.clear();
            _MeshUniPhyID.clear();
            _MeshUniPhyDim.clear();

            int elmtid,phyid,geoid,ntags,elmttype,vtktype;
            int nodes,dim,maxdim;
            vector<int> tempconn;
            MeshType meshtype,bcmeshtype;

            _nNodesPerBulkElmt=-1;
            _nNodesPerLineElmt=0;
            _nNodesPerSurfaceElmt=0;
            maxdim=-1;
            for(int e=0;e<_nElmts;++e){
                in>>elmtid>>elmttype>>ntags>>phyid>>geoid;

                nodes=GetNodesNumViaGmshElmtType(elmttype);
                dim=GetElmtDimViaGmshElmtType(elmttype);
                vtktype=GetElmtVTKCellTypeViaGmshElmtType(elmttype);
                meshtype=GetElmtTypeViaGmshElmtType(elmttype);
                bcmeshtype=GetBCElmtTypeViaGmshElmtType(elmttype);

                if(dim==1){
                    _nNodesPerLineElmt=nodes;
                    _nNodesPerBulkElmt=nodes;
                    _LineElmtType=meshtype;
                    _BulkElmtType=meshtype;
                }
                if(dim==2){
                    _nNodesPerSurfaceElmt=nodes;
                    _nNodesPerBulkElmt=nodes;
                    _LineElmtType=bcmeshtype;
                    _SurfaceElmtType=meshtype;
                    _BulkElmtType=meshtype;
                }
                if(dim==3){
                    _nNodesPerBulkElmt=nodes;
                    _SurfaceElmtType=bcmeshtype;
                    _BulkElmtType=meshtype;
                }
                
                if(dim>maxdim) maxdim=dim;

               

                if(dim>_nDim) _nDim=dim;
                if(dim<_nDimMin) _nDimMin=dim;


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
                _CellVTKType[elmtid-1]=vtktype;

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
    if(_nPhyGroups<1){
        //if phy group is empty; use UniPhyID as the phy group information
        //thus the physical name will be given as "number", the number to be the physical name
        _BulkPhyGroupNameList.clear();
        _BounPhyGroupNameList.clear();
        int maxid=-1;
        for(int i=0;i<int(_MeshUniPhyID.size());++i){
            // only the max dim related id will be put into PhyGroup
            // for lower dimension? No, since you don't define the physical information for them
            // they will be ignored, if and only if you remember to define them in Gmsh
            if(_MeshUniPhyDim[i]==_nDim){
                _PhyGroupDimVec.push_back(_MeshUniPhyDim[i]);
                _PhyGroupArray.push_back(make_pair(_MeshUniPhyID[i],to_string(_MeshUniPhyID[i])));
                if(_MeshUniPhyID[i]>maxid) maxid=_MeshUniPhyID[i];
                _BulkPhyGroupNameList.push_back(to_string(_MeshUniPhyID[i]));
            }
        }
        _nBulkElmts=0;
        vector<int> tempconn;
        vector<vector<int>> phyconn(_PhyGroupArray.size(),vector<int>(0));
        tempconn.clear();
        for(int e=0;e<_nElmts;++e){
            if(_ElmtDimVec[e]==_nDim){
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
        _PhyGroupNameList.clear();
        _PhyIDToNameList.clear();
        _PhyNameToIDList.clear();
        _PhyNameToElmtIndexList.clear();
        for(int j=0;j<int(phyconn.size());++j){
            _PhyGroupNameList.push_back(_PhyGroupArray[j].second);
            _PhyIDToNameList.push_back(make_pair(_PhyGroupArray[j].first,_PhyGroupArray[j].second));
            _PhyNameToIDList.push_back(make_pair(_PhyGroupArray[j].second,_PhyGroupArray[j].first));
            _PhyNameToElmtIndexList.push_back(make_pair(_PhyGroupArray[j].second,phyconn[j]));
        }
        _PhyGroupNameList.push_back("alldomain");
        _PhyIDToNameList.push_back(make_pair(maxid+1,"alldomain"));
        _PhyNameToIDList.push_back(make_pair("alldomain",maxid+1));
        _PhyNameToElmtIndexList.push_back(make_pair("alldomain",tempconn));
        _BulkPhyGroupNameList.push_back("alldomain");

        return (_nBulkElmts>0);
    }
    else if(_nPhyGroups==1&&_MeshUniPhyID.size()==1){
        // only the bulk elmts is defined (in this case, on BC elmts is defined!!!)
        // then all the elmts should be stored into "alldomain"
        string blockname=_PhyGroupArray[0].second;
        _PhyGroupNameList.push_back(blockname);
        _PhyIDToNameList.push_back(make_pair(_PhyGroupArray[0].first,blockname));
        _PhyNameToIDList.push_back(make_pair(blockname,_PhyGroupArray[0].first));

        _BulkPhyGroupNameList.clear();
        _BounPhyGroupNameList.clear();

        _BulkPhyGroupNameList.push_back(blockname);

        _nBulkElmts=0;
        vector<int> tempconn;
        tempconn.clear();
        for(int e=0;e<_nElmts;++e){
            _nBulkElmts+=1;
            tempconn.push_back(e+1);// store all the elmts for "alldomain"
        }
        _PhyGroupNameList.push_back("alldomain");
        _PhyIDToNameList.push_back(make_pair(_PhyGroupArray[0].first,"alldomain"));
        _PhyNameToIDList.push_back(make_pair("alldomain",_PhyGroupArray[0].first));

        _PhyNameToElmtIndexList.push_back(make_pair(blockname,tempconn));
        _PhyNameToElmtIndexList.push_back(make_pair("alldomain",tempconn));

        _BulkPhyGroupNameList.push_back("alldomain");

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
            if(_ElmtDimVec[e]==_nDim){
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
        _PhyGroupNameList.clear();
        _PhyIDToNameList.clear();
        _PhyNameToIDList.clear();
        _PhyNameToElmtIndexList.clear();
        _BounPhyGroupNameList.clear();
        _BulkPhyGroupNameList.clear();
        for(int j=0;j<int(phyconn.size());++j){
            _PhyGroupNameList.push_back(_PhyGroupArray[j].second);
            _PhyIDToNameList.push_back(make_pair(_PhyGroupArray[j].first,_PhyGroupArray[j].second));
            _PhyNameToIDList.push_back(make_pair(_PhyGroupArray[j].second,_PhyGroupArray[j].first));
            _PhyNameToElmtIndexList.push_back(make_pair(_PhyGroupArray[j].second,phyconn[j]));
            if(_PhyGroupDimVec[j]==_nDim){
                _BulkPhyGroupNameList.push_back(_PhyGroupArray[j].second);
            }
            
            if(_nDim==1){
                if(_PhyGroupDimVec[j]==0){
                    _BounPhyGroupNameList.push_back(_PhyGroupArray[j].second);
                }
            }
            else if(_nDim==2){
                if(_PhyGroupDimVec[j]<=1){
                    _BounPhyGroupNameList.push_back(_PhyGroupArray[j].second);
                }
            }
            else if(_nDim==3){
                if(_PhyGroupDimVec[j]<=2){
                    _BounPhyGroupNameList.push_back(_PhyGroupArray[j].second);
                }
            }
        }
        _PhyGroupNameList.push_back("alldomain");
        _PhyIDToNameList.push_back(make_pair(maxid+1,"alldomain"));
        _PhyNameToIDList.push_back(make_pair("alldomain",maxid+1));
        _PhyNameToElmtIndexList.push_back(make_pair("alldomain",tempconn));

        _BulkPhyGroupNameList.push_back("alldomain");

        return (_nBulkElmts>0);
    }
    return false;
}