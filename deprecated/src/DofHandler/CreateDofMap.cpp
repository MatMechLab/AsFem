#include "DofHandler/DofHandler.h"

void DofHandler::CreateDofMap(Mesh &mesh,BCSystem &bcSystem)
{
    int iblock;
    int DofIndex;
    string bcname,bctype;
    BCBlock bcBlock;
    int i,ii,j,k,e,ee,iInd,nDim;
    int nNodesPerBCElmt=0;


    _NodalDofFlag.resize(_nNodes,vector<int>(_nDofsPerNode,0));
    _NodalDofActiveFlag.resize(_nNodes,vector<double>(_nDofsPerNode,0.0));
    _ElmtalDofIndex.resize(_nBulkElmts,vector<int>(_nDofsPerNode*_nNodesPerBulkElmt,0));
    _ElmtalDofActiveFlags.resize(_nBulkElmts,vector<double>(_nDofsPerNode*_nNodesPerBulkElmt,-1.0));
    


    _nActiveDofs=0;
    


    vector<int> elConn,elDofsIndex;
    elConn.resize(mesh.GetNodesNumPerBulkElmt(),0);

    for(e=1;e<=mesh.GetBulkElmtsNum();++e){
        mesh.GetIthBulkElmtConn(e,elConn);
        elDofsIndex=mesh.GetIthBulkElmtDofsIndex(e);
        for(i=1;i<=mesh.GetNodesNumPerBulkElmt();++i){
            j=elConn[i-1];
            for(k=1;k<=(int)elDofsIndex.size();++k){
                _NodalDofFlag[j-1][elDofsIndex[k-1]-1]=1;
                _NodalDofActiveFlag[j-1][elDofsIndex[k-1]-1]=1.0;
            }
        }
    }
    

    // now we start to count for the active dofs
    _nActiveDofs=0;
    for(i=0;i<_nNodes;++i)
    {
        for(j=0;j<_nDofsPerNode;++j)
        {
            if(_NodalDofFlag[i][j]!=0)
            {
                _nActiveDofs+=1;
                _NodalDofFlag[i][j]=_nActiveDofs;
                _NodalDofActiveFlag[i][j]=1.0;
            }
            else{
                continue;
            }
        }
    }
    // modify the nodal dof active flags according to related boundary condition
    nDim=mesh.GetDim();
    if(nDim==1){
        nNodesPerBCElmt=1;
    }
    else if(nDim==2){
        nNodesPerBCElmt=mesh.GetNodesNumPerLineElmt();
    }
    else if(nDim==3){
        nNodesPerBCElmt=mesh.GetNodesNumPerSurfaceElmt();
    }
    for(iblock=1;iblock<=bcSystem.GetBCBlockNum();++iblock){
        bcBlock=bcSystem.GetIthBCBlock(iblock);
        DofIndex=GetDofIndexViaName(bcBlock._DofName);
        
        bctype=bcBlock._BCElmtName;
        bcname=bcBlock._BoundaryName;
        //cout<<"bcname="<<bcname<<", dofindex="<<DofIndex<<endl; 
        if(bctype=="dirichlet"){
            for(e=1;e<=mesh.GetElmtsNumViaPhyName(bcname);++e){
                ee=mesh.GetElmtIndexViaPhyName(bcname,e);
                for(i=1;i<=nNodesPerBCElmt;++i){
                    j=mesh.GetIthElmtJthConn(ee,i);
                    // iInd=GetIthNodeJthDofIndex(j,DofIndex)-1;// nodal's dof index(global one)
                    // now we need to find the related local one
                    _NodalDofActiveFlag[j-1][DofIndex-1]=0.0;
                    // cout<<"j="<<j<<", iInd="<<iInd<<endl;
                }
            }
        }
        else {
            continue;
        }
    }

    // this has worse performance!
    // for(auto &it:_NodalDofFlag){
    //     for(j=0;j<_nDofsPerNode;++j){
    //         if(it[j]!=0){
    //             _nActiveDofs+=1;
    //             it[j]=_nActiveDofs;
    //         }
    //         else{
    //             continue;
    //         }
    //     }
    // }

    // int iInd,ii;
    
    for(e=0;e<_nBulkElmts;++e){
        for(j=1;j<=GetNodesNumPerBulkElmt();++j){
            iInd=mesh.GetIthBulkElmtJthConn(e+1,j);
            for(k=1;k<=GetDofsNumPerNode();++k){
                ii=(j-1)*GetDofsNumPerNode()+k-1;
                
                if(_NodalDofFlag[iInd-1][k-1]!=0){
                    _ElmtalDofIndex[e][ii]=_NodalDofFlag[iInd-1][k-1];
                }
                
                if(_NodalDofActiveFlag[iInd-1][k-1]>0.0){
                    
                    _ElmtalDofActiveFlags[e][ii]=1.0;
                }
                else{
                    _ElmtalDofActiveFlags[e][ii]=0.0;
                }
            }
        }
        _ElmtalDofIndex[e].erase(
                remove(_ElmtalDofIndex[e].begin(),_ElmtalDofIndex[e].end(),0),
                _ElmtalDofIndex[e].end());
        _ElmtalDofIndex[e].shrink_to_fit();

        _ElmtalDofActiveFlags[e].erase(
                remove(_ElmtalDofActiveFlags[e].begin(),_ElmtalDofActiveFlags[e].end(),-1.0),
                _ElmtalDofActiveFlags[e].end());
        _ElmtalDofActiveFlags[e].shrink_to_fit();
    }
    
}