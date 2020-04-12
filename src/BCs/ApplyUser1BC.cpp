//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "BCs/BCSystem.h"

void BCSystem::ApplyUser1BC(Mesh &mesh,DofHandler &dofHandler,FE &fe,
            const PetscInt &DofIndex,const PetscReal &bcvalue,vector<string> bclist,const PetscReal (&ctan)[2],
            Mat &K,Vec &RHS,Vec &U){
    PetscInt i,j,e,ee,gpInd;
    PetscInt iInd;
    PetscScalar value,uqp;
    string bcname;
    int rankne,eStart,eEnd;

    MPI_Comm_size(PETSC_COMM_WORLD,&_size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);

    // we can get the correct value on the ghosted node!
    VecScatterCreateToAll(U,&_scatteru,&_Useq);
    VecScatterBegin(_scatteru,U,_Useq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatteru,U,_Useq,INSERT_VALUES,SCATTER_FORWARD);

    for(unsigned int ibc=0;ibc<bclist.size();++ibc){
        bcname=bclist[ibc];

        rankne=mesh.GetElmtsNumViaPhyName(bcname)/_size;
        eStart=_rank*rankne;
        eEnd=(_rank+1)*rankne;
        if(_rank==_size-1) eEnd=mesh.GetElmtsNumViaPhyName(bcname);

        _nDim=mesh.GetDimViaPhyName(bcname);
        _nNodesPerBCElmt=mesh.GetBCElmtNodesNumViaPhyName(bcname);

        // cout<<"bcname="<<bcname<<", nDim="<<_nDim<<", nNodesPerBCElmt="<<_nNodesPerBCElmt<<endl;

        for(e=eStart;e<eEnd;++e){
            ee=mesh.GetIthElmtIndexViaPhyName(bcname,e+1);
            // cout<<"bcname="<<bcname<<", bcvalue="<<bcvalue<<", e="<<ee<<":eEnd="<<eEnd<<endl;
            // for(i=1;i<=mesh.GetIthElmtNodesNumViaPhyName(bcname,e+1);++i){
            //     j=mesh.GetIthElmtJthConn(ee,i);
            //     iInd=dofHandler.GetIthNodeJthDofIndex(j,DofIndex)-1;
            //     VecSetValues(U,1,&iInd,&bcvalue,INSERT_VALUES);
            //     VecSetValues(RHS,1,&iInd,&fix,INSERT_VALUES);
            //     MatSetValues(K,1,&iInd,1,&iInd,&_PenaltyFactor,INSERT_VALUES);
            //     // cout<<iInd+1<<" ";
            // }
            // cout<<endl;
            if(_nDim==0){
                for(i=1;i<=_nNodesPerBCElmt;++i){
                    j=mesh.GetIthElmtJthConn(ee,i);
                    iInd=dofHandler.GetIthNodeJthDofIndex(j,DofIndex)-1;
                    VecGetValues(_Useq,1,&iInd,&uqp);
                    if(uqp<0.9){
                        VecSetValue(RHS,iInd,bcvalue,ADD_VALUES);
                    }
                }
            }
            else if(_nDim==1){
                mesh.GetIthElmtNodes(ee,_elNodes);
                // cout<<"e="<<ee<<":";
                for(gpInd=1;gpInd<=fe._qp_line.GetQpPointsNum();++gpInd){
                    _xi=fe._qp_line(gpInd,1);
                    fe._shp_line.Calc(_xi,_elNodes,true);
                    _JxW=fe._shp_line.GetDetJac()*fe._qp_line(gpInd,0);
                    // cout<<"gp="<<gpInd<<", xi="<<xi<<", w="<<fe._qp_line(gpInd,0)<<",J="<<fe._shp_line.GetDetJac()<<endl;
                    for(i=1;i<=_nNodesPerBCElmt;++i){
                        j=mesh.GetIthElmtJthConn(ee,i);
                        iInd=dofHandler.GetIthNodeJthDofIndex(j,DofIndex)-1;
                        VecGetValues(_Useq,1,&iInd,&uqp);
                        if(uqp<0.9){
                            value=fe._shp_line.shape_value(i)*bcvalue*_JxW;
                            VecSetValue(RHS,iInd,value,ADD_VALUES);   
                        }
                        // cout<<"\tiInd="<<iInd<<",x="<<_elNodes(i,1)<<",y="<<_elNodes(i,2)<<",";
                        // cout<<"j="<<j<<", gp value="<<fe._shp_line.shape_value(i)*bcvalue*JxW<<",";
                        // cout<<fe._shp_line.shape_value(i)<<endl;
                    }
                }
            }
            else if(_nDim==2){
                mesh.GetIthElmtNodes(ee,_elNodes);
                for(gpInd=1;gpInd<=fe._qp_surface.GetQpPointsNum();++gpInd){
                    _xi=fe._qp_surface(gpInd,1);
                    _eta=fe._qp_surface(gpInd,2);
                    fe._shp_surface.Calc(_xi,_eta,_elNodes,true);
                    _JxW=fe._shp_surface.GetDetJac()*fe._qp_surface(gpInd,0);
                    for(i=1;i<=_nNodesPerBCElmt;++i){
                        j=mesh.GetIthElmtJthConn(ee,i);
                        iInd=dofHandler.GetIthNodeJthDofIndex(j,DofIndex)-1;
                        VecGetValues(_Useq,1,&iInd,&uqp);
                        if(uqp<0.9){
                            value=fe._shp_surface.shape_value(i)*bcvalue*_JxW;
                            VecSetValue(RHS,iInd,value,ADD_VALUES);
                        }
                    }
                }
            }
        }
    }
    if(ctan[0]){}
    MatGetSize(K,&iInd,&iInd);
}