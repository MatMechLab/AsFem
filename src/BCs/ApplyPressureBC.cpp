//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#include "BCs/BCSystem.h"

void BCSystem::ApplyPressureBC(Mesh &mesh,DofHandler &dofHandler,FE &fe,
                    const PetscInt &DofIndex,const PetscReal &bcvalue,vector<string> bclist,
                    Vec &RHS){
    PetscInt i,j,e,ee,gpInd;
    PetscInt iInd;
    PetscScalar value=0.0;
    string bcname;
    int rankne,eStart,eEnd;

    // int eStart,eEnd;


    for(unsigned int ibc=0;ibc<bclist.size();++ibc){
        bcname=bclist[ibc];

        MPI_Comm_size(PETSC_COMM_WORLD,&_size);
        MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);

        rankne=mesh.GetElmtsNumViaPhyName(bcname)/_size;
        // eStart=_rank*rankne;
        // eEnd=(_rank+1)*rankne;
        // if(_rank==_size-1) eEnd=mesh.GetElmtsNumViaPhyName(bcname);
        if(rankne){}
        eStart=0;
        eEnd=mesh.GetElmtsNumViaPhyName(bcname);

        _nDim=mesh.GetDimViaPhyName(bcname);
        _nNodesPerBCElmt=mesh.GetBCElmtNodesNumViaPhyName(bcname);

        if(_rank==0){
            // TODO: make it parallel, the current version is quite stupid!!!
        for(e=eStart;e<eEnd;++e){
            ee=mesh.GetIthElmtIndexViaPhyName(bcname,e+1);
            _normals=0.0;
           
            if(_nDim==0){
                // for point case,(bulk dim=1, bc dim=0)
                for(i=1;i<=_nNodesPerBCElmt;++i){
                    j=mesh.GetIthElmtJthConn(ee,i);
                    iInd=dofHandler.GetIthNodeJthDofIndex(j,DofIndex)-1;
                    VecSetValue(RHS,iInd,bcvalue,ADD_VALUES);
                }
            }
            else if(_nDim==1){
                // for line case (bulk dim>=2, bc dim=1)
                mesh.GetIthElmtNodes(ee,_elNodes);
                // cout<<"e="<<ee<<":";
                for(gpInd=1;gpInd<=fe._qp_line.GetQpPointsNum();++gpInd){
                    _xi=fe._qp_line(gpInd,1);
                    fe._shp_line.Calc(_xi,_elNodes,false);
                    _JxW=fe._shp_line.GetDetJac()*fe._qp_line(gpInd,0);
                    _normals=0.0;
                    _xs[0][0]=0.0;// dx/dxi
                    _xs[1][0]=0.0;// dy/dxi
                    _xs[2][0]=0.0;// dz/dxi
                    for(i=1;i<=_nNodesPerBCElmt;++i){
                        _xs[0][0]+=fe._shp_line.shape_grad(i)(1)*_elNodes(i,1);
                        _xs[1][0]+=fe._shp_line.shape_grad(i)(1)*_elNodes(i,2);
                        // dont need the following line, now we dont know how to define
                        // a normal vector of 3 D line
                        //_xs[2][0]+=fe._shp_line.shape_grad(i)(1)*_elNodes(i,2);
                        // cout<<fe._shp_line.shape_grad(i)(1)<<" ";
                    }
                    // cout<<endl;
                    _dist=sqrt(_xs[0][0]*_xs[0][0]+_xs[1][0]*_xs[1][0]);
                    _normals(1)= _xs[1][0]/_dist;// dy/dxi
                    _normals(2)=-_xs[0][0]/_dist;// dx/dxi
                    _normals(3)= 0.0;            // dz/dxi==0    
                    // cout<<"normal="<<_normals(1)<<" "<<_normals(2)<<endl;
                    // cout<<"|normal|="<<sqrt(_normals(1)*_normals(1)+_normals(2)*_normals(2))<<endl;
                    
                    for(i=1;i<=_nNodesPerBCElmt;++i){
                        j=mesh.GetIthElmtJthConn(ee,i);
                        iInd=dofHandler.GetIthNodeJthDofIndex(j,DofIndex)-1;
                        // here Traction=P*normal
                        if(DofIndex==1){
                            // for Ux
                            // cout<<"Ux="<<bcvalue*_normals(1)<<endl;
                            value=fe._shp_line.shape_value(i)*bcvalue*_normals(1)*_JxW;
                            // cout<<"i="<<i<<","<<"Ux="<<value<<endl;
                            VecSetValue(RHS,iInd,value,ADD_VALUES);
                        }
                        else if(DofIndex==2){
                            // for Uy
                            // cout<<"Uy="<<bcvalue*_normals(2)<<endl;
                            value=fe._shp_line.shape_value(i)*bcvalue*_normals(2)*_JxW;
                            VecSetValue(RHS,iInd,value,ADD_VALUES);
                        }
                    }
                }
            }
            else if(_nDim==2){
                // for surface case (bulk dim=3, bc dim=2)
                mesh.GetIthElmtNodes(ee,_elNodes);
                for(gpInd=1;gpInd<=fe._qp_surface.GetQpPointsNum();++gpInd){
                    _xi=fe._qp_surface(gpInd,1);
                    _eta=fe._qp_surface(gpInd,2);
                    fe._shp_surface.Calc(_xi,_eta,_elNodes,false);
                    _JxW=fe._shp_surface.GetDetJac()*fe._qp_surface(gpInd,0);

                    _xs[0][0]=0.0;// dx/dxi
                    _xs[0][1]=0.0;// dx/deta

                    _xs[1][0]=0.0;// dy/dxi
                    _xs[1][1]=0.0;// dy/deta

                    _xs[2][0]=0.0;// dz/dxi
                    _xs[2][1]=0.0;// dz/deta
                    for(i=1;i<=_nNodesPerBCElmt;++i){
                        _xs[0][0]+=fe._shp_surface.shape_grad(i)(1)*_elNodes(i,1);
                        _xs[0][1]+=fe._shp_surface.shape_grad(i)(2)*_elNodes(i,1);
                        
                        _xs[1][0]+=fe._shp_surface.shape_grad(i)(1)*_elNodes(i,2);
                        _xs[1][1]+=fe._shp_surface.shape_grad(i)(2)*_elNodes(i,2);

                        _xs[2][0]+=fe._shp_surface.shape_grad(i)(1)*_elNodes(i,3);
                        _xs[2][1]+=fe._shp_surface.shape_grad(i)(2)*_elNodes(i,3);
                    }
                    _normals(1) = _xs[2-1][1-1]*_xs[3-1][2-1]-_xs[3-1][1-1]*_xs[2-1][2-1];
                    _normals(2) = _xs[3-1][1-1]*_xs[1-1][2-1]-_xs[1-1][1-1]*_xs[3-1][2-1];
                    _normals(3) = _xs[1-1][1-1]*_xs[2-1][2-1]-_xs[2-1][1-1]*_xs[1-1][2-1];

                    _dist=sqrt(_normals(1)*_normals(1)+_normals(2)*_normals(2)*_normals(2)+_normals(3)*_normals(3));

                    _normals(1)=_normals(1)/_dist;
                    _normals(2)=_normals(2)/_dist;
                    _normals(3)=_normals(3)/_dist;

                    for(i=1;i<=_nNodesPerBCElmt;++i){
                        j=mesh.GetIthElmtJthConn(ee,i);
                        iInd=dofHandler.GetIthNodeJthDofIndex(j,DofIndex)-1;
                        if(DofIndex==1){
                            value=fe._shp_surface.shape_value(i)*bcvalue*_normals(1)*_JxW;
                        }
                        else if(DofIndex==2){
                            value=fe._shp_surface.shape_value(i)*bcvalue*_normals(2)*_JxW;
                        }
                        else if(DofIndex==3){
                            value=fe._shp_surface.shape_value(i)*bcvalue*_normals(3)*_JxW;
                        }
                        VecSetValue(RHS,iInd,value,ADD_VALUES);
                    }
                }
            }
        }

        }//----------->end if
        // VecAssemblyBegin(RHS);
        // VecAssemblyEnd(RHS);
    }
}