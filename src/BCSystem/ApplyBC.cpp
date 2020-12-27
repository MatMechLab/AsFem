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
//+++ Date   : 2020.12.26
//+++ Purpose: here we apply the boundary conditions we defined in
//+++          each [bcs] sub blocks
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "BCSystem/BCSystem.h"
#include "DofHandler/DofHandler.h"

void BCSystem::ApplyBC(const Mesh &mesh,const DofHandler &dofHandler,FE &fe,const FECalcType &calctype,const double &t,const double (&ctan)[2],Vec &U,Mat &AMATRIX,Vec &RHS){
    double bcvalue;
    vector<string> bcnamelist;
    int DofIndex;
    if(ctan[0]){}
    for(auto it:_BCBlockList){
        bcvalue=it._BCValue;
        if(it._IsTimeDependent) bcvalue=it._BCValue*t;
        DofIndex=it._DofID;
        bcnamelist=it._BoundaryNameList;
        if(it._BCType==BCType::DIRICHLETBC){
            ApplyDirichletBC(mesh,dofHandler,calctype,DofIndex,bcvalue,bcnamelist,U,AMATRIX,RHS);
        }
        else if(it._BCType==BCType::NEUMANNBC){
            if(calctype==FECalcType::ComputeResidual){
                ApplyNeumannBC(mesh,dofHandler,fe,DofIndex,bcvalue,bcnamelist,RHS);
            }
        }
        else if(it._BCType==BCType::NULLBC){
            continue;
        }
        else{
            // for other type boundary conditions
            int rankne,eStart,eEnd;
            int e,ee,i,j,k,iInd,jInd,gpInd;
            double value;

            MPI_Comm_size(PETSC_COMM_WORLD,&_size);
            MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);

            // we can get the correct value on the ghosted node!
            VecScatterCreateToAll(U,&_scatteru,&_Useq);
            VecScatterBegin(_scatteru,U,_Useq,INSERT_VALUES,SCATTER_FORWARD);
            VecScatterEnd(_scatteru,U,_Useq,INSERT_VALUES,SCATTER_FORWARD);

            for(auto bcname:bcnamelist){

                rankne=mesh.GetBulkMeshElmtsNumViaPhysicalName(bcname)/_size;
                eStart=_rank*rankne;
                eEnd=(_rank+1)*rankne;
                if(_rank==_size-1) eEnd=mesh.GetBulkMeshElmtsNumViaPhysicalName(bcname);

                _nDim=mesh.GetBulkMeshDimViaPhyName(bcname);
                _nNodesPerBCElmt=mesh.GetBulkMeshNodesNumPerElmtViaPhysicalName(bcname);
                _normals=0.0;
                for(e=eStart;e<eEnd;++e){
                    ee=mesh.GetBulkMeshIthElmtIDViaPhyName(bcname,e+1);

                    if(_nDim==0){
                        // for point case,(bulk dim=1, bc dim=0)
                        for(i=1;i<=_nNodesPerBCElmt;++i){
                            j=mesh.GetBulkMeshIthElmtJthNodeID(ee,i);
                            iInd=dofHandler.GetIthNodeJthDofIndex(j,DofIndex)-1;
                            if(calctype==FECalcType::ComputeResidual){
                                VecSetValue(RHS,iInd,bcvalue,ADD_VALUES);
                            }
                        }
                    }
                    else if(_nDim==1){
                        // for line case (bulk dim>=2, bc dim=1)
                        mesh.GetBulkMeshIthElmtNodes(ee,_elNodes);
                        for(gpInd=1;gpInd<=fe._LineQPoint.GetQpPointsNum();++gpInd){
                            _xi=fe._LineQPoint(gpInd,1);
                            fe._LineShp.Calc(_xi,_elNodes,false);// we calculate the derivatives on local coordinate!!!
                            _JxW=fe._LineShp.GetDetJac()*fe._LineQPoint(gpInd,0);
                            _normals=0.0;
                            _xs[0][0]=0.0;// dx/dxi
                            _xs[1][0]=0.0;// dy/dxi
                            _xs[2][0]=0.0;// dz/dxi
                            for(i=1;i<=_nNodesPerBCElmt;++i){
                                _xs[0][0]+=fe._LineShp.shape_grad(i)(1)*_elNodes(i,1);
                                _xs[1][0]+=fe._LineShp.shape_grad(i)(1)*_elNodes(i,2);
                                // dont need the following line, now we dont know how to define
                                // a normal vector of 3 D line
                                //_xs[2][0]+=fe._shp_line.shape_grad(i)(1)*_elNodes(i,2);
                                // cout<<fe._shp_line.shape_grad(i)(1)<<" ";
                            }
                            _dist=sqrt(_xs[0][0]*_xs[0][0]+_xs[1][0]*_xs[1][0]);
                            _normals(1)= _xs[1][0]/_dist;// dy/dxi
                            _normals(2)=-_xs[0][0]/_dist;// dx/dxi
                            _normals(3)= 0.0;

                            //*****************************************************
                            //*** calculate the quantities on current gauss point
                            _xi=fe._LineQPoint(gpInd,1);
                            fe._LineShp.Calc(_xi,_elNodes,true);
                            _JxW=fe._LineShp.GetDetJac()*fe._LineQPoint(gpInd,0);
                            _gpU=0.0;_gpGradU.setZero();_gpCoord.setZero();
                            for(i=1;i<=_nNodesPerBCElmt;++i){
                                j=mesh.GetBulkMeshIthElmtJthNodeID(ee,i);
                                iInd=dofHandler.GetIthNodeJthDofIndex(j,DofIndex)-1;
                                VecGetValues(_Useq,1,&iInd,&value);

                                _gpU+=fe._SurfaceShp.shape_value(i)*bcvalue*_JxW;

                                _gpGradU(1)+=value*fe._LineShp.shape_grad(i)(1);
                                _gpGradU(2)+=value*fe._LineShp.shape_grad(i)(2);
                                _gpGradU(3)+=value*fe._LineShp.shape_grad(i)(3);

                                _gpCoord(1)+=_elNodes(i,1)*fe._LineShp.shape_value(i);
                                _gpCoord(2)+=_elNodes(i,2)*fe._LineShp.shape_value(i);
                                _gpCoord(3)+=_elNodes(i,3)*fe._LineShp.shape_value(i);
                            }

                            if(calctype==FECalcType::ComputeResidual){
                                for(i=1;i<=_nNodesPerBCElmt;++i){
                                    RunBCLibs(it._BCType,calctype,_normals,_gpU,_gpGradU,
                                    bcvalue,
                                    fe._LineShp.shape_value(i),fe._LineShp.shape_value(i),
                                    fe._LineShp.shape_grad(i),fe._LineShp.shape_grad(i),_localK,_localR);

                                    j=mesh.GetBulkMeshIthElmtJthNodeID(ee,i);
                                    iInd=dofHandler.GetIthNodeJthDofIndex(j,DofIndex)-1;

                                    VecSetValue(RHS,iInd,_localR*_JxW,ADD_VALUES);
                                }
                            }
                            else if(calctype==FECalcType::ComputeJacobian){
                                for(i=1;i<=_nNodesPerBCElmt;++i){
                                    k=mesh.GetBulkMeshIthElmtJthNodeID(ee,i);
                                    iInd=dofHandler.GetIthNodeJthDofIndex(k,DofIndex)-1;
                                    for(j=1;j<=_nNodesPerBCElmt;++j){
                                        RunBCLibs(it._BCType,calctype,_normals,_gpU,_gpGradU,
                                        bcvalue,
                                        fe._LineShp.shape_value(i),fe._LineShp.shape_value(j),
                                        fe._LineShp.shape_grad(i),fe._LineShp.shape_grad(j),_localK,_localR);

                                        k=mesh.GetBulkMeshIthElmtJthNodeID(ee,j);
                                        jInd=dofHandler.GetIthNodeJthDofIndex(k,DofIndex)-1;

                                        MatSetValue(AMATRIX,iInd,jInd,_localK*_JxW,ADD_VALUES);
                                    }
                                }
                            }
                        }
                    }
                    else if(_nDim==2){
                        // for surface case (bulk dim=3, bc dim=2)
                        mesh.GetBulkMeshIthElmtNodes(ee,_elNodes);
                        for(gpInd=1;gpInd<=fe._SurfaceQPoint.GetQpPointsNum();++gpInd){
                            _xi=fe._SurfaceQPoint(gpInd,1);
                            _eta=fe._SurfaceQPoint(gpInd,2);
                            fe._SurfaceShp.Calc(_xi,_eta,_elNodes,false);
                            _JxW=fe._SurfaceShp.GetDetJac()*fe._SurfaceQPoint(gpInd,0);
                            
                            _xs[0][0]=0.0;// dx/dxi
                            _xs[0][1]=0.0;// dx/deta

                            _xs[1][0]=0.0;// dy/dxi
                            _xs[1][1]=0.0;// dy/deta

                            _xs[2][0]=0.0;// dz/dxi
                            _xs[2][1]=0.0;// dz/deta
                            for(i=1;i<=_nNodesPerBCElmt;++i){
                                _xs[0][0]+=fe._SurfaceShp.shape_grad(i)(1)*_elNodes(i,1);
                                _xs[0][1]+=fe._SurfaceShp.shape_grad(i)(2)*_elNodes(i,1);
                        
                                _xs[1][0]+=fe._SurfaceShp.shape_grad(i)(1)*_elNodes(i,2);
                                _xs[1][1]+=fe._SurfaceShp.shape_grad(i)(2)*_elNodes(i,2);

                                _xs[2][0]+=fe._SurfaceShp.shape_grad(i)(1)*_elNodes(i,3);
                                _xs[2][1]+=fe._SurfaceShp.shape_grad(i)(2)*_elNodes(i,3);
                            }
                            _normals(1) = _xs[2-1][1-1]*_xs[3-1][2-1]-_xs[3-1][1-1]*_xs[2-1][2-1];
                            _normals(2) = _xs[3-1][1-1]*_xs[1-1][2-1]-_xs[1-1][1-1]*_xs[3-1][2-1];
                            _normals(3) = _xs[1-1][1-1]*_xs[2-1][2-1]-_xs[2-1][1-1]*_xs[1-1][2-1];

                            _dist=sqrt(_normals(1)*_normals(1)+_normals(2)*_normals(2)*_normals(2)+_normals(3)*_normals(3));

                            _normals(1)=_normals(1)/_dist;
                            _normals(2)=_normals(2)/_dist;
                            _normals(3)=_normals(3)/_dist;

                            _xi=fe._SurfaceQPoint(gpInd,1);
                            _eta=fe._SurfaceQPoint(gpInd,2);
                            fe._SurfaceShp.Calc(_xi,_eta,_elNodes,true);
                            _JxW=fe._SurfaceShp.GetDetJac()*fe._SurfaceQPoint(gpInd,0);
                            
                            //*****************************************************
                            //*** calculate the quantities on current gauss point
                            _xi=fe._LineQPoint(gpInd,1);
                            fe._LineShp.Calc(_xi,_elNodes,true);
                            _JxW=fe._SurfaceShp.GetDetJac()*fe._SurfaceQPoint(gpInd,0);
                            _gpU=0.0;_gpGradU.setZero();_gpCoord.setZero();
                            for(i=1;i<=_nNodesPerBCElmt;++i){
                                j=mesh.GetBulkMeshIthElmtJthNodeID(ee,i);
                                iInd=dofHandler.GetIthNodeJthDofIndex(j,DofIndex)-1;
                                VecGetValues(_Useq,1,&iInd,&value);

                                _gpU+=fe._SurfaceShp.shape_value(i)*bcvalue*_JxW;

                                _gpGradU(1)+=value*fe._SurfaceShp.shape_grad(i)(1);
                                _gpGradU(2)+=value*fe._SurfaceShp.shape_grad(i)(2);
                                _gpGradU(3)+=value*fe._SurfaceShp.shape_grad(i)(3);

                                _gpCoord(1)+=_elNodes(i,1)*fe._SurfaceShp.shape_value(i);
                                _gpCoord(2)+=_elNodes(i,2)*fe._SurfaceShp.shape_value(i);
                                _gpCoord(3)+=_elNodes(i,3)*fe._SurfaceShp.shape_value(i);
                            }
                            if(calctype==FECalcType::ComputeResidual){
                                for(i=1;i<=_nNodesPerBCElmt;++i){
                                    RunBCLibs(it._BCType,calctype,_normals,_gpU,_gpGradU,
                                    bcvalue,
                                    fe._LineShp.shape_value(i),fe._LineShp.shape_value(i),
                                    fe._LineShp.shape_grad(i),fe._LineShp.shape_grad(i),_localK,_localR);

                                    j=mesh.GetBulkMeshIthElmtJthNodeID(ee,i);
                                    iInd=dofHandler.GetIthNodeJthDofIndex(j,DofIndex)-1;

                                    VecSetValue(RHS,iInd,_localR*_JxW,ADD_VALUES);
                                }
                            }
                            else if(calctype==FECalcType::ComputeJacobian){
                                for(i=1;i<=_nNodesPerBCElmt;++i){
                                    k=mesh.GetBulkMeshIthElmtJthNodeID(ee,i);
                                    iInd=dofHandler.GetIthNodeJthDofIndex(k,DofIndex)-1;
                                    for(j=1;j<=_nNodesPerBCElmt;++j){
                                        RunBCLibs(it._BCType,calctype,_normals,_gpU,_gpGradU,
                                        bcvalue,
                                        fe._SurfaceShp.shape_value(i),fe._SurfaceShp.shape_value(j),
                                        fe._SurfaceShp.shape_grad(i),fe._SurfaceShp.shape_grad(j),_localK,_localR);

                                        k=mesh.GetBulkMeshIthElmtJthNodeID(ee,j);
                                        jInd=dofHandler.GetIthNodeJthDofIndex(k,DofIndex)-1;

                                        MatSetValue(AMATRIX,iInd,jInd,_localK*_JxW,ADD_VALUES);
                                    }
                                }
                            }
                        }
                    }
                }

            }//===> end-of-boundary-name-list-loop
        }
    }
}
//****************************************************
void BCSystem::ApplyInitialBC(const Mesh &mesh,const DofHandler &dofHandler,const double &t,Vec &U){
    PetscReal bcvalue;
    vector<string> bcnamelist;
    PetscInt DofIndex;
    PetscInt i,j,e,ee,iInd;
    int rankne,eStart,eEnd;

    MPI_Comm_size(PETSC_COMM_WORLD,&_size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);
    
    for(auto it:_BCBlockList){
        bcvalue=it._BCValue;
        if(it._IsTimeDependent) bcvalue=t*it._BCValue;
        bcnamelist=it._BoundaryNameList;
        DofIndex=it._DofID;
        if(it._BCType==BCType::DIRICHLETBC){
            for(auto bcname:bcnamelist){
                rankne=mesh.GetBulkMeshElmtsNumViaPhysicalName(bcname)/_size;
                eStart=_rank*rankne;
                eEnd=(_rank+1)*rankne;
                if(_rank==_size-1) eEnd=mesh.GetBulkMeshElmtsNumViaPhysicalName(bcname);
                for(e=eStart;e<eEnd;++e){
                    ee=mesh.GetBulkMeshIthElmtIDViaPhyName(bcname,e+1);
                    for(i=1;i<=mesh.GetBulkMeshIthElmtNodesNumViaPhyName(bcname,e+1);++i){
                        j=mesh.GetBulkMeshIthElmtJthNodeID(ee,i);
                        iInd=dofHandler.GetIthNodeJthDofIndex(j,DofIndex)-1;
                        VecSetValues(U,1,&iInd,&bcvalue,INSERT_VALUES);
                    }
                }
            }
        }
        else{
            continue;
        }
    }

    VecAssemblyBegin(U);
    VecAssemblyEnd(U);
}