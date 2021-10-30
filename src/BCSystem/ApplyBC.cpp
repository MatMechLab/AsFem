//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2021
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

void BCSystem::ApplyBC(const Mesh &mesh,const DofHandler &dofHandler,FE &fe,const FECalcType &calctype,const double &t,const double (&ctan)[3],Vec &U,Vec &V,Mat &AMATRIX,Vec &RHS){
    double bcvalue;
    vector<string> bcnamelist;
    vector<int> DofsIndex;
    if(ctan[0]){}

    _elmtinfo.t=t;
    _elmtinfo.dt=0.0;
    for(auto it:_BCBlockList){
        bcvalue=it._BCValue;
        if(it._IsTimeDependent) bcvalue=it._BCValue*t;
        DofsIndex=it._DofIDs;
        bcnamelist=it._BoundaryNameList;
        if(it._BCType==BCType::DIRICHLETBC||
           it._BCType==BCType::USER1DIRICHLETBC||
           it._BCType==BCType::USER2DIRICHLETBC||
           it._BCType==BCType::USER3DIRICHLETBC||
           it._BCType==BCType::USER4DIRICHLETBC||
           it._BCType==BCType::USER5DIRICHLETBC){
            ApplyDirichletBC(calctype,it._BCType,bcnamelist,DofsIndex,bcvalue,it._Parameters,mesh,dofHandler,U,AMATRIX,RHS);
        }
        else if(it._BCType==BCType::NODALDIRICHLETBC){
            ApplyNodalDirichletBC(calctype,it._BCType,bcnamelist,DofsIndex,bcvalue,it._Parameters,mesh,dofHandler,U,AMATRIX,RHS);
        }
        else if(it._BCType==BCType::NODALNEUMANNBC){
            if(calctype==FECalcType::ComputeResidual){
                ApplyNodalNeumannBC(mesh,dofHandler,fe,DofsIndex,bcvalue,bcnamelist,RHS);
            }
        }
        else if(it._BCType==BCType::NULLBC){
            continue;
        }
        else{
            // for other type boundary conditions
            int rankne,eStart,eEnd;
            int e,ee,i,j,ii,jj,k,iInd,jInd,gpInd,ki,kj;
            double value;
            vector<int> dofids; // start from 0, not 1 !!!

            dofids.resize(DofsIndex.size(),0);

            _elmtinfo.nDofs=static_cast<int>(DofsIndex.size());

            MPI_Comm_size(PETSC_COMM_WORLD,&_size);
            MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);

            // we can get the correct value on the ghosted node!
            VecScatterCreateToAll(U,&_scatteru,&_Useq);
            VecScatterBegin(_scatteru,U,_Useq,INSERT_VALUES,SCATTER_FORWARD);
            VecScatterEnd(_scatteru,U,_Useq,INSERT_VALUES,SCATTER_FORWARD);

            VecScatterCreateToAll(V,&_scatterv,&_Vseq);
            VecScatterBegin(_scatterv,V,_Vseq,INSERT_VALUES,SCATTER_FORWARD);
            VecScatterEnd(_scatterv,V,_Vseq,INSERT_VALUES,SCATTER_FORWARD);

            for(auto bcname:bcnamelist){
                rankne=mesh.GetBulkMeshElmtsNumViaPhysicalName(bcname)/_size;
                eStart=_rank*rankne;
                eEnd=(_rank+1)*rankne;
                if(_rank==_size-1) eEnd=mesh.GetBulkMeshElmtsNumViaPhysicalName(bcname);

                _nDim=mesh.GetBulkMeshDimViaPhyName(bcname);
                _nNodesPerBCElmt=mesh.GetBulkMeshNodesNumPerElmtViaPhysicalName(bcname);
                _elmtinfo.nDim=_nDim;
                _elmtinfo.nNodes=_nNodesPerBCElmt;

                _normals=0.0;
                for(e=eStart;e<eEnd;++e){
                    ee=mesh.GetBulkMeshIthElmtIDViaPhyName(bcname,e+1);
                    
                    if(_nDim==0){
                        // for point case,(bulk dim=1, bc dim=0)
                        for(i=1;i<=_nNodesPerBCElmt;++i){
                            j=mesh.GetBulkMeshIthElmtJthNodeID(ee,i);
                            _elmtinfo.gpCoords(1)=mesh.GetBulkMeshIthNodeJthCoord(j,1);
                            _elmtinfo.gpCoords(2)=mesh.GetBulkMeshIthNodeJthCoord(j,2);
                            _elmtinfo.gpCoords(3)=mesh.GetBulkMeshIthNodeJthCoord(j,3);
                            for(k=0;k<_elmtinfo.nDofs;k++){
                                iInd=dofHandler.GetIthNodeJthDofIndex(j,DofsIndex[k])-1;
                                dofids[k]=iInd;
                            }
                        }
                        for(k=1;k<=_elmtinfo.nDofs;k++){
                            _soln.gpU[k]=0.0;
                            _soln.gpGradU[k]=0.0;
                            _soln.gpV[k]=0.0;
                            _soln.gpGradV[k]=0.0;
                        }
                        for(i=1;i<=_nNodesPerBCElmt;++i){
                            j=mesh.GetBulkMeshIthElmtJthNodeID(ee,i);
                            for(k=1;k<=_elmtinfo.nDofs;k++){
                                iInd=dofHandler.GetIthNodeJthDofIndex(j,DofsIndex[k-1])-1;
                                dofids[k-1]=iInd;
                                VecGetValues(_Useq,1,&iInd,&value);
                                _soln.gpU[k]=value;
                                _soln.gpGradU[k]=0.0;
                            }
                        }
                        if(calctype==FECalcType::ComputeResidual){
                            for(i=1;i<=_nNodesPerBCElmt;++i){
                                _shp.test=1.0;
                                _shp.grad_test=0.0;_shp.grad_test(1)=1.0;

                                _shp.trial=0.0;
                                _shp.grad_trial=0.0;_shp.grad_trial(1)=1.0;
                                
                                RunBCLibs(calctype,it._BCType,bcvalue,it._Parameters,_normals,ctan,_elmtinfo,_soln,_shp,_localR,_localK);
                                // assemble the localR to the global one
                                j=mesh.GetBulkMeshIthElmtJthNodeID(ee,i);
                                for(k=0;k<_elmtinfo.nDofs;k++){
                                    iInd=dofHandler.GetIthNodeJthDofIndex(j,DofsIndex[k])-1;
                                    VecSetValue(RHS,iInd,_localR(k+1)*1.0,ADD_VALUES);
                                }
                            }
                        }
                        else if(calctype==FECalcType::ComputeJacobian){
                            for(i=1;i<=_nNodesPerBCElmt;++i){
                                _shp.test=1.0;
                                _shp.grad_test=0.0;_shp.grad_test(1)=1.0;
                                ii=mesh.GetBulkMeshIthElmtJthNodeID(ee,i);
                                for(j=1;j<=_nNodesPerBCElmt;++j){
                                    _shp.trial=1.0;
                                    _shp.grad_trial=0.0;_shp.grad_trial(1)=1.0;
                                    jj=mesh.GetBulkMeshIthElmtJthNodeID(ee,j);
                                    RunBCLibs(calctype,it._BCType,bcvalue,it._Parameters,_normals,ctan,_elmtinfo,_soln,_shp,_localR,_localK);
                                    // assemble localK to the global one
                                    for(ki=0;ki<_elmtinfo.nDofs;ki++){
                                        iInd=dofHandler.GetIthNodeJthDofIndex(ii,DofsIndex[ki])-1;
                                        for(kj=0;kj<_elmtinfo.nDofs;kj++){
                                            jInd=dofHandler.GetIthNodeJthDofIndex(jj,DofsIndex[kj])-1;
                                            MatSetValue(AMATRIX,iInd,jInd,_localK(ki+1,kj+1)*1.0,ADD_VALUES);
                                        }
                                    }//===> end-of-localK-assemble
                                }//==>end-of-J-index-loop
                            }//===>end-of-I-index-loop
                        }//===>end-of-computeJacobian
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
                            
                            for(k=1;k<=_elmtinfo.nDofs;k++){
                                _soln.gpU[k]=0.0;
                                _soln.gpGradU[k]=0.0;
                                _soln.gpV[k]=0.0;
                                _soln.gpGradV[k]=0.0;
                            }
                            // calculate the quantities on current gauss point
                            _elmtinfo.gpCoords.setZero();
                            for(i=1;i<=_nNodesPerBCElmt;++i){
                                j=mesh.GetBulkMeshIthElmtJthNodeID(ee,i);
                                for(k=1;k<=_elmtinfo.nDofs;k++){
                                    iInd=dofHandler.GetIthNodeJthDofIndex(j,DofsIndex[k-1])-1;
                                    dofids[k-1]=iInd;
                                    VecGetValues(_Useq,1,&iInd,&value);
                                    _soln.gpU[k]+=fe._SurfaceShp.shape_value(i)*value;
                                    
                                    _soln.gpGradU[k](1)+=value*fe._LineShp.shape_grad(i)(1);
                                    _soln.gpGradU[k](2)+=value*fe._LineShp.shape_grad(i)(2);
                                    _soln.gpGradU[k](3)+=value*fe._LineShp.shape_grad(i)(3);
                                    
                                    VecGetValues(_Vseq,1,&iInd,&value);
                                    _soln.gpV[k]+=fe._SurfaceShp.shape_value(i)*value;
                                    
                                    _soln.gpGradV[k](1)+=value*fe._LineShp.shape_grad(i)(1);
                                    _soln.gpGradV[k](2)+=value*fe._LineShp.shape_grad(i)(2);
                                    _soln.gpGradV[k](3)+=value*fe._LineShp.shape_grad(i)(3);

                                }
                                _elmtinfo.gpCoords(1)+=fe._LineShp.shape_value(i)*_elNodes(i,1);
                                _elmtinfo.gpCoords(2)+=fe._LineShp.shape_value(i)*_elNodes(i,2);
                                _elmtinfo.gpCoords(3)+=fe._LineShp.shape_value(i)*_elNodes(i,3);
                            }
                            
                            if(calctype==FECalcType::ComputeResidual){
                                for(i=1;i<=_nNodesPerBCElmt;++i){
                                    _shp.test=fe._LineShp.shape_value(i);
                                    _shp.grad_test=fe._LineShp.shape_grad(i);

                                    _shp.trial=fe._LineShp.shape_value(i);
                                    _shp.grad_trial=fe._LineShp.shape_grad(i);

                                    RunBCLibs(calctype,it._BCType,bcvalue,it._Parameters,_normals,ctan,_elmtinfo,_soln,_shp,_localR,_localK);
                                    // assemble the localR to the global one
                                    j=mesh.GetBulkMeshIthElmtJthNodeID(ee,i);
                                    for(k=0;k<_elmtinfo.nDofs;k++){
                                        iInd=dofHandler.GetIthNodeJthDofIndex(j,DofsIndex[k])-1;
                                        VecSetValue(RHS,iInd,_localR(k+1)*_JxW,ADD_VALUES);
                                    }
                                }
                            }
                            else if(calctype==FECalcType::ComputeJacobian){
                                for(i=1;i<=_nNodesPerBCElmt;++i){
                                    _shp.test=fe._LineShp.shape_value(i);
                                    _shp.grad_test=fe._LineShp.shape_grad(i);
                                    ii=mesh.GetBulkMeshIthElmtJthNodeID(ee,i);
                                    for(j=1;j<=_nNodesPerBCElmt;++j){
                                        _shp.trial=fe._LineShp.shape_value(j);
                                        _shp.grad_trial=fe._LineShp.shape_grad(j);
                                        jj=mesh.GetBulkMeshIthElmtJthNodeID(ee,j);
                                        RunBCLibs(calctype,it._BCType,bcvalue,it._Parameters,_normals,ctan,_elmtinfo,_soln,_shp,_localR,_localK);
                                        // assemble localK to the global one
                                        for(ki=0;ki<_elmtinfo.nDofs;ki++){
                                            iInd=dofHandler.GetIthNodeJthDofIndex(ii,DofsIndex[ki])-1;
                                            for(kj=0;kj<_elmtinfo.nDofs;kj++){
                                                jInd=dofHandler.GetIthNodeJthDofIndex(jj,DofsIndex[kj])-1;
                                                MatSetValue(AMATRIX,iInd,jInd,_localK(ki+1,kj+1)*_JxW,ADD_VALUES);
                                            }
                                        }
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
                            //*****************************************************
                            //*** calculate the quantities on current gauss point
                            _xi=fe._SurfaceQPoint(gpInd,1);
                            _eta=fe._SurfaceQPoint(gpInd,2);
                            fe._SurfaceShp.Calc(_xi,_eta,_elNodes,true);
                            _JxW=fe._SurfaceShp.GetDetJac()*fe._SurfaceQPoint(gpInd,0);

                            for(k=1;k<=_elmtinfo.nDofs;k++){
                                _soln.gpU[k]=0.0;
                                _soln.gpGradU[k]=0.0;
                                _soln.gpV[k]=0.0;
                                _soln.gpGradV[k]=0.0;
                            }
                            // calculate the solution and other quantities on current gauss point
                            _elmtinfo.gpCoords.setZero();
                            for(i=1;i<=_nNodesPerBCElmt;++i){
                                j=mesh.GetBulkMeshIthElmtJthNodeID(ee,i);
                                for(k=1;k<=_elmtinfo.nDofs;k++){
                                    iInd=dofHandler.GetIthNodeJthDofIndex(j,DofsIndex[k-1])-1;
                                    dofids[k-1]=iInd;
                                    
                                    VecGetValues(_Useq,1,&iInd,&value);
                                    _soln.gpU[k]+=fe._SurfaceShp.shape_value(i)*value;
                                    
                                    _soln.gpGradU[k](1)+=value*fe._LineShp.shape_grad(i)(1);
                                    _soln.gpGradU[k](2)+=value*fe._LineShp.shape_grad(i)(2);
                                    _soln.gpGradU[k](3)+=value*fe._LineShp.shape_grad(i)(3);
                                    
                                    VecGetValues(_Vseq,1,&iInd,&value);
                                    _soln.gpV[k]+=fe._SurfaceShp.shape_value(i)*value;
                                    
                                    _soln.gpGradV[k](1)+=value*fe._LineShp.shape_grad(i)(1);
                                    _soln.gpGradV[k](2)+=value*fe._LineShp.shape_grad(i)(2);
                                    _soln.gpGradV[k](3)+=value*fe._LineShp.shape_grad(i)(3);
                                }
                                _elmtinfo.gpCoords(1)+=fe._SurfaceShp.shape_value(i)*_elNodes(i,1);
                                _elmtinfo.gpCoords(2)+=fe._SurfaceShp.shape_value(i)*_elNodes(i,2);
                                _elmtinfo.gpCoords(3)+=fe._SurfaceShp.shape_value(i)*_elNodes(i,3);
                            }
                            
                            // assemble localR to the global one 
                            if(calctype==FECalcType::ComputeResidual){
                                for(i=1;i<=_nNodesPerBCElmt;++i){
                                    _shp.test=fe._SurfaceShp.shape_value(i);
                                    _shp.grad_test=fe._SurfaceShp.shape_grad(i);

                                    _shp.trial=fe._SurfaceShp.shape_value(i);
                                    _shp.grad_trial=fe._SurfaceShp.shape_grad(i);

                                    RunBCLibs(calctype,it._BCType,bcvalue,it._Parameters,_normals,ctan,_elmtinfo,_soln,_shp,_localR,_localK);
                                    // assemble the localR to the global one
                                    j=mesh.GetBulkMeshIthElmtJthNodeID(ee,i);
                                    for(k=0;k<_elmtinfo.nDofs;k++){
                                        iInd=dofHandler.GetIthNodeJthDofIndex(j,DofsIndex[k])-1;
                                        VecSetValue(RHS,iInd,_localR(k+1)*_JxW,ADD_VALUES);
                                    }
                                }
                            }
                            else if(calctype==FECalcType::ComputeJacobian){
                                for(i=1;i<=_nNodesPerBCElmt;++i){
                                    _shp.test=fe._SurfaceShp.shape_value(i);
                                    _shp.grad_test=fe._SurfaceShp.shape_grad(i);
                                    ii=mesh.GetBulkMeshIthElmtJthNodeID(ee,i);
                                    for(j=1;j<=_nNodesPerBCElmt;++j){
                                        _shp.trial=fe._SurfaceShp.shape_value(j);
                                        _shp.grad_trial=fe._SurfaceShp.shape_grad(j);
                                        jj=mesh.GetBulkMeshIthElmtJthNodeID(ee,j);
                                        RunBCLibs(calctype,it._BCType,bcvalue,it._Parameters,_normals,ctan,_elmtinfo,_soln,_shp,_localR,_localK);
                                        // assemble localK to the global one
                                        for(ki=0;ki<_elmtinfo.nDofs;ki++){
                                            iInd=dofHandler.GetIthNodeJthDofIndex(ii,DofsIndex[ki])-1;
                                            for(kj=0;kj<_elmtinfo.nDofs;kj++){
                                                jInd=dofHandler.GetIthNodeJthDofIndex(jj,DofsIndex[kj])-1;
                                                MatSetValue(AMATRIX,iInd,jInd,_localK(ki+1,kj+1)*_JxW,ADD_VALUES);
                                            }
                                        }//===> end-of-localK-assemble-loop
                                    }//===> end-of-local-J-node-loop
                                }//===> end-of-local-I-node-loop
                            }//===> end-of-jacobian-calculation
                        }//===> end-of-gauss-point-loop
                    }//===> end-of-dim=2-case-calculation
                }//===> end-of-boundary-element-loop

            }//===> end-of-boundary-name-list-loop

            // now we can destroy all the temporary array
            // delete scatter
            VecScatterDestroy(&_scatteru);
            VecScatterDestroy(&_scatterv);
            // delete Useq and Vseq
            VecDestroy(&_Useq);
            VecDestroy(&_Vseq);
        }//===> end-of-boundary-type-if-else-condition
    }//===> end-of-bcblock-loop

    VecAssemblyBegin(U);
    VecAssemblyEnd(U);
    VecAssemblyBegin(RHS);
    VecAssemblyEnd(RHS);

    MatAssemblyBegin(AMATRIX,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(AMATRIX,MAT_FINAL_ASSEMBLY);
}
//****************************************************
void BCSystem::ApplyInitialBC(const Mesh &mesh,const DofHandler &dofHandler,const double &t,Vec &U){
    PetscReal bcvalue;
    vector<string> bcnamelist;
    vector<int> DofsIndex;
    PetscInt i,j,e,ee,iInd;
    int rankne,eStart,eEnd;

    MPI_Comm_size(PETSC_COMM_WORLD,&_size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);
    
    for(auto it:_BCBlockList){
        bcvalue=it._BCValue;
        if(it._IsTimeDependent) bcvalue=t*it._BCValue;
        bcnamelist=it._BoundaryNameList;
        DofsIndex=it._DofIDs;
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
                        for(const auto id:DofsIndex){
                            iInd=dofHandler.GetIthNodeJthDofIndex(j,id)-1;
                            VecSetValues(U,1,&iInd,&bcvalue,INSERT_VALUES);
                        }
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
