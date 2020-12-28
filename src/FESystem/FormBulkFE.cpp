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
//+++ Date   : 2020.11.30
//+++ Purpose: Loop all the bulk elements, where we can calculate
//+++          residual, k matrix as well as projection quantities
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FESystem/FESystem.h"

void FESystem::FormBulkFE(const FECalcType &calctype,const double &t,const double &dt,const double (&ctan)[2],
                Mesh &mesh,const DofHandler &dofHandler,FE &fe,
                ElmtSystem &elmtSystem,MateSystem &mateSystem,
                const Vec &U,const Vec &V,
                Vec &Hist,const Vec &HistOld,Vec &Proj,
                Mat &AMATRIX,Vec &RHS){
    
    if(calctype==FECalcType::ComputeResidual){
        VecSet(RHS,0.0);
    }
    else if(calctype==FECalcType::ComputeJacobian){
        MatZeroEntries(AMATRIX);
    }
    else if(calctype==FECalcType::InitHistoryVariable){
        VecSet(Hist,0.0);
    }
    else if(calctype==FECalcType::Projection){
        VecSet(Proj,0.0);
    }
    else{
        MessagePrinter::PrintErrorTxt("unsupported calculation type in FormBulkFE, please check your code");
        MessagePrinter::AsFem_Exit();
    }

    MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&_size);

    // we can get the correct value on the ghosted node!
    VecScatterCreateToAll(U,&_scatteru,&_Useq);
    VecScatterBegin(_scatteru,U,_Useq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatteru,U,_Useq,INSERT_VALUES,SCATTER_FORWARD);

    VecScatterCreateToAll(V,&_scatterv,&_Vseq);
    VecScatterBegin(_scatterv,V,_Vseq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterv,V,_Vseq,INSERT_VALUES,SCATTER_FORWARD);
    // //*** for proj and hist variables
    // VecScatterCreateToAll(Proj,&_scatterproj,&_ProjSeq);
    // VecScatterBegin(_scatterproj,Proj,_ProjSeq,INSERT_VALUES,SCATTER_FORWARD);
    // VecScatterEnd(_scatterproj,Proj,_ProjSeq,INSERT_VALUES,SCATTER_FORWARD);

    VecScatterCreateToAll(Hist,&_scatterhist,&_HistSeq);
    VecScatterBegin(_scatterhist,Hist,_HistSeq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterhist,Hist,_HistSeq,INSERT_VALUES,SCATTER_FORWARD);

    VecScatterCreateToAll(HistOld,&_scatterhistold,&_HistOldSeq);
    VecScatterBegin(_scatterhistold,HistOld,_HistOldSeq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterhistold,HistOld,_HistOldSeq,INSERT_VALUES,SCATTER_FORWARD);


    int rankne=mesh.GetBulkMeshBulkElmtsNum()/_size;
    int eStart=_rank*rankne;
    int eEnd=(_rank+1)*rankne;
    if(_rank==_size-1) eEnd=mesh.GetBulkMeshBulkElmtsNum();

    PetscInt nDofs,nNodes,nDofsPerNode,nDofsPerSubElmt,e;
    PetscInt i,j,jj;
    PetscInt nDim,gpInd;
    PetscReal xi,eta,zeta,w,JxW,DetJac,elVolume;
    nDim=mesh.GetDim();

    _BulkVolumes=0.0;
    
    for(int ee=eStart;ee<eEnd;++ee){
        e=ee+1;
        mesh.GetBulkMeshIthBulkElmtNodes(e,_elNodes);
        mesh.GetBulkMeshIthBulkElmtConn(e,_elConn);
        dofHandler.GetIthBulkElmtDofIndex0(e,_elDofs,_elDofsActiveFlag);
        nDofs=dofHandler.GetIthBulkElmtDofsNum(e);
        nNodes=mesh.GetBulkMeshIthBulkElmtNodesNum(e);
        nDofsPerNode=nDofs/nNodes;
        
        VecGetValues(_Useq,nDofs,_elDofs.data(),_elU.data());
        VecGetValues(_Vseq,nDofs,_elDofs.data(),_elV.data());
        
        if(calctype==FECalcType::ComputeResidual){
            fill(_R.begin(),_R.end(),0.0);
        }
        else if(calctype==FECalcType::ComputeJacobian){
            fill(_K.begin(),_K.end(),0.0);
        }
        
        xi=0.0;eta=0.0;zeta=0.0;DetJac=1.0;w=1.0;
        elVolume=0.0;
        for(gpInd=1;gpInd<=fe._BulkQPoint.GetQpPointsNum();++gpInd){
            // init all the local K&R array/matrix
            // get local history(old) value on each gauss point
            for(i=1;i<=_nHist;++i){
                _gpHist[i-1]=0.0;
                j=(e-1)*_nHist*fe._BulkQPoint.GetQpPointsNum()+(gpInd-1)*_nHist+i-1;
                VecGetValues(_HistOldSeq,1,&j,&_gpHistOld[i-1]);
            }
            // calculate the current shape funs on each gauss point
            if(nDim==1){
                w =fe._BulkQPoint.GetIthQpPointJthCoord(gpInd,0);
                xi=fe._BulkQPoint.GetIthQpPointJthCoord(gpInd,1);
                fe._BulkShp.Calc(xi,_elNodes,true);
                DetJac=fe._BulkShp.GetDetJac();
            }
            else if(nDim==2){
                w  =fe._BulkQPoint.GetIthQpPointJthCoord(gpInd,0);
                xi =fe._BulkQPoint.GetIthQpPointJthCoord(gpInd,1);
                eta=fe._BulkQPoint.GetIthQpPointJthCoord(gpInd,2);
                fe._BulkShp.Calc(xi,eta,_elNodes,true);
                DetJac=fe._BulkShp.GetDetJac();
            }
            else if(nDim==3){
                w   =fe._BulkQPoint.GetIthQpPointJthCoord(gpInd,0);
                xi  =fe._BulkQPoint.GetIthQpPointJthCoord(gpInd,1);
                eta =fe._BulkQPoint.GetIthQpPointJthCoord(gpInd,2);
                zeta=fe._BulkQPoint.GetIthQpPointJthCoord(gpInd,3);
                fe._BulkShp.Calc(xi,eta,zeta,_elNodes,true);
                DetJac=fe._BulkShp.GetDetJac();
            }
            JxW=w*DetJac;
            elVolume+=1.0*JxW;
            // calculate the coordinate of current gauss point
            _gpCoord(1)=0.0;_gpCoord(2)=0.0;_gpCoord(3)=0.0;
            for(i=1;i<=nNodes;++i){
                _gpCoord(1)+=_elNodes(i,1)*fe._BulkShp.shape_value(i);
                _gpCoord(2)+=_elNodes(i,2)*fe._BulkShp.shape_value(i);
                _gpCoord(3)+=_elNodes(i,3)*fe._BulkShp.shape_value(i);
            }
            
            _localK.setZero();_localR.setZero();
            // now we do the loop for local element, *local element could have multiple contributors according
            // to your model, i.e. one element (or one domain) can be assigned by multiple [elmt] sub block in your input file !!!
            for(int ielmt=1;ielmt<=static_cast<int>(dofHandler.GetIthElmtElmtMateTypePair(e).size());ielmt++){
                elmttype=dofHandler.GetIthElmtJthKernelElmtType(e,ielmt);
                matetype=dofHandler.GetIthElmtJthKernelMateType(e,ielmt);
                localDofIndex=dofHandler.GetIthBulkElmtJthKernelDofIndex(e,ielmt);
                mateindex=dofHandler.GetIthBulkElmtJthKernelMateIndex(e,ielmt);
                nDofsPerSubElmt=static_cast<int>(localDofIndex.size());
                
                // now we calculate the local dofs and their derivatives
                // *this is only the local one, which means, i.e., if current element use dofs=u v
                // then we only calculate u v and their derivatives on each gauss point,
                // then for the next loop, the same element may use 'dofs=u v w', then we will calculate
                // u v w and their derivatives, and so on!!!
                // In short, here we dont offer the localK and localR for the whole element, instead, the quantities
                // of a single gauss point according to each sub [elmt] block
                for(j=1;j<=nDofsPerSubElmt;j++){
                    // !!!: the index starts from 1, not 0, please following the same way in your UEL !!!
                    _gpU[j]=0.0;_gpV[j]=0.0;
                    _gpGradU[j](1)=0.0;_gpGradU[j](2)=0.0;_gpGradU[j](3)=0.0;
                    _gpGradV[j](1)=0.0;_gpGradV[j](2)=0.0;_gpGradV[j](3)=0.0;
                    jj=localDofIndex[j-1];
                    for(i=1;i<=nNodes;++i){
                        _gpU[j]+=_elU[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_value(i);
                        _gpV[j]+=_elV[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_value(i);

                        _gpGradU[j](1)+=_elU[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(1);
                        _gpGradU[j](2)+=_elU[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(2);
                        _gpGradU[j](3)+=_elU[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(3);

                        _gpGradV[j](1)+=_elV[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(1);
                        _gpGradV[j](2)+=_elV[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(2);
                        _gpGradV[j](3)+=_elV[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(3);
                    }
                }

                if(calctype==FECalcType::ComputeResidual){
                    _subR.setZero();
                }
                else if(calctype==FECalcType::ComputeJacobian){
                    _subK.setZero();
                }
                else if(calctype==FECalcType::Projection){
                    fill(_gpProj.begin(),_gpProj.end(),0.0);
                }

                //*****************************************************
                //*** For user material calculation(UMAT)
                //*****************************************************
                mateSystem.RunBulkMateLibs(matetype,mateindex,nDim,t,dt,_gpCoord,_gpU,_gpV,_gpGradU,_gpGradV,_gpHist,_gpHistOld);

                //*****************************************************
                //*** For user element calculation(UEL)
                //*****************************************************
                if(calctype==FECalcType::ComputeResidual){
                    for(i=1;i<=nNodes;i++){
                        elmtSystem.RunBulkElmtLibs(calctype,elmttype,nDim,nNodes,t,dt,ctan,
                            _gpCoord,_gpU,_gpV,_gpGradU,_gpGradV,
                            fe._BulkShp.shape_value(i),fe._BulkShp.shape_value(i),// for Residual, we only need test fun
                            fe._BulkShp.shape_grad(i),fe._BulkShp.shape_grad(i),
                            mateSystem.GetScalarMatePtr(),
                            mateSystem.GetVectorMatePtr(),
                            mateSystem.GetRank2MatePtr(),
                            mateSystem.GetRank4MatePtr(),
                            _gpHist,_gpHistOld,_gpProj,_subK,_subR);

                        AssembleSubResidualToLocalResidual(nDofsPerNode,nDofsPerSubElmt,i,_subR,_localR);
                    }
                }
                else if(calctype==FECalcType::ComputeJacobian){
                    for(i=1;i<=nNodes;i++){
                        for(j=1;j<=nNodes;j++){
                            elmtSystem.RunBulkElmtLibs(calctype,elmttype,nDim,nNodes,t,dt,ctan,
                            _gpCoord,_gpU,_gpV,_gpGradU,_gpGradV,
                            fe._BulkShp.shape_value(i),fe._BulkShp.shape_value(j),
                            fe._BulkShp.shape_grad(i),fe._BulkShp.shape_grad(j),
                            mateSystem.GetScalarMatePtr(),
                            mateSystem.GetVectorMatePtr(),
                            mateSystem.GetRank2MatePtr(),
                            mateSystem.GetRank4MatePtr(),
                            _gpHist,_gpHistOld,_gpProj,_subK,_subR);
                            
                            AssembleSubJacobianToLocalJacobian(nDofsPerNode,i,j,_subK,_localK);
                        }
                    }
                }
                else if(calctype==FECalcType::Projection){
                    for(i=1;i<=nNodes;i++){
                        elmtSystem.RunBulkElmtLibs(calctype,elmttype,nDim,nNodes,t,dt,ctan,
                            _gpCoord,_gpU,_gpV,_gpGradU,_gpGradV,
                            fe._BulkShp.shape_value(i),fe._BulkShp.shape_value(i),// for Residual, we only need test fun
                            fe._BulkShp.shape_grad(i),fe._BulkShp.shape_grad(i),
                            mateSystem.GetScalarMatePtr(),
                            mateSystem.GetVectorMatePtr(),
                            mateSystem.GetRank2MatePtr(),
                            mateSystem.GetRank4MatePtr(),
                            _gpHist,_gpHistOld,_gpProj,_subK,_subR);
                        
                        AssembleLocalProjectionToGlobal(_elConn[i-1]-1,_nProj,JxW,fe._BulkShp.shape_value(i),_gpProj,Proj);
                    }
                }
            }//=====> end-of-sub-element-loop

            //***********************************************
            //*** accumulate all the local contribution inside gauss loop
            if(calctype==FECalcType::ComputeResidual){
                AccumulateLocalResidual(nDofs,_elDofsActiveFlag,JxW,_localR,_R);
            }
            else if(calctype==FECalcType::ComputeJacobian){
                AccumulateLocalJacobian(nDofs,_elDofsActiveFlag,JxW,_localK,_K);
            }
        }//----->end of gauss point loop
        mesh.SetBulkMeshIthBulkElmtVolume(e,elVolume);
        _BulkVolumes+=elVolume;

        
        if(calctype==FECalcType::ComputeResidual){
            AssembleLocalResidualToGlobalResidual(nDofs,_elDofs,_R,RHS);
        }
        else if(calctype==FECalcType::ComputeJacobian){
            AssembleLocalJacobianToGlobalJacobian(nDofs,_elDofs,_K,AMATRIX);
        }
    }//------>end of element loop


    //********************************************************************
    //*** finish all the final assemble for different matrix and array
    //********************************************************************
    if(calctype==FECalcType::ComputeResidual){
        VecAssemblyBegin(RHS);
        VecAssemblyEnd(RHS);
    }
    else if(calctype==FECalcType::ComputeJacobian){
        MatAssemblyBegin(AMATRIX,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(AMATRIX,MAT_FINAL_ASSEMBLY);
    }
    else if(calctype==FECalcType::Projection){
        Projection(mesh.GetBulkMeshNodesNum(),_nProj,Proj);
    }
    else if(calctype==FECalcType::InitHistoryVariable||calctype==FECalcType::UpdateHistoryVariable){
        VecAssemblyBegin(Hist);
        VecAssemblyEnd(Hist);
    }



    // delete scatter
    VecScatterDestroy(&_scatteru);
    VecScatterDestroy(&_scatterv);
    // VecScatterDestroy(&_scatterproj);
    VecScatterDestroy(&_scatterhist);
    VecScatterDestroy(&_scatterhistold);
    VecDestroy(&_Useq);
    VecDestroy(&_Vseq);
    // VecDestroy(&_ProjSeq);
    VecDestroy(&_HistSeq);
    VecDestroy(&_HistOldSeq);
}