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
//+++ Date   : 2020.11.30
//+++ Purpose: Loop all the bulk elements, where we can calculate
//+++          residual, k matrix as well as projection quantities
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FESystem/FESystem.h"

void FESystem::FormBulkFE(const FECalcType &calctype,const double &t,const double &dt,const double (&ctan)[2],
                Mesh &mesh,const DofHandler &dofHandler,FE &fe,
                ElmtSystem &elmtSystem,MateSystem &mateSystem,
                SolutionSystem &solutionSystem,
                Mat &AMATRIX,Vec &RHS){
    
    if(calctype==FECalcType::ComputeResidual){
        VecSet(RHS,0.0);
    }
    else if(calctype==FECalcType::ComputeJacobian){
        MatZeroEntries(AMATRIX);
    }
    else if(calctype==FECalcType::Projection){
        VecSet(solutionSystem._Proj,0.0);
        VecSet(solutionSystem._ProjScalarMate,0.0);
        VecSet(solutionSystem._ProjVectorMate,0.0);
        VecSet(solutionSystem._ProjRank2Mate,0.0);
        VecSet(solutionSystem._ProjRank4Mate,0.0);
    }
    else if(calctype==FECalcType::UpdateHistoryVariable||
            calctype==FECalcType::InitHistoryVariable||
            calctype==FECalcType::InitMaterialAndProjection){
        // do nothing!
    }
    else{
        MessagePrinter::PrintErrorTxt("unsupported calculation type in FormBulkFE, please check your code");
        MessagePrinter::AsFem_Exit();
    }

    MPI_Comm_rank(PETSC_COMM_WORLD,&_rank);
    MPI_Comm_size(PETSC_COMM_WORLD,&_size);

    // we can get the correct value on the ghosted node!
    // please keep in mind, we will always use Unew and V in SNES and TS !!!
    VecScatterCreateToAll(solutionSystem._Unew,&_scatteru,&_Useq);
    VecScatterBegin(_scatteru,solutionSystem._Unew,_Useq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatteru,solutionSystem._Unew,_Useq,INSERT_VALUES,SCATTER_FORWARD);


    VecScatterCreateToAll(solutionSystem._V,&_scatterv,&_Vseq);
    VecScatterBegin(_scatterv,solutionSystem._V,_Vseq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterv,solutionSystem._V,_Vseq,INSERT_VALUES,SCATTER_FORWARD);


    // for the disp and velocity in the previous step
    VecScatterCreateToAll(solutionSystem._Uold,&_scatteruold,&_Uoldseq);
    VecScatterBegin(_scatteruold,solutionSystem._Uold,_Uoldseq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatteruold,solutionSystem._Uold,_Uoldseq,INSERT_VALUES,SCATTER_FORWARD);

    VecScatterCreateToAll(solutionSystem._Vold,&_scattervold,&_Voldseq);
    VecScatterBegin(_scattervold,solutionSystem._Vold,_Voldseq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scattervold,solutionSystem._Vold,_Voldseq,INSERT_VALUES,SCATTER_FORWARD);

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

        // for the disp and velocity in current time step
        VecGetValues(_Useq,nDofs,_elDofs.data(),_elU.data());
        VecGetValues(_Vseq,nDofs,_elDofs.data(),_elV.data());
        // for the disp and velocity in the previous time step
        VecGetValues(_Uoldseq,nDofs,_elDofs.data(),_elUold.data());
        VecGetValues(_Voldseq,nDofs,_elDofs.data(),_elVold.data());
        
        if(calctype==FECalcType::ComputeResidual){
            fill(_R.begin(),_R.end(),0.0);
        }
        else if(calctype==FECalcType::ComputeJacobian){
            fill(_K.begin(),_K.end(),0.0);
        }
        else if(calctype==FECalcType::Projection){
            for(auto &it:_gpProj) it.second=0.0;
        }
        
        xi=0.0;eta=0.0;zeta=0.0;DetJac=1.0;w=1.0;
        elVolume=0.0;
        for(gpInd=1;gpInd<=fe._BulkQPoint.GetQpPointsNum();++gpInd){
            // init all the local K&R array/matrix
            // get local history(old) value on each gauss point
            if(calctype!=FECalcType::InitMaterialAndProjection){
                // the scalar/vector/rank-2/rank-4 materials in MateSystem is used by each quadrature point, so it is only used for one single gauss point. The materials of the whole system is stored in solution's materials array!!!
                mateSystem.GetScalarMateOldPtr()=solutionSystem._ScalarMaterialsOld[(e-1)*fe._BulkQPoint.GetQpPointsNum()+gpInd-1];
                mateSystem.GetVectorMateOldPtr()=solutionSystem._VectorMaterialsOld[(e-1)*fe._BulkQPoint.GetQpPointsNum()+gpInd-1];
                mateSystem.GetRank2MateOldPtr()=solutionSystem._Rank2TensorMaterialsOld[(e-1)*fe._BulkQPoint.GetQpPointsNum()+gpInd-1];
                mateSystem.GetRank4MateOldPtr()=solutionSystem._Rank4TensorMaterialsOld[(e-1)*fe._BulkQPoint.GetQpPointsNum()+gpInd-1];
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
            
            
            if(calctype==FECalcType::ComputeResidual){
                _localR.setZero();
            }
            else if(calctype==FECalcType::ComputeJacobian){
                _localK.setZero();
            }
            else if(calctype==FECalcType::Projection){
                for(auto &it:_gpProj) it.second=0.0;
            }
            else if(calctype==FECalcType::InitHistoryVariable||calctype==FECalcType::UpdateHistoryVariable){
                for(auto &it:solutionSystem._ScalarMaterials[(e-1)*fe._BulkQPoint.GetQpPointsNum()+gpInd-1]) it.second=0.0;
                for(auto &it:solutionSystem._VectorMaterials[(e-1)*fe._BulkQPoint.GetQpPointsNum()+gpInd-1]) it.second=0.0;
                for(auto &it:solutionSystem._Rank2TensorMaterials[(e-1)*fe._BulkQPoint.GetQpPointsNum()+gpInd-1]) it.second=0.0;
                for(auto &it:solutionSystem._Rank4TensorMaterials[(e-1)*fe._BulkQPoint.GetQpPointsNum()+gpInd-1]) it.second=0.0;
            }
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
                    _gpU[j]=0.0;_gpV[j]=0.0;_gpUOld[j]=0.0;_gpVOld[j]=0.0;
                    _gpGradU[j](1)=0.0;_gpGradU[j](2)=0.0;_gpGradU[j](3)=0.0;
                    _gpGradUOld[j](1)=0.0;_gpGradUOld[j](2)=0.0;_gpGradUOld[j](3)=0.0;
                    _gpGradV[j](1)=0.0;_gpGradV[j](2)=0.0;_gpGradV[j](3)=0.0;
                    _gpGradVOld[j](1)=0.0;_gpGradVOld[j](2)=0.0;_gpGradVOld[j](3)=0.0;
                    jj=localDofIndex[j-1];
                    for(i=1;i<=nNodes;++i){
                        _gpU[j]+=_elU[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_value(i);
                        _gpUOld[j]+=_elUold[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_value(i);

                        _gpV[j]+=_elV[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_value(i);
                        _gpVOld[j]+=_elVold[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_value(i);

                        _gpGradU[j](1)+=_elU[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(1);
                        _gpGradU[j](2)+=_elU[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(2);
                        _gpGradU[j](3)+=_elU[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(3);

                        _gpGradUOld[j](1)+=_elUold[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(1);
                        _gpGradUOld[j](2)+=_elUold[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(2);
                        _gpGradUOld[j](3)+=_elUold[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(3);

                        _gpGradV[j](1)+=_elV[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(1);
                        _gpGradV[j](2)+=_elV[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(2);
                        _gpGradV[j](3)+=_elV[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(3);

                        _gpGradVOld[j](1)+=_elVold[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(1);
                        _gpGradVOld[j](2)+=_elVold[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(2);
                        _gpGradVOld[j](3)+=_elVold[(i-1)*nDofsPerNode+jj-1]*fe._BulkShp.shape_grad(i)(3);
                    }
                }
                // for local soln structure
                _elmtsoln.gpGradU=_gpGradU;
                _elmtsoln.gpGradUold=_gpGradUOld;
                _elmtsoln.gpGradV=_gpGradV;
                _elmtsoln.gpGradVold=_gpGradVOld;
                _elmtsoln.gpU=_gpU;
                _elmtsoln.gpUold=_gpUOld;
                _elmtsoln.gpV=_gpV;
                _elmtsoln.gpVold=_gpVOld;

                if(calctype==FECalcType::ComputeResidual){
                    _subR.setZero();
                }
                else if(calctype==FECalcType::ComputeJacobian){
                    _subK.setZero();
                }
                //*****************************************************
                //*** For user material calculation(UMAT)
                //*****************************************************
                if(calctype==FECalcType::InitHistoryVariable){
                    mateSystem.InitBulkMateLibs(matetype,mateindex,nDim,_gpCoord,_gpU,_gpV,_gpGradU,_gpGradV);
                }
                else{
                    mateSystem.RunBulkMateLibs(matetype,mateindex,nDim,t,dt,_gpCoord,_gpU,_gpUOld,_gpV,_gpVOld,
                                               _gpGradU,_gpGradUOld,_gpGradV,_gpGradVOld);
                }
                //*****************************************************
                //*** For user element calculation(UEL)
                //*****************************************************
                // we set up the local information data structure before we go into each UEL
                _elmtinfo.dt=dt;_elmtinfo.t=t;
                _elmtinfo.gpCoords=_gpCoord;
                _elmtinfo.nDim=nDim;
                _elmtinfo.nDofs=nDofsPerSubElmt;
                _elmtinfo.nNodes=nNodes;
                if(calctype==FECalcType::ComputeResidual){
                    for(i=1;i<=nNodes;i++){
                        // for local shape function
                        _elmtshp.test=fe._BulkShp.shape_value(i);
                        _elmtshp.trial=fe._BulkShp.shape_value(i);
                        _elmtshp.grad_test=fe._BulkShp.shape_grad(i);
                        _elmtshp.grad_trial=fe._BulkShp.shape_grad(i);
                        
                        elmtSystem.RunBulkElmtLibs(calctype,elmttype,ctan,_elmtinfo,_elmtsoln,_elmtshp,mateSystem.GetMaterialsPtr(),mateSystem.GetMaterialsOldPtr(),_gpProj,_subK,_subR);
                        AssembleSubResidualToLocalResidual(nDofsPerNode,nDofsPerSubElmt,i,_subR,_localR);
                    }
                }
                else if(calctype==FECalcType::ComputeJacobian){
                    for(i=1;i<=nNodes;i++){
                        for(j=1;j<=nNodes;j++){
                            // for local shape function
                            _elmtshp.test=fe._BulkShp.shape_value(i);
                            _elmtshp.trial=fe._BulkShp.shape_value(j);
                            _elmtshp.grad_test=fe._BulkShp.shape_grad(i);
                            _elmtshp.grad_trial=fe._BulkShp.shape_grad(j);

                            elmtSystem.RunBulkElmtLibs(calctype,elmttype,ctan,_elmtinfo,_elmtsoln,_elmtshp,mateSystem.GetMaterialsPtr(),mateSystem.GetMaterialsOldPtr(),_gpProj,_subK,_subR);
                            
                            AssembleSubJacobianToLocalJacobian(nDofsPerNode,i,j,_subK,_localK);
                        }
                    }
                }
                else if(calctype==FECalcType::Projection){
                    for(i=1;i<=nNodes;i++){
                        // for local shape function
                        _elmtshp.test=fe._BulkShp.shape_value(i);
                        _elmtshp.trial=fe._BulkShp.shape_value(i);
                        _elmtshp.grad_test=fe._BulkShp.shape_grad(i);
                        _elmtshp.grad_trial=fe._BulkShp.shape_grad(i);
                        
                        elmtSystem.RunBulkElmtLibs(calctype,elmttype,ctan,_elmtinfo,_elmtsoln,_elmtshp,mateSystem.GetMaterialsPtr(),mateSystem.GetMaterialsOldPtr(),_gpProj,_subK,_subR);
                    }
                    // here we should not assemble the local projection, because the JxW should not be accumulated
                    // inside the element-loop, but the gpProj should be.
                    // therefore, each sub element should use its own place of gpProj, in short, the gpProj is shared
                    // between different elements
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
            else if(calctype==FECalcType::Projection){
                AssembleLocalProjectionToGlobal(nNodes,JxW,fe._BulkShp,_gpProj,
                                                mateSystem.GetScalarMatePtr(),
                                                mateSystem.GetVectorMatePtr(),
                                                mateSystem.GetRank2MatePtr(),
                                                mateSystem.GetRank4MatePtr(),
                                                solutionSystem);
            }
            else if(calctype==FECalcType::InitHistoryVariable||calctype==FECalcType::UpdateHistoryVariable){
                AssembleSubHistToLocal(e,fe._BulkQPoint.GetQpPointsNum(),gpInd,mateSystem.GetMaterialsPtr(),solutionSystem);
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
        else if(calctype==FECalcType::InitHistoryVariable||calctype==FECalcType::UpdateHistoryVariable){
            AssembleLocalHistToGlobal(e,fe._BulkQPoint.GetQpPointsNum(),solutionSystem);
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
        Projection(mesh.GetBulkMeshNodesNum(),solutionSystem);
    }




    // delete scatter
    VecScatterDestroy(&_scatteru);
    VecScatterDestroy(&_scatteruold);
    VecScatterDestroy(&_scatterv);
    VecScatterDestroy(&_scattervold);

    VecDestroy(&_Useq);
    VecDestroy(&_Uoldseq);

    VecDestroy(&_Vseq);
    VecDestroy(&_Voldseq);

}
