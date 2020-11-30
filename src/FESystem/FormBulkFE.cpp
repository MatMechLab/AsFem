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
                Mesh &mesh,DofHandler &dofHandler,FE &fe,
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
    //*** for proj and hist variables
    VecScatterCreateToAll(Proj,&_scatterproj,&_ProjSeq);
    VecScatterBegin(_scatterproj,Proj,_ProjSeq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterproj,Proj,_ProjSeq,INSERT_VALUES,SCATTER_FORWARD);

    VecScatterCreateToAll(Hist,&_scatterhist,&_HistSeq);
    VecScatterBegin(_scatterhist,Hist,_HistSeq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterhist,Hist,_HistSeq,INSERT_VALUES,SCATTER_FORWARD);

    VecScatterCreateToAll(HistOld,&_scatterhistold,&_HistOldSeq);
    VecScatterBegin(_scatterhistold,HistOld,_HistOldSeq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(_scatterhistold,HistOld,_HistOldSeq,INSERT_VALUES,SCATTER_FORWARD);

    

    int rankne=mesh.GetBulkElmtsNum()/_size;
    int eStart=_rank*rankne;
    int eEnd=(_rank+1)*rankne;
    if(_rank==_size-1) eEnd=mesh.GetBulkElmtsNum();

    PetscInt nDofs,nNodes,nDofsPerNode,e;
    PetscInt i,j;
    PetscInt nDim,gpInd;
    PetscReal xi,eta,zeta,w,JxW,DetJac,elVolume;
    ElmtType iuel;MateType iumate;int mateindex;
    nDim=mesh.GetDim();

    for(int ee=eStart;ee<eEnd;++ee){
        e=ee+1;
        mesh.GetIthBulkElmtNodes(e,_elNodes);
        mesh.GetIthBulkElmtConn(e,_elConn);
        // dofHandler.GetIthElmtDofIndex(e,_elDofs);
        dofHandler.GetIthBulkElmtDofIndex0(e,_elDofs,_elDofsActiveFlag);
        nDofs=dofHandler.GetIthBulkElmtDofsNum(e);
        nNodes=mesh.GetIthBulkElmtNodesNum(e);
        nDofsPerNode=nDofs/nNodes;
        //PetscPrintf(PETSC_COMM_WORLD,"************* work before !!!\n");
        VecGetValues(_Useq,nDofs,_elDofs.data(),_elU.data());
        VecGetValues(_Vseq,nDofs,_elDofs.data(),_elV.data());
        //PetscPrintf(PETSC_COMM_WORLD,"************* work ater !!!\n");

        // get current element's uel ID and material ID
        // iuel=mesh.GetIthBulkElmtElmtID(e);
        // imate=mesh.GetIthBulkElmtMateID(e);
        // mateindex=mesh.GetIthBulkElmtMateIDIndex(e);

        // iuel=mesh.GetIthBulkElmtElmtType(e);
        // iumate=mesh.GetIthBulkElmtMateType(e);
        // mateindex=mesh.GetIthBulkElmtMateIndex(e);
        

        // if(iumate==MateType::LinearThermalMechanicsMate){
        //     cout<<"work, mateindex="<<mateindex<<endl;
        // }

        
        if(calctype==FECalcType::ComputeResidual){
            fill(_R.begin(),_R.end(),0.0);
        }
        else if(calctype==FECalcType::ComputeJacobian){
            fill(_K.begin(),_K.end(),0.0);
        }
        
        xi=0.0;eta=0.0;zeta=0.0;DetJac=1.0;w=1.0;
        elVolume=0.0;
        for(gpInd=1;gpInd<=fe._BulkQPoint.GetQpPointsNum();++gpInd){
            if(isw==3||isw==6){
                _localR.setZero();
                if(isw==6){
                    _localK.setZero();
                }
            }
            for(i=1;i<=_nHist;++i){
                _gpHist[i-1]=0.0;
                j=(e-1)*_nHist*fe._BulkQPoint.GetQpPointsNum()+(gpInd-1)*_nHist+i-1;
                VecGetValues(_HistOldSeq,1,&j,&_gpHistOld[i-1]);
            }

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
            _Volumes+=1.0*JxW;
            elVolume+=1.0*JxW;
            
            for(j=1;j<=nDofsPerNode;++j){
                _gpU[j-1]=0.0;_gpV[j-1]=0.0;
                _gpGradU[j-1](1)=0.0;_gpGradU[j-1](2)=0.0;_gpGradU[j-1](3)=0.0;
                _gpGradV[j-1](1)=0.0;_gpGradV[j-1](2)=0.0;_gpGradV[j-1](3)=0.0;
                _gpCoord(1)=0.0;_gpCoord(2)=0.0;_gpCoord(3)=0.0;
                for(i=1;i<=nNodes;++i){
                    _gpU[j-1]+=_elU[(i-1)*nDofsPerNode+j-1]*fe._shp_bulk.shape_value(i);
                    _gpV[j-1]+=_elV[(i-1)*nDofsPerNode+j-1]*fe._shp_bulk.shape_value(i);

                    _gpGradU[j-1](1)+=_elU[(i-1)*nDofsPerNode+j-1]*fe._shp_bulk.shape_grad(i)(1);
                    _gpGradU[j-1](2)+=_elU[(i-1)*nDofsPerNode+j-1]*fe._shp_bulk.shape_grad(i)(2);
                    _gpGradU[j-1](3)+=_elU[(i-1)*nDofsPerNode+j-1]*fe._shp_bulk.shape_grad(i)(3);

                    _gpGradV[j-1](1)+=_elV[(i-1)*nDofsPerNode+j-1]*fe._shp_bulk.shape_grad(i)(1);
                    _gpGradV[j-1](2)+=_elV[(i-1)*nDofsPerNode+j-1]*fe._shp_bulk.shape_grad(i)(2);
                    _gpGradV[j-1](3)+=_elV[(i-1)*nDofsPerNode+j-1]*fe._shp_bulk.shape_grad(i)(3);

                    _gpCoord(1)+=_elNodes(i,1)*fe._shp_bulk.shape_value(i);
                    _gpCoord(2)+=_elNodes(i,2)*fe._shp_bulk.shape_value(i);
                    _gpCoord(3)+=_elNodes(i,3)*fe._shp_bulk.shape_value(i);
                }
            }
            

            //*****************************************************
            //*** For user material calculation(UMAT)
            //*****************************************************
            mateSystem.RunMateLib(iumate,mateindex,nDim,t,dt,_gpCoord,_gpU,_gpV,_gpGradU,_gpGradV,_gpHist,_gpHistOld);
            
            //*****************************************************
            //*** For user element calculation(UEL)
            //*****************************************************
            elmtSystem.RunElmtLib(isw,iuel,nDim,nNodes,t,dt,ctan,_gpCoord,_gpU,_gpV,_gpGradU,_gpGradV,
                                  fe._shp_bulk,
                                  mateSystem._ScalarMaterials,
                                  mateSystem._VectorMaterials,
                                  mateSystem._Rank2Materials,
                                  mateSystem._Rank4Materials,
                                  _gpHist,_gpHistOld,_gpProj,_localK,_localR);
            
            if(isw==3||isw==6){
                // accumulate the local contribution
                for(i=1;i<=nDofs;++i){
                    _R[i-1]+=_localR(i)*_elDofsActiveFlag[i-1]*JxW*_KMatrixFactor;
                    if(isw==6){
                        for(j=1;j<=nDofs;++j){
                            _K[(i-1)*nDofs+j-1]+=_localK(i,j)*_elDofsActiveFlag[i-1]*_elDofsActiveFlag[j-1]*JxW*_KMatrixFactor;
                            if(abs(_localK(i,j))>_MaxKMatrixValue) _MaxKMatrixValue=abs(_localK(i,j));
                        }
                    }
                }
            }
            else if(isw==9){
                Projection(nNodes,Proj,_gpProj,fe._shp_bulk,JxW);
            }
            else if(isw==4||isw==8){
                // assemble local hist to global one(for gauss point)
                AssembleLocalHistToGlobal(e,gpInd,_gpHist,Hist);
            }
        }//----->end of gauss point loop
        mesh.SetIthBulkElmtVolume(e,elVolume);

        
        if(isw==3||isw==6){
            // for(i=1;i<=nDofs;++i){
            //     for(j=1;j<=nDofs;++j){
            //         PetscPrintf(PETSC_COMM_WORLD,"%14.6e ",_K[(i-1)*nDofs+j-1]);
            //     }
            //     PetscPrintf(PETSC_COMM_WORLD,"\n");
            // }
            // PetscPrintf(PETSC_COMM_WORLD,"\n\n");
            // assemble local to global contribution
            AssembleLocalToGlobal(isw,nDofs,_elDofs,_K,_R,AMATRIX,RHS);
        }
    }//------>end of element loop

    if(isw==3||isw==6){
        VecAssemblyBegin(RHS);
        VecAssemblyEnd(RHS);
        if(isw==6){
            MatAssemblyBegin(AMATRIX,MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(AMATRIX,MAT_FINAL_ASSEMBLY);
        }
    }
    else if(isw==4||isw==8){
        VecAssemblyBegin(Hist);
        VecAssemblyEnd(Hist);
    }
    else if(isw==9){
        VecAssemblyBegin(Proj);
        VecAssemblyEnd(Proj);
    }

    // delete scatter
    VecScatterDestroy(&_scatteru);
    VecScatterDestroy(&_scatterv);
    VecScatterDestroy(&_scatterproj);
    VecScatterDestroy(&_scatterhist);
    VecScatterDestroy(&_scatterhistold);
    VecDestroy(&_Useq);
    VecDestroy(&_Vseq);
    VecDestroy(&_ProjSeq);
    VecDestroy(&_HistSeq);
    VecDestroy(&_HistOldSeq);
}