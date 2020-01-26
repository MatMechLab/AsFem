#include "FESystem/FESystem.h"

void FESystem::FormFE(const int &isw,const double &t,const double &dt,const double (&ctan)[2],
                      Mesh &mesh,DofHandler &dofHandler,
                      ElmtSystem &elmtSystem,MaterialSystem &mateSystem,
                      const Eigen::VectorXd &U,const Eigen::VectorXd &V,
                      Eigen::MatrixXd &Hist,const Eigen::MatrixXd &HistOld,Eigen::MatrixXd &Proj,
                      Eigen::SparseMatrix<double,Eigen::RowMajor> &AMATRIX,Eigen::VectorXd &RHS){
    // loop over all the bulk elmts
    if(isw==3||isw==6){
        RHS.setZero();
        if(isw==6){
            AMATRIX*=0.0;
            // AMATRIX.setZero();
            // fill(AMATRIX.operator*=)
        }
    }
    else if(isw==4){
        // init history value
        Hist.setZero();
    }
    else if(isw==8)
    {
        Hist=HistOld;
    }
    else if(isw==9)
    {
        Proj.setZero();
    }
    int e,i,j,gpInd,nDofs,nNodes,nDofsPerNode;
    int nDim=mesh.GetDim();
    double xi,eta,zeta,w,JxW,DetJac;
    Eigen::Vector3d _gpCoord;
    MaterialType imate;
    int mateindex;
    ElmtType iuel;
    _Volumes=0.0;

    for(e=1;e<=mesh.GetBulkElmtsNum();++e){
        mesh.GetIthBulkElmtNodes(e,_elNodes);
        mesh.GetIthBulkElmtConn(e,_elConn);
        // dofHandler.GetIthElmtDofIndex(e,_elDofs);
        dofHandler.GetIthElmtDofIndex(e,_elDofs,_elDofsActiveFlag);
        nDofs=dofHandler.GetIthElmtDofsNum(e);
        nNodes=mesh.GetIthBulkElmtNodesNum(e);
        nDofsPerNode=nDofs/nNodes;

        // cout<<"*** e="<<e<<":"<<endl;
        // cout<<"***   coords=";
        // for(i=1;i<=nNodes;++i){
        //     if(i>1)cout<<"***          ";
        //     cout<<_elNodes(i,1)<<" "<<_elNodes(i,2)<<" "<<_elNodes(i,3)<<endl;
        // }
        // cout<<endl;
        // cout<<"***   ndofs="<<nDofs;
        // for(i=1;i<=nDofs;++i){
        //     cout<<_elDofsActiveFlag[i-1]<<" ";
        // }
        // cout<<endl;
        
        // cout<<"***   elU=";
        for(i=0;i<nDofs;++i){
            _elU[i]=U.coeff(_elDofs[i]-1);
            _elV[i]=V.coeff(_elDofs[i]-1);
            // cout<<_elU[i]<<" ";
        }
        // cout<<endl;
        
        // get current element's uel ID and material ID
        iuel=mesh.GetIthBulkElmtElmtID(e);
        imate=mesh.GetIthBulkElmtMateID(e);
        mateindex=mesh.GetIthBulkElmtMateIDIndex(e);
        
        //*******************************************************************
        // please, do not zero them inside the uel
        // this is very important, otherwise you will have trouble!!!
        //*******************************************************************
        if(isw==3||isw==6){
            _localR.setZero();
            if(isw==6){
                _localK.setZero();
            }
        }
        xi=0.0;eta=0.0;zeta=0.0;DetJac=1.0;w=1.0;
        for(gpInd=1;gpInd<=fe._qp_bulk.GetQpPointsNum();++gpInd){
            for(i=1;i<=_nHist;++i){
                _gpHist[i-1]=0.0;
                _gpHistOld[i-1]=HistOld.coeff(e-1,(gpInd-1)*_nHist+i-1);
            }

            if(nDim==1){
                w =fe._qp_bulk.GetIthQpPointJthCoord(gpInd,0);
                xi=fe._qp_bulk.GetIthQpPointJthCoord(gpInd,1);
                fe._shp_bulk.Calc(xi,_elNodes);
                DetJac=fe._shp_bulk.GetDetJac();
            }
            else if(nDim==2){
                w  =fe._qp_bulk.GetIthQpPointJthCoord(gpInd,0);
                xi =fe._qp_bulk.GetIthQpPointJthCoord(gpInd,1);
                eta=fe._qp_bulk.GetIthQpPointJthCoord(gpInd,2);
                fe._shp_bulk.Calc(xi,eta,_elNodes);
                DetJac=fe._shp_bulk.GetDetJac();
            }
            else if(nDim==3){
                w   =fe._qp_bulk.GetIthQpPointJthCoord(gpInd,0);
                xi  =fe._qp_bulk.GetIthQpPointJthCoord(gpInd,1);
                eta =fe._qp_bulk.GetIthQpPointJthCoord(gpInd,2);
                zeta=fe._qp_bulk.GetIthQpPointJthCoord(gpInd,3);
                fe._shp_bulk.Calc(xi,eta,zeta,_elNodes);
                DetJac=fe._shp_bulk.GetDetJac();
            }
            JxW=w*DetJac;
            _Volumes+=1.0*JxW;

            // cout<<"***qp="<<gpInd-1<<endl;
            // if(e==1){
            // cout<<"***qp="<<gpInd-1<<", w="<<w
            //     <<", xi="<<xi
            //     <<", eta="<<eta
            //     <<", zeta="<<zeta
            //     <<", DetJac="<<DetJac<<", JxW="<<JxW<<endl;
            // }
            
            
            // cout<<"***          gpU=";
            //*************************************************
            for(j=1;j<=nDofsPerNode;++j){
                _gpU[j-1]=0.0;_gpV[j-1]=0.0;
                _gpGradU[j-1].coeffRef(0)=0.0;_gpGradU[j-1].coeffRef(1)=0.0;_gpGradU[j-1].coeffRef(2)=0.0;
                _gpGradV[j-1].coeffRef(0)=0.0;_gpGradV[j-1].coeffRef(1)=0.0;_gpGradV[j-1].coeffRef(2)=0.0;
                _gpCoord.coeffRef(0)=0.0;_gpCoord.coeffRef(1)=0.0;_gpCoord.coeffRef(2)=0.0;
                for(i=1;i<=nNodes;++i){
                    _gpU[j-1]+=_elU[(i-1)*nDofsPerNode+j-1]*fe._shp_bulk.shape_value(i);
                    _gpV[j-1]+=_elV[(i-1)*nDofsPerNode+j-1]*fe._shp_bulk.shape_value(i);

                    _gpGradU[j-1].coeffRef(0)+=_elU[(i-1)*nDofsPerNode+j-1]*fe._shp_bulk.shape_grad(i).coeff(0);
                    _gpGradU[j-1].coeffRef(1)+=_elU[(i-1)*nDofsPerNode+j-1]*fe._shp_bulk.shape_grad(i).coeff(1);
                    _gpGradU[j-1].coeffRef(2)+=_elU[(i-1)*nDofsPerNode+j-1]*fe._shp_bulk.shape_grad(i).coeff(2);

                    _gpGradV[j-1].coeffRef(0)+=_elV[(i-1)*nDofsPerNode+j-1]*fe._shp_bulk.shape_grad(i).coeff(0);
                    _gpGradV[j-1].coeffRef(1)+=_elV[(i-1)*nDofsPerNode+j-1]*fe._shp_bulk.shape_grad(i).coeff(1);
                    _gpGradV[j-1].coeffRef(2)+=_elV[(i-1)*nDofsPerNode+j-1]*fe._shp_bulk.shape_grad(i).coeff(2);

                    _gpCoord.coeffRef(0)+=_elNodes(i,1)*fe._shp_bulk.shape_value(i);
                    _gpCoord.coeffRef(1)+=_elNodes(i,2)*fe._shp_bulk.shape_value(i);
                    _gpCoord.coeffRef(2)+=_elNodes(i,3)*fe._shp_bulk.shape_value(i);
                }
                // cout<<_gpU[j-1]<<" "
                //     <<_gpGradU[j-1](0)<<" "
                //     <<_gpGradU[j-1](1)<<" "
                //     <<_gpGradU[j-1](2)<<endl;
            }
            // cout<<"***          gpCoords="<<_gpCoord(0)<<" "<<_gpCoord(1)<<" "<<_gpCoord(2)<<endl;
            // cout<<"***          shp=";
            // for(i=1;i<=nNodes;++i){
            //     cout<<"i="<<i-1<<":";
            //     // cout<<fe._shp_bulk.shape_value(i)<<" ";
            //     cout<<fe._shp_bulk.shape_grad(i)(0)<<" "
            //         <<fe._shp_bulk.shape_grad(i)(1)<<" "
            //         <<fe._shp_bulk.shape_grad(i)(2)<<endl;
            // }

            for(j=1;j<=_nHist;++j){
                _gpHist[j-1]=0.0;
            }

            //*****************************************************
            //*** For user material calculation(UMAT)
            //*****************************************************
            mateSystem.RunMateLib(imate,mateindex,nDim,t,dt,_gpCoord,_gpU,_gpV,_gpGradU,_gpGradV,_gpHist,_gpHistOld);
            
            
            //*****************************************************
            //*** For user element calculation(UEL)
            //*****************************************************
            elmtSystem.RunElmtLib(isw,iuel,nDim,nNodes,JxW,t,dt,ctan,_gpCoord,_gpU,_gpV,_gpGradU,_gpGradV,
                                  fe._shp_bulk,
                                  mateSystem._MateValues,
                                  mateSystem._Rank2TensorMateValues,
                                  mateSystem._Rank4TensorMateValues,
                                  _gpHist,_gpHistOld,_gpProj,_localK,_localR);
            
            // cout<<"qp="<<gpInd<<":"<<endl;
            
            if(isw==9){
                // do the projection assemble(local to global)
                // Projection(nNodes,Proj,_gpProj,fe._shp_bulk,DetJac);// value is not correct(too small)
                Projection(nNodes,Proj,_gpProj,fe._shp_bulk,JxW);// still smaller
                // Projection(nNodes,Proj,_gpProj,fe._shp_bulk,w);
            }
            else if(isw==4||isw==8){
                // assemble local hist to global one(for gauss point)
                AssembleLocalHistToGlobal(e,gpInd,_gpHist,Hist);
            }
        }
        //*****************************************************
        //*** assemble the contribution from local to global
        //*****************************************************
        if(isw==3||isw==6){
            // AssembleLocalToGlobal(isw,nDofs,_elDofs,_localK,_localR,AMATRIX,RHS);
            // cout<<"print:"<<endl;
            AssembleLocalToGlobal(isw,nDofs,_elDofs,_elDofsActiveFlag,_localK,_localR,AMATRIX,RHS);
        }
    }
}