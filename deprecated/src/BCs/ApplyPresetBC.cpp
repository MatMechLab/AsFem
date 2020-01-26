#include "BCs/BCSystem.h"

void BCSystem::ApplyPresetBC(Mesh &mesh,DofHandler &dofHandler,FE &fe,const double (&ctan)[2],
                        const int &DofIndex,const double &bcvalue,string bcname,
                        Eigen::SparseMatrix<double,Eigen::RowMajor> &K,
                        Eigen::VectorXd &RHS,
                        Eigen::VectorXd &U){

    int e,ee,i,ii,j,jj,iInd,jInd;
    int gpInd;
    double gpU;
    double elU[27];
    for(e=1;e<=mesh.GetElmtsNumViaPhyName(bcname);++e){
        ee=mesh.GetElmtIndexViaPhyName(bcname,e);
        // get U on each nodal point
        if(_nDim==1){
            j=mesh.GetIthElmtJthConn(ee,1);
            iInd=dofHandler.GetIthNodeJthDofIndex(j,DofIndex)-1;
            RHS.coeffRef(iInd)+=(U.coeff(iInd)-bcvalue)*_MaxKMatrixValue;
            K.coeffRef(iInd,iInd)+=-1.0*ctan[0]*_MaxKMatrixValue;
        }
        else
        {
            // get u on each nodal point
        for(i=1;i<=_nNodesPerBCElmt;++i){
            j=mesh.GetIthElmtJthConn(ee,i);
            iInd=dofHandler.GetIthNodeJthDofIndex(j,DofIndex)-1;
            elU[i-1]=U.coeff(iInd);
        }
        mesh.GetIthElmtNodes(ee,_elNodes);
        if(_nDim==2){
            // cout<<"e="<<ee<<":";
            for(gpInd=1;gpInd<=fe._qp_line.GetQpPointsNum();++gpInd){
                _xi=fe._qp_line(gpInd,1);
                fe._shp_line.Calc(_xi,_elNodes);
                _JxW=fe._shp_line.GetDetJac()*fe._qp_line(gpInd,0);
                // cout<<"gp="<<gpInd<<", xi="<<xi<<", w="<<fe._qp_line(gpInd,0)<<",J="<<fe._shp_line.GetDetJac()<<endl;
                gpU=0.0;
                for(i=1;i<=_nNodesPerBCElmt;++i){
                   gpU+=elU[i-1]*fe._shp_line.shape_value(i);
                }
                // now we apply the dirichlet bc
                for(i=1;i<=_nNodesPerBCElmt;++i){
                    ii=mesh.GetIthElmtJthConn(ee,i);
                    iInd=dofHandler.GetIthNodeJthDofIndex(ii,DofIndex)-1;
                    RHS.coeffRef(iInd)+=fe._shp_line.shape_value(i)*(gpU-bcvalue)*_MaxKMatrixValue*_JxW;
                    // cout<<"\tiInd="<<iInd<<",x="<<_elNodes(i,1)<<",y="<<_elNodes(i,2)<<",";
                    // cout<<"j="<<j<<", gp value="<<fe._shp_line.shape_value(i)*bcvalue*JxW<<",";
                    // cout<<fe._shp_line.shape_value(i)<<endl;
                    for(j=1;j<=_nNodesPerBCElmt;++j){
                        jj=mesh.GetIthElmtJthConn(ee,j);
                        jInd=dofHandler.GetIthNodeJthDofIndex(jj,DofIndex)-1;
                        K.coeffRef(iInd,jInd)+=-fe._shp_line.shape_value(j)*fe._shp_line.shape_value(i)*_MaxKMatrixValue*_JxW*ctan[0];
                    }
                }
            }
        }
        else if(_nDim==3){
            // cout<<"e="<<ee<<":";
            for(gpInd=1;gpInd<=fe._qp_surface.GetQpPointsNum();++gpInd){
                _xi=fe._qp_surface(gpInd,1);
                _eta=fe._qp_surface(gpInd,2);
                fe._shp_surface.Calc(_xi,_eta,_elNodes);
                _JxW=fe._shp_surface.GetDetJac()*fe._qp_surface(gpInd,0);
                // cout<<"gp="<<gpInd<<", xi="<<xi<<", w="<<fe._qp_line(gpInd,0)<<",J="<<fe._shp_line.GetDetJac()<<endl;
                gpU=0.0;
                for(i=1;i<=_nNodesPerBCElmt;++i){
                   gpU+=elU[i-1]*fe._shp_line.shape_value(i);
                }
                // now we apply the dirichlet bc
                for(i=1;i<=_nNodesPerBCElmt;++i){
                    ii=mesh.GetIthElmtJthConn(ee,i);
                    iInd=dofHandler.GetIthNodeJthDofIndex(ii,DofIndex)-1;
                    RHS.coeffRef(iInd)+=fe._shp_surface.shape_value(i)*(gpU-bcvalue)*_MaxKMatrixValue*_JxW;
                    // cout<<"\tiInd="<<iInd<<",x="<<_elNodes(i,1)<<",y="<<_elNodes(i,2)<<",";
                    // cout<<"j="<<j<<", gp value="<<fe._shp_line.shape_value(i)*bcvalue*JxW<<",";
                    // cout<<fe._shp_line.shape_value(i)<<endl;
                    for(j=1;j<=_nNodesPerBCElmt;++j){
                        jj=mesh.GetIthElmtJthConn(ee,j);
                        jInd=dofHandler.GetIthNodeJthDofIndex(jj,DofIndex)-1;
                        K.coeffRef(iInd,jInd)+=-fe._shp_surface.shape_value(j)*fe._shp_surface.shape_value(i)*_MaxKMatrixValue*_JxW*ctan[0];
                    }
                }
            }
        }
        }// for else
    }
}