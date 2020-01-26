#include "BCs/BCSystem.h"

void BCSystem::ApplyNeumannBC(Mesh &mesh,DofHandler &dofHandler,FE &fe,
                        const int &DofIndex,const double &bcvalue,string bcname,
                        Eigen::VectorXd &RHS){
    int e,ee,i,j,iInd;
    int gpInd;
    
    for(e=1;e<=mesh.GetElmtsNumViaPhyName(bcname);++e){
        ee=mesh.GetElmtIndexViaPhyName(bcname,e);
        if(_nDim==1){
            for(i=1;i<=_nNodesPerBCElmt;++i){
                j=mesh.GetIthElmtJthConn(ee,i);
                iInd=dofHandler.GetIthNodeJthDofIndex(j,DofIndex)-1;
                RHS.coeffRef(iInd)+=bcvalue;
            }
        }
        else if(_nDim==2){
            mesh.GetIthElmtNodes(ee,_elNodes);
            // cout<<"e="<<ee<<":";
            for(gpInd=1;gpInd<=fe._qp_line.GetQpPointsNum();++gpInd){
                _xi=fe._qp_line(gpInd,1);
                fe._shp_line.Calc(_xi,_elNodes);
                _JxW=fe._shp_line.GetDetJac()*fe._qp_line(gpInd,0);
                // cout<<"gp="<<gpInd<<", xi="<<xi<<", w="<<fe._qp_line(gpInd,0)<<",J="<<fe._shp_line.GetDetJac()<<endl;
                for(i=1;i<=_nNodesPerBCElmt;++i){
                    j=mesh.GetIthElmtJthConn(ee,i);
                    iInd=dofHandler.GetIthNodeJthDofIndex(j,DofIndex)-1;
                    RHS.coeffRef(iInd)+=fe._shp_line.shape_value(i)*bcvalue*_JxW;
                    // cout<<"\tiInd="<<iInd<<",x="<<_elNodes(i,1)<<",y="<<_elNodes(i,2)<<",";
                    // cout<<"j="<<j<<", gp value="<<fe._shp_line.shape_value(i)*bcvalue*JxW<<",";
                    // cout<<fe._shp_line.shape_value(i)<<endl;
                }
            }
        }
        else if(_nDim==3){
            mesh.GetIthElmtNodes(ee,_elNodes);
            for(gpInd=1;gpInd<=fe._qp_surface.GetQpPointsNum();++gpInd){
                _xi=fe._qp_surface(gpInd,1);
                _eta=fe._qp_surface(gpInd,2);
                fe._shp_surface.Calc(_xi,_eta,_elNodes);
                _JxW=fe._shp_surface.GetDetJac()*fe._qp_surface(gpInd,0);
                for(i=1;i<=_nNodesPerBCElmt;++i){
                    j=mesh.GetIthElmtJthConn(ee,i);
                    iInd=dofHandler.GetIthNodeJthDofIndex(j,DofIndex)-1;
                    RHS.coeffRef(iInd)+=fe._shp_surface.shape_value(i)*bcvalue*_JxW;
                }
            }
        }
    }
}