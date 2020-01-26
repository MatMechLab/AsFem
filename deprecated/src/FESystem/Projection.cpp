#include "FESystem/FESystem.h"

void FESystem::Projection(const int &nNodes,Eigen::MatrixXd &Proj,const vector<double> &elProj,
                          ShapeFun &shp,const double &DetJac){
    double w;
    int i,j,k,iInd;
    // for(i=0;i<nNodes;++i){
    //     iInd=_elConn[i]-1;
    //     for(j=1;j<=nNodes;++j){
    //         w=DetJac*shp.shape_value(j);
    //         Proj.coeffRef(iInd,0)+=w;
    //         for(k=1;k<=_nProj;++k){
    //             Proj.coeffRef(iInd,k)+=w*elProj[k-1]*shp.shape_value(j);
    //         }
    //     }
    // }
    if(i){};
    for(j=1;j<=nNodes;++j){
        iInd=_elConn[j-1]-1;
        w=DetJac*shp.shape_value(j);
        Proj.coeffRef(iInd,0)+=w;
        for(k=1;k<=_nProj;++k){
            Proj.coeffRef(iInd,k)+=w*elProj[k-1];
        }
    }

}