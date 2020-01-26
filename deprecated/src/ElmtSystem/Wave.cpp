#include "ElmtSystem/ElmtSystem.h"

void ElmtSystem::Wave(const int &isw,const int &nDim,const int &nNodes,const double &JxW,
                 const double &t,const double &dt,const double (&ctan)[2],
                 const Eigen::Vector3d &gpCoord,
                 const vector<double> &gpU,const vector<double> &gpV,
                 const vector<Eigen::Vector3d> &gpGradU,const vector<Eigen::Vector3d> &gpGradV,
                 const ShapeFun &shp,
                 const vector<double> &MaterialValues,
                 const vector<RankTwoTensor> &Rank2MaterialValues,
                 const vector<RankFourTensor> &Rank4MaterialValues,
                 vector<double> &Hist,const vector<double> &HistOld,vector<double> &Proj,
                 Eigen::MatrixXd &K,Eigen::VectorXd &rhs){

    // just to get rid of complain from compiler(since for developer, warning is treated as errors!)
    // for normal users, you don't need the following two lines,
    // for you, warnings will always be warnings
    // waring will be error only when you are a developer(we suggest developer to do this!!!)
    if(nDim){};
    if(Rank2MaterialValues.size()){};
    if(Rank4MaterialValues.size()){};
    if(gpGradV[0].coeff(0)){};

    
    // In CahnHilliard material:
    // MateValues[0]--->M (velocity of wave)
    // MateValues[1]--->dM/du
    // MateValues[2]--->F(source term of wave equation)
    // MateValues[3]--->dF/du(its derivative)
    double M=MaterialValues[0];
    double dMdu=MaterialValues[1];
    double F=MaterialValues[2];
    double dFdu=MaterialValues[3];

    if(isw==3||isw==6){
        for(int i=1;i<=nNodes;++i){
            // R_u
            rhs.coeffRef(2*(i-1)+0)+=gpV[0]*shp.shape_value(i)*JxW
                                    -gpU[1]*shp.shape_value(i)*JxW;
            // R_v
            rhs.coeffRef(2*(i-1)+1)+=gpV[1]*shp.shape_value(i)*JxW
                                    +M*(gpGradU[0]*shp.shape_grad(i))*JxW
                                    -F*shp.shape_value(i)*JxW;
            if(isw==6){
                for(int j=1;j<=nNodes;++j){
                    // Ku,udot
                    K.coeffRef(2*(i-1),2*(j-1))+=shp.shape_value(j)*shp.shape_value(i)*ctan[1]*JxW;
                    // Ku,u
                    K.coeffRef(2*(i-1),2*(j-1))+=0.0;
                    // Ku,v
                    K.coeffRef(2*(i-1),2*(j-1)+1)+=-shp.shape_value(j)*shp.shape_value(i)*JxW*ctan[0];


                    // Kv,u
                    K.coeffRef(2*(i-1)+1,2*(j-1))+=M*(shp.shape_grad(j)*shp.shape_grad(i))*ctan[0]*JxW
                                                  +dMdu*shp.shape_value(j)*(gpGradU[0]*shp.shape_grad(i))*ctan[0]*JxW
                                                  -dFdu*(shp.shape_grad(j)*shp.shape_grad(i))*ctan[0]*JxW;
                    // Kv,vdot
                    K.coeffRef(2*(i-1)+1,2*(j-1)+1)+=shp.shape_value(j)*shp.shape_value(i)*ctan[1]*JxW;
                }
            }
        }
    }
    else if(isw==4){
        // init history variable
        fill(Hist.begin(),Hist.end(),0.0);
    }
    else if(isw==8){
        // update hist
        Hist=HistOld;
    }
    else if(isw==9){
        // do projection
        Proj[0]=gpU[0];
        Proj[1]=gpU[1];
        Proj[2]=gpV[0];
        Proj[3]=gpV[1];
        Proj[4]=gpGradU[0](1);
        Proj[5]=gpGradU[0](2);
        Proj[6]=t*dt;
        Proj[7]=gpCoord(0);
        Proj[8]=gpV[0];
    }
}