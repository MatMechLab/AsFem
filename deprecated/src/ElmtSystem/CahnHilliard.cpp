#include "ElmtSystem/ElmtSystem.h"

void ElmtSystem::CahnHilliard(const int &isw,const int &nDim,const int &nNodes,const double &JxW,
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

    
    // In CahnHilliard material:
    // MateValues[0]--->M
    // MateValues[1]--->dM/dc
    // MateValues[2]--->Chi
    // MateValues[3]--->Kappa
    // MateValues[4]--->Free energy
    // MateValues[5]--->mu=dF/dc
    // MateValues[6]--->dMu/dc
    double M=MaterialValues[0];
    double dMdc=MaterialValues[1];
    // double Chi=MaterialValues[2];
    double Kappa=MaterialValues[3];
    double F=MaterialValues[4];
    double dFdc=MaterialValues[5];
    double d2Fdc2=MaterialValues[6];

    // Kappa=0.0;
    // Chi=2.5;
    // M=gpU[0]*(1-gpU[0]);
    // dMdc=1-2*gpU[0];
    // dFdc=log(gpU[0])-log(1-gpU[0])+Chi*(1-2*gpU[0]);
    // d2Fdc2=1.0/gpU[0]+1.0/(1.0-gpU[0])-Chi*2.0;
    // if(MaterialValues[1]){}
    // double D=2.0;
    // double dDdc=0.0;
    if(isw==3||isw==6){
        for(int i=1;i<=nNodes;++i){
            // R_c
            rhs.coeffRef(2*(i-1))+=gpV[0]*shp.shape_value(i)*JxW
                                  +M*(gpGradU[1]*shp.shape_grad(i))*JxW;
            // R_mu
            rhs.coeffRef(2*(i-1)+1)+=gpU[1]*shp.shape_value(i)*JxW
                                    -dFdc*shp.shape_value(i)*JxW
                                    -Kappa*(gpGradU[0]*shp.shape_grad(i))*JxW;
            if(isw==6){
                for(int j=1;j<=nNodes;++j){
                    // Kc,cdot
                    K.coeffRef(2*(i-1),2*(j-1))+=shp.shape_value(j)*shp.shape_value(i)*ctan[1]*JxW;
                    // Kc,c
                    K.coeffRef(2*(i-1),2*(j-1))+=dMdc*shp.shape_value(j)*(gpGradU[1]*shp.shape_grad(i))*ctan[0]*JxW;
                    // Kc,mu
                    K.coeffRef(2*(i-1),2*(j-1)+1)+=M*(shp.shape_grad(j)*shp.shape_grad(i))*JxW*ctan[0];

                    // Kmu,c
                    K.coeffRef(2*(i-1)+1,2*(j-1))+=-d2Fdc2*shp.shape_value(j)*shp.shape_value(i)*ctan[0]*JxW
                                                   -Kappa*(shp.shape_grad(j)*shp.shape_grad(i))*ctan[0]*JxW;
                    // Kmu,mu
                    K.coeffRef(2*(i-1)+1,2*(j-1)+1)+=shp.shape_value(j)*shp.shape_value(i)*ctan[0]*JxW;
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
        Proj[0]=F;
        Proj[1]=gpU[0];
        Proj[2]=gpU[1];
        Proj[3]=gpGradU[0](0);
        Proj[4]=gpGradU[0](1);
        Proj[5]=gpGradU[0](2);
        Proj[6]=t*dt;
        Proj[7]=gpCoord(0);
        Proj[8]=gpV[0];
        Proj[9]=gpGradV[0](0);
    }
}