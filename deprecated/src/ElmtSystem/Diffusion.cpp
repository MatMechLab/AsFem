#include "ElmtSystem/ElmtSystem.h"

void ElmtSystem::Diffusion(const int &isw,const int &nDim,const int &nNodes,const double &JxW,
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
    
    // In constpoisson material:
    // _MateValues[0]=1.0;// D
    // _MateValues[1]=0.0;// dD/dc
    double D=MaterialValues[0];
    double dDdc=MaterialValues[1];
    // if(MaterialValues[1]){}
    // double D=2.0;
    // double dDdc=0.0;
    if(isw==3||isw==6){
        for(int i=1;i<=nNodes;++i){
            rhs.coeffRef(i-1)+=gpV[0]*shp.shape_value(i)*JxW
                              +D*(gpGradU[0]*shp.shape_grad(i))*JxW;
            if(isw==6){
                for(int j=1;j<=nNodes;++j){
                    K.coeffRef(i-1,j-1)+=shp.shape_value(j)*shp.shape_value(i)*JxW*ctan[1]
                                        +dDdc*shp.shape_value(j)*(gpGradU[0]*shp.shape_grad(i))*JxW*ctan[0]
                                        +D*(shp.shape_grad(j)*shp.shape_grad(i))*JxW*ctan[0];
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
        Proj[1]=gpGradU[0](0);
        Proj[2]=gpGradU[0](1);
        Proj[3]=gpGradU[0](2);
        Proj[4]=t*dt;
        Proj[5]=gpCoord(0);
        Proj[6]=gpV[0];
        Proj[7]=gpGradV[0](0);
    }
}